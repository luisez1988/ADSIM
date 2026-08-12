"""
    write_probes.jl

Module for writing node probe time series for ADSIM.

A probe follows a single node through the whole simulation and records every
nodal state variable at a fixed sampling interval, producing one CSV file per
node. The fields and their names match the point data of the VTK output, so a
probe series and a VTK snapshot of the same node agree column for column.

The sampling interval is the "time interval between data saved" set in GiD under
Calculation Data > Simulation probing, which reaches the solver as
`data_saving_interval`. It is independent of the VTK output cadence, so nodes can
be sampled far more finely than the mesh is dumped.

Sampling never alters the time stepping. A sample is taken at the first step at
or past each scheduled probe time and stamped with the true current time, rather
than shortening dt to land on it exactly. Snapping dt to probe times would change
the integration and therefore the results.

Note: This module requires read_mesh.jl to be included first to access the
MeshData type.
"""

module WriteProbes

using Printf
using ..ADSIMVersion: get_version

export ProbeWriter, init_node_probes, probe_sample!, close_probes!

# Rows written between flushes. Keeping the streams open and flushing in blocks
# avoids a syscall per sample, while bounding how much of the tail is lost if the
# solver dies part way through a run.
const FLUSH_EVERY = 25

"""
ProbeWriter

Handle for a set of open node probes.

# Fields
- `node_ids::Vector{Int}`: Probed node ids, in file order
- `streams::Vector{IOStream}`: One open CSV stream per probed node
- `paths::Vector{String}`: Path of each probe file, for reporting
- `interval::Float64`: Sampling interval [time unit]
- `next_time::Float64`: Time at which the next sample is due [time unit]
- `rows_since_flush::Int`: Rows written since the last flush
- `nsamples::Int`: Total rows written per node
- `ngases::Int`: Number of gas species, fixing the per-gas column count
- `ndim::Int`: Mesh dimensionality, fixing the velocity column count
"""
mutable struct ProbeWriter
    node_ids::Vector{Int}
    streams::Vector{IOStream}
    paths::Vector{String}
    interval::Float64
    next_time::Float64
    rows_since_flush::Int
    nsamples::Int
    ngases::Int
    ndim::Int
end

"""
    init_node_probes(project_name, stage, mesh, gas_names, probing_nodes,
                     interval, units, start_time, log_print) -> Union{ProbeWriter, Nothing}

Open one CSV probe file per requested node and write its header.

Returns `nothing` when probing is not active, which is the case when no nodes
were requested or the sampling interval is not positive. Callers should treat a
`nothing` writer as a no-op; [`probe_sample!`](@ref) and [`close_probes!`](@ref)
both accept it.

Node ids outside the mesh are reported and skipped rather than aborting the run,
so a stale probe list in an input file does not cost a simulation.

# Arguments
- `project_name`: Project name, used in the file name
- `stage::Int`: Calculation stage, used in the file name
- `mesh`: Mesh data structure, for the node count and coordinates
- `gas_names::Vector{String}`: Gas species names, in solver column order
- `probing_nodes`: Dictionary holding `"nodes_to_probe"` and `"number_of_nodes"`
- `interval::Float64`: Sampling interval [time unit]
- `units::Dict`: Unit labels, recorded in the file preamble
- `start_time::Float64`: Time of the first sample [time unit]
- `log_print`: Logging function

# Returns
- `Union{ProbeWriter, Nothing}`: Probe handle, or `nothing` when inactive
"""
function init_node_probes(project_name, stage::Int, mesh, gas_names::Vector{String},
                          probing_nodes, interval::Float64, units,
                          start_time::Float64, log_print)
    requested = Int[Int(n) for n in get(probing_nodes, "nodes_to_probe", Int[])]

    # The count field and the id list are written independently by the GiD
    # interface, so a mismatch means the input file is inconsistent.
    declared = Int(get(probing_nodes, "number_of_nodes", length(requested)))
    if declared != length(requested)
        log_print("Warning: [probing] declares number_of_nodes = $(declared) but lists " *
                  "$(length(requested)) node(s) in nodes_to_probe. Using the listed nodes.")
    end

    isempty(requested) && return nothing

    if interval <= 0.0
        log_print("Warning: $(length(requested)) node(s) requested for probing but " *
                  "data_saving_interval = $(interval) is not positive. No probe files " *
                  "will be written.")
        return nothing
    end

    # Drop ids the mesh cannot resolve, and duplicates, keeping the requested order
    node_ids = Int[]
    for nid in requested
        if nid < 1 || nid > mesh.num_nodes
            log_print("Warning: probe node $(nid) is outside the mesh " *
                      "(1..$(mesh.num_nodes)) and will be skipped.")
        elseif nid in node_ids
            log_print("Warning: probe node $(nid) is listed more than once; " *
                      "keeping the first occurrence.")
        else
            push!(node_ids, nid)
        end
    end

    isempty(node_ids) && return nothing

    output_dir = "output"
    isdir(output_dir) || mkdir(output_dir)

    ngases = length(gas_names)
    ndim = size(mesh.coordinates, 2)

    streams = IOStream[]
    paths = String[]
    for nid in node_ids
        path = joinpath(output_dir, "$(project_name)_node$(nid)_stage$(stage).csv")
        io = open(path, "w")
        write_probe_header(io, project_name, stage, nid, mesh, gas_names, interval,
                           units, ndim)
        push!(streams, io)
        push!(paths, path)
    end

    log_print("   ✓ Probing $(length(node_ids)) node(s) every $(interval) " *
              "$(get(units, "time_unit", "s")): $(join(node_ids, ", "))")

    return ProbeWriter(node_ids, streams, paths, interval, start_time, 0, 0, ngases, ndim)
end

"""
    write_probe_header(io, project_name, stage, node_id, mesh, gas_names, interval, units, ndim)

Write the commented preamble and the column header row of one probe file.

The preamble records what was simulated and in which units, so a probe file
remains interpretable on its own once separated from the input deck. Every line
of it starts with `#`, which the common readers skip by default.
"""
function write_probe_header(io::IO, project_name, stage::Int, node_id::Int, mesh,
                            gas_names::Vector{String}, interval::Float64, units, ndim::Int)
    coords = mesh.coordinates[node_id, :]
    coord_str = join([@sprintf("%.6g", c) for c in coords], ", ")

    println(io, "# ADSIM v$(get_version()) - node probe time series")
    println(io, "# project = $(project_name)")
    println(io, "# stage = $(stage)")
    println(io, "# node_id = $(node_id)")
    println(io, "# coordinates = ($(coord_str))")
    println(io, "# probe_interval = $(interval) $(get(units, "time_unit", "s"))")
    println(io, "# units: time=$(get(units, "time_unit", "s")), " *
                "length=$(get(units, "geometry_unit", "m")), " *
                "mass=$(get(units, "mass_unit", "kg")), " *
                "pressure=$(get(units, "pressure_unit", "Pa")), " *
                "temperature=$(get(units, "temperature_unit", "K"))")
    println(io, "# concentrations are in mol/m^3, rates per unit time")

    # Column order mirrors the VTK point data, so the two outputs line up
    cols = String["time"]
    push!(cols, "Total_Concentration")
    push!(cols, "Absolute_Pressure")
    push!(cols, "Reaction_Rate")
    push!(cols, "Lime_Concentration")
    push!(cols, "CaCO3_Concentration")
    push!(cols, "Degree_of_Carbonation")
    push!(cols, "Volumetric_Binder_Content")
    push!(cols, "Temperature")
    push!(cols, "Temperature_Rate")
    for g in gas_names
        push!(cols, "$(g)_Concentration")
    end
    for g in gas_names
        push!(cols, "$(g)_Concentration_Rate")
    end
    axis = ("X", "Y", "Z")
    for d in 1:ndim
        push!(cols, "Gas_Seepage_Velocity_$(axis[d])")
    end

    println(io, join(cols, ","))
    return nothing
end

"""
    probe_sample!(pw, time, C_g, dC_g_dt, total_concentration, P, T, dT_dt,
                  C_lime, C_caco3, dC_lime_dt, degree_of_carbonation,
                  binder_content, v; force=false)

Record one row per probed node if a sample is due at `time`.

Does nothing when `pw` is `nothing`, or when `time` has not yet reached the next
scheduled probe time and `force` is not set. Pass `force=true` to record the
initial state regardless of schedule.

The schedule advances in whole multiples of the interval rather than from the
sample time, so the series does not drift when a step overshoots a probe time. If
several probe times fall inside one step, the schedule skips past them: the
solver cannot report state it never computed.

# Arguments
- `pw::Union{ProbeWriter, Nothing}`: Probe handle
- `time::Float64`: Current simulation time [time unit]
- remaining arguments: the nodal state arrays, as passed to the VTK writer
- `force::Bool`: Record unconditionally, used for the first sample
"""
function probe_sample!(pw::Union{ProbeWriter, Nothing}, time::Float64,
                       C_g, dC_g_dt, total_concentration, P, T, dT_dt,
                       C_lime, C_caco3, dC_lime_dt, degree_of_carbonation,
                       binder_content, v; force::Bool=false)
    pw === nothing && return nothing
    (force || time >= pw.next_time) || return nothing

    for (k, nid) in enumerate(pw.node_ids)
        io = pw.streams[k]
        @printf(io, "%.8e", time)
        @printf(io, ",%.8e", total_concentration[nid])
        @printf(io, ",%.8e", P[nid])
        @printf(io, ",%.8e", dC_lime_dt[nid])
        @printf(io, ",%.8e", C_lime[nid])
        @printf(io, ",%.8e", C_caco3[nid])
        @printf(io, ",%.8e", degree_of_carbonation[nid])
        @printf(io, ",%.8e", binder_content[nid])
        @printf(io, ",%.8e", T[nid])
        @printf(io, ",%.8e", dT_dt[nid])
        for g in 1:pw.ngases
            @printf(io, ",%.8e", C_g[nid, g])
        end
        for g in 1:pw.ngases
            @printf(io, ",%.8e", dC_g_dt[nid, g])
        end
        for d in 1:pw.ndim
            @printf(io, ",%.8e", v[nid, d])
        end
        println(io)
    end

    pw.nsamples += 1
    pw.rows_since_flush += 1
    if pw.rows_since_flush >= FLUSH_EVERY
        for io in pw.streams
            flush(io)
        end
        pw.rows_since_flush = 0
    end

    # Advance the schedule past `time` in one step, so a tiny interval relative to
    # dt costs a division rather than a long loop
    if time >= pw.next_time
        missed = floor((time - pw.next_time) / pw.interval) + 1
        pw.next_time += missed * pw.interval
    end

    return nothing
end

"""
    close_probes!(pw, log_print) -> Nothing

Flush and close every open probe file, reporting what was written.

Accepts `nothing`, so callers need not test whether probing was active.
"""
function close_probes!(pw::Union{ProbeWriter, Nothing}, log_print)
    pw === nothing && return nothing

    for io in pw.streams
        flush(io)
        close(io)
    end

    log_print("   ✓ Wrote $(length(pw.paths)) node probe file(s), " *
              "$(pw.nsamples) sample(s) each:")
    for path in pw.paths
        log_print("      $(path)")
    end

    return nothing
end

end # module WriteProbes
