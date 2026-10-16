"""
    write_probes.jl

Node probe time series output for ADSIM.

One CSV file per probed node with every nodal state variable at a fixed sampling
interval. Field names and order match the VTK point data.

Sampling interval: "time interval between data saved" in GiD (Calculation Data >
Simulation probing), passed to the solver as `data_saving_interval`. Independent of
the VTK output interval.

Sampling does not modify the time stepping. A sample is taken at the first step at or
after each scheduled probe time, with the actual current time. Snapping dt to probe
times would change the integration and the results.

Requires read_mesh.jl to be included first (MeshData type).
"""

module WriteProbes

using Printf
using ..ADSIMVersion: get_version

export ProbeWriter, init_node_probes, probe_sample!, close_probes!

# Rows between flushes. Block flushing avoids a syscall per sample and bounds the
# data lost if the solver stops mid-run.
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
- `ngases::Int`: Number of gas species (per-gas column count)
- `ndim::Int`: Mesh dimensionality (velocity column count)
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

Returns `nothing` when probing is inactive (no nodes requested or non-positive
interval). [`probe_sample!`](@ref) and [`close_probes!`](@ref) accept `nothing` as a no-op.

Node ids outside the mesh are reported and skipped instead of aborting the run.

# Arguments
- `project_name`, `stage`: Used in the file name
- `gas_names`: Gas species names, solver column order
- `probing_nodes`: Dictionary with `"nodes_to_probe"` and `"number_of_nodes"`
- `interval`: Sampling interval [time unit]
- `units`: Unit labels for the file preamble
- `start_time`: Time of the first sample [time unit]
"""
function init_node_probes(project_name, stage::Int, mesh, gas_names::Vector{String},
                          probing_nodes, interval::Float64, units,
                          start_time::Float64, log_print)
    requested = Int[Int(n) for n in get(probing_nodes, "nodes_to_probe", Int[])]

    # Count and id list written independently by GiD; mismatch = inconsistent input file
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

    # Drop ids outside the mesh and duplicates, requested order kept
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

Preamble with project, node and units, such that the file is self-describing. Lines
start with `#` (skipped by common CSV readers).
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

    # Same column order as the VTK point data
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

No-op when `pw` is `nothing`, or when `time` is before the next probe time without
`force`. `force=true` records the initial state.

Schedule advances in whole multiples of the interval (no drift when a step overshoots a
probe time). Multiple probe times within one step: skipped, no interpolation.

# Arguments
- `time`: Current simulation time [time unit]
- Remaining positional arguments: nodal state arrays, as passed to the VTK writer
- `force`: Record unconditionally (first sample)
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

    # Advance past `time` with one division instead of a loop (small interval vs. dt)
    if time >= pw.next_time
        missed = floor((time - pw.next_time) / pw.interval) + 1
        pw.next_time += missed * pw.interval
    end

    return nothing
end

"""
    close_probes!(pw, log_print) -> Nothing

Flush and close every open probe file and log the written files. No-op for `nothing`.
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
