#______________________________________________________
# probe_reader.jl
#
# Reads the node probe time series ADSIM writes when the
# [probing] block of the calculation file names one or more
# nodes. A probe is one node sampled over time, which is the
# natural output for a 0-D case where the whole point is the
# time history rather than a spatial profile.
#______________________________________________________

module ProbeReader

export ProbeData, read_probe, probe_series

"""
    ProbeData

One node's probe time series.

# Fields
- `node_id::Int`: node the probe follows
- `columns::Vector{String}`: column names in file order
- `time::Vector{Float64}`: the `time` column
- `values::Dict{String,Vector{Float64}}`: every column keyed by name
"""
struct ProbeData
    node_id::Int
    columns::Vector{String}
    time::Vector{Float64}
    values::Dict{String,Vector{Float64}}
end

"""
    read_probe(path) -> ProbeData

Parse one `<project>_node<id>_stage<n>.csv` written by `WriteProbes`.

The file carries a `#`-prefixed preamble (project, stage, node id, coordinates,
sampling interval, units), then a single header row, then one comma-separated
row per sample written in `%.8e`.
"""
function read_probe(path::AbstractString)
    isfile(path) || error("probe file not found: $path")

    node_id = 0
    columns = String[]
    rows = Vector{Vector{Float64}}()

    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue

        if startswith(s, "#")
            m = match(r"node_id\s*=\s*(\d+)", s)
            m !== nothing && (node_id = parse(Int, m.captures[1]))
            continue
        end

        if isempty(columns)
            columns = String.(strip.(split(s, ",")))
            continue
        end

        push!(rows, [parse(Float64, f) for f in split(s, ",")])
    end

    isempty(columns) && error("no header row found in probe file: $path")

    values = Dict{String,Vector{Float64}}()
    for (j, name) in enumerate(columns)
        values[name] = [r[j] for r in rows if length(r) >= j]
    end

    haskey(values, "time") || error("probe file has no 'time' column: $path")
    return ProbeData(node_id, columns, values["time"], values)
end

"""
    probe_series(probe, field) -> (time, value)

Pull one named column out of a probe alongside its time column.

Errors listing the available columns if the field is absent, so a typo in a
case file reports what it should have said.
"""
function probe_series(probe::ProbeData, field::AbstractString)
    haskey(probe.values, field) ||
        error("Field '$field' not found in probe for node $(probe.node_id). " *
              "Available: $(probe.columns)")
    return probe.time, probe.values[field]
end

end # module
