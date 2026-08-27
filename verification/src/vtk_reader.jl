#______________________________________________________
# vtk_reader.jl
#
# Minimal reader for the ASCII legacy VTK files that ADSIM
# writes (see src/write_vtk.jl). This is what replaces the
# manual "Plot Over Line -> export CSV" step in ParaView:
# every field ADSIM saves is already in the .vtk file, so we
# parse it directly and build the 1D profile ourselves.
#______________________________________________________

module VTKReader

export VTKData, read_vtk, line_profile

"""
    VTKData

Parsed contents of one ADSIM `.vtk` snapshot.

- `time`        : physical time read from the header (`Time = ...`)
- `time_step`   : step index read from the header
- `points`      : Nnodes x 3 matrix of node coordinates
- `scalars`     : Dict field_name => Vector (one value per node)
- `vectors`     : Dict field_name => Nnodes x 3 matrix
"""
struct VTKData
    time::Float64
    time_step::Int
    points::Matrix{Float64}
    scalars::Dict{String,Vector{Float64}}
    vectors::Dict{String,Matrix{Float64}}
end

# Pull "Time = 300.0" and "Time Step 1" out of the header comment line.
function _parse_header(line::AbstractString)
    t = 0.0
    step = 0
    m = match(r"Time\s*=\s*([-+0-9.eE]+)", line)
    m !== nothing && (t = parse(Float64, m.captures[1]))
    ms = match(r"Time Step\s*(\d+)", line)
    ms !== nothing && (step = parse(Int, ms.captures[1]))
    return t, step
end

"""
    read_vtk(path) -> VTKData

Read an ADSIM legacy-ASCII UNSTRUCTURED_GRID `.vtk` file.
Only POINTS, SCALARS and VECTORS blocks are needed for verification;
CELLS / CELL_TYPES are skipped.
"""
function read_vtk(path::AbstractString)
    lines = readlines(path)
    time, step = length(lines) >= 2 ? _parse_header(lines[2]) : (0.0, 0)

    npoints = 0
    points = zeros(Float64, 0, 3)
    scalars = Dict{String,Vector{Float64}}()
    vectors = Dict{String,Matrix{Float64}}()

    i = 1
    n = length(lines)
    while i <= n
        line = strip(lines[i])
        toks = split(line)

        if !isempty(toks) && toks[1] == "POINTS"
            npoints = parse(Int, toks[2])
            points = zeros(Float64, npoints, 3)
            for k in 1:npoints
                i += 1
                v = split(strip(lines[i]))
                points[k, 1] = parse(Float64, v[1])
                points[k, 2] = parse(Float64, v[2])
                points[k, 3] = parse(Float64, v[3])
            end

        elseif !isempty(toks) && toks[1] == "SCALARS"
            name = toks[2]
            i += 1  # LOOKUP_TABLE line
            data = zeros(Float64, npoints)
            for k in 1:npoints
                i += 1
                data[k] = parse(Float64, strip(lines[i]))
            end
            scalars[name] = data

        elseif !isempty(toks) && toks[1] == "VECTORS"
            name = toks[2]
            data = zeros(Float64, npoints, 3)
            for k in 1:npoints
                i += 1
                v = split(strip(lines[i]))
                data[k, 1] = parse(Float64, v[1])
                data[k, 2] = parse(Float64, v[2])
                data[k, 3] = parse(Float64, v[3])
            end
            vectors[name] = data
        end
        i += 1
    end

    return VTKData(time, step, points, scalars, vectors)
end

# Resolve a per-node value array for `field`: a scalar directly, or a chosen
# component / magnitude of a vector field.
function _field_values(vtk::VTKData, field::AbstractString, component::Symbol)
    if haskey(vtk.scalars, field)
        return vtk.scalars[field]
    elseif haskey(vtk.vectors, field)
        v = vtk.vectors[field]
        component === :x && return v[:, 1]
        component === :y && return v[:, 2]
        component === :z && return v[:, 3]
        component === :magnitude && return vec(sqrt.(sum(abs2, v; dims=2)))
        error("Unknown vector component '$component' (use :x, :y, :z or :magnitude)")
    else
        error("Field '$field' not found in VTK. Scalars: $(collect(keys(vtk.scalars))); Vectors: $(collect(keys(vtk.vectors)))")
    end
end

"""
    line_profile(vtk, field; axis=:y, origin=0.0, component=:magnitude, center=(0.0,0.0))

Build a 1D profile of `field` along coordinate `axis` (`:x`, `:y` or `:z`),
which is exactly what a ParaView "Plot Over Line" produced. `position`
is measured as `abs(coord - origin)`, so set `origin` to the location of
the driving boundary to match an analytical solution written with x=0 at
that boundary.

Returns `(position, value)` sorted by position. Nodes that share the same
position (e.g. the two columns of a 1D strip mesh) are averaged.

`axis = :radius` instead measures each node's distance from `center` and scores
an axisymmetric solution. That path deliberately keeps every node as its own
sample rather than averaging a ring: on a mesh that is not itself axisymmetric,
the spread between nodes at equal radius IS the error such a solution must not
have, and averaging it away would hide a direction-dependent operator — which is
the whole reason a 2D case is in the suite.
"""
function line_profile(vtk::VTKData, field::AbstractString; axis::Symbol=:y,
                      origin::Float64=0.0, component::Symbol=:magnitude,
                      center::Tuple{Float64,Float64}=(0.0, 0.0))
    value = _field_values(vtk, field, component)

    if axis === :radius
        r = [hypot(vtk.points[k, 1] - center[1], vtk.points[k, 2] - center[2])
             for k in axes(vtk.points, 1)]
        order = sortperm(r)
        return r[order], value[order]
    end

    col = axis === :x ? 1 : axis === :y ? 2 : 3
    coord = vtk.points[:, col]

    # Average duplicate coordinates (round to tame float noise like 0.0999999).
    acc = Dict{Float64,Vector{Float64}}()
    for k in eachindex(coord)
        key = round(coord[k]; digits=9)
        push!(get!(acc, key, Float64[]), value[k])
    end
    keys_sorted = sort(collect(keys(acc)))
    pos = [abs(c - origin) for c in keys_sorted]
    val = [sum(acc[c]) / length(acc[c]) for c in keys_sorted]

    order = sortperm(pos)
    return pos[order], val[order]
end

end # module
