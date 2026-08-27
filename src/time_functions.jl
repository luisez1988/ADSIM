#------------------------------------------------------------------------------
# ADSIM Time Function Module
# Named, reusable curves f(t) that drive boundary conditions in time.
#
# A time function is defined once in the [[time_function]] section of the
# calculation TOML and referenced by integer id from the *_tf blocks of the mesh
# file. Id 0 is reserved and means "no curve": the boundary keeps the constant
# value typed in the mesh file, so every model written before this feature
# existed runs unchanged.
#
# Four kinds are supported:
#   :constant  a flat value, mostly useful as a multiplier of 1.0
#   :ramp      linear between (t_start, v_start) and (t_end, v_end), held flat outside
#   :periodic  mean + amplitude * wave(2π (t - t_start)/period + phase)
#   :table     lookup in a two-column CSV of time and value
#
# The mode decides how the curve combines with the value in the mesh file:
#   :absolute    the curve IS the boundary value in physical units, and the
#                number in the mesh file is ignored
#   :multiplier  applied BC = (mesh file value) * f(t), with f dimensionless
#
# A table is a continuous function of t, not a sequence of steps to march
# through: the lookup is a stateless binary search, so nothing here couples to
# the solver's dt, and a checkpoint restart at an arbitrary time is exact.
#------------------------------------------------------------------------------

"""
    TimeFunction

One named curve f(t). Only the fields relevant to `kind` carry meaning; the rest
hold their defaults. Kept as a single concrete struct rather than an abstract
hierarchy because a handful of unused Float64 fields costs less than a dynamic
dispatch inside the time loop.

# Fields
- `id::Int`: 1-based identifier referenced from the mesh file. Id 0 is never stored.
- `name::String`: Name shown in the GiD tree, used only for error messages
- `kind::Symbol`: `:constant`, `:ramp`, `:periodic` or `:table`
- `mode::Symbol`: `:absolute` or `:multiplier`
- `value::Float64`: Constant value (`:constant`)
- `t_start, t_end, v_start, v_end::Float64`: Ramp endpoints (`:ramp`)
- `wave::Symbol`: `:sine`, `:square`, `:triangle` or `:sawtooth` (`:periodic`)
- `mean, amplitude, period, phase::Float64`: Periodic parameters
- `times, values::Vector{Float64}`: Table samples, times strictly increasing (`:table`)
- `interp::Symbol`: `:linear` or `:step` (zero order hold) between table samples
- `before, after::Symbol`: Behaviour outside the table range, `:hold`, `:zero`,
  `:repeat` or `:error`
- `average::Bool`: Table only. Return the mean of f over the step rather than a
  point sample, so a flux keeps the right total mass when dt over-steps the file
- `file::String`: Source path of the table, for error messages
"""
struct TimeFunction
    id::Int
    name::String
    kind::Symbol
    mode::Symbol

    # :constant
    value::Float64

    # :ramp
    t_start::Float64
    t_end::Float64
    v_start::Float64
    v_end::Float64

    # :periodic
    wave::Symbol
    mean::Float64
    amplitude::Float64
    period::Float64
    phase::Float64

    # :table
    times::Vector{Float64}
    values::Vector{Float64}
    interp::Symbol
    before::Symbol
    after::Symbol
    average::Bool
    file::String
end

#------------------------------------------------------------------------------
# CSV reading
#------------------------------------------------------------------------------
"""
    read_time_series_csv(path::AbstractString) -> (times, values)

Read a two-column time series. Blank lines and lines opening with `#` are
skipped, and a leading header row is detected by its first field failing to
parse as a number, so both of these are accepted:

    # measured chamber pressure
    time,pressure
    0.0,101325.0
    3600.0,102000.0

and the same file without the header. Times must be strictly increasing; the
check is done here rather than at the lookup so a malformed file is reported
once, with its name and row, instead of silently returning wrong values every
time step.

Comma, semicolon, tab and whitespace all work as separators, which covers the
usual spreadsheet exports without pulling in a CSV dependency.
"""
function read_time_series_csv(path::AbstractString)
    isfile(path) || error("time series file not found: $(path)")

    times = Float64[]
    values = Float64[]
    header_checked = false

    for (row, line) in enumerate(eachline(path))
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue

        fields = filter(!isempty, strip.(split(s, r"[,;\t ]+")))
        length(fields) >= 2 ||
            error("time series $(path), line $(row): expected two columns, got \"$(s)\"")

        t = tryparse(Float64, fields[1])
        v = tryparse(Float64, fields[2])

        if t === nothing || v === nothing
            # The first unparseable row is taken to be the header. A later one is
            # a genuine error and must not be swallowed.
            if !header_checked && isempty(times)
                header_checked = true
                continue
            end
            error("time series $(path), line $(row): could not parse \"$(s)\" as two numbers")
        end
        header_checked = true

        if !isempty(times) && t <= times[end]
            error("time series $(path), line $(row): time $(t) does not increase past " *
                  "the previous value $(times[end]). Times must be strictly increasing.")
        end

        push!(times, t)
        push!(values, v)
    end

    isempty(times) && error("time series $(path) contains no data rows")
    return times, values
end

#------------------------------------------------------------------------------
# Reading the definitions
#------------------------------------------------------------------------------
"""
    read_time_functions(calc_data::Dict, data_dir::AbstractString) -> Dict{Int, TimeFunction}

Build the curve table from the `[[time_function]]` array of the calculation
TOML. The section is optional: a file written without it simply yields an empty
dictionary and every boundary stays constant.

`file` entries of table curves are resolved relative to `data_dir`, the folder
holding the calculation file, so a CSV travels with the input set.
"""
function read_time_functions(calc_data::Dict, data_dir::AbstractString = ".")
    functions = Dict{Int, TimeFunction}()

    entries = get(calc_data, "time_function", Any[])
    isempty(entries) && return functions

    for (position, entry) in enumerate(entries)
        id = Int(get(entry, "id", position))
        name = String(get(entry, "name", "time function $(id)"))

        id > 0 || error("time function \"$(name)\": id must be positive, id 0 is " *
                        "reserved to mean a constant boundary condition")
        haskey(functions, id) &&
            error("time function id $(id) is defined twice (\"$(name)\" and " *
                  "\"$(functions[id].name)\")")

        kind = Symbol(lowercase(String(get(entry, "type", "constant"))))
        kind in (:constant, :ramp, :periodic, :table) ||
            error("time function \"$(name)\": unknown type \"$(kind)\", expected " *
                  "constant, ramp, periodic or table")

        mode = Symbol(lowercase(String(get(entry, "mode", "absolute"))))
        mode in (:absolute, :multiplier) ||
            error("time function \"$(name)\": unknown mode \"$(mode)\", expected " *
                  "absolute or multiplier")

        # Ramp
        t_start = Float64(get(entry, "t_start", 0.0))
        t_end = Float64(get(entry, "t_end", 1.0))
        v_start = Float64(get(entry, "v_start", 0.0))
        v_end = Float64(get(entry, "v_end", 0.0))

        # Periodic
        wave = Symbol(lowercase(String(get(entry, "wave", "sine"))))
        mean = Float64(get(entry, "mean", 0.0))
        amplitude = Float64(get(entry, "amplitude", 0.0))
        period = Float64(get(entry, "period", 1.0))
        phase = Float64(get(entry, "phase", 0.0))

        # Table
        interp = Symbol(lowercase(String(get(entry, "interpolation", "linear"))))
        before = Symbol(lowercase(String(get(entry, "before", "hold"))))
        after = Symbol(lowercase(String(get(entry, "after", "hold"))))
        average = Bool(get(entry, "average", false))
        file = String(get(entry, "file", ""))

        times = Float64[]
        values = Float64[]

        if kind == :ramp && t_end <= t_start
            error("time function \"$(name)\": t_end ($(t_end)) must be greater than " *
                  "t_start ($(t_start))")
        end

        if kind == :periodic
            wave in (:sine, :square, :triangle, :sawtooth) ||
                error("time function \"$(name)\": unknown wave \"$(wave)\", expected " *
                      "sine, square, triangle or sawtooth")
            period > 0.0 ||
                error("time function \"$(name)\": period must be positive, got $(period)")
        end

        if kind == :table
            isempty(file) &&
                error("time function \"$(name)\": type is table but no file was given")
            interp in (:linear, :step) ||
                error("time function \"$(name)\": unknown interpolation \"$(interp)\", " *
                      "expected linear or step")
            for (label, edge) in (("before", before), ("after", after))
                edge in (:hold, :zero, :repeat, :error) ||
                    error("time function \"$(name)\": unknown $(label) behaviour " *
                          "\"$(edge)\", expected hold, zero, repeat or error")
            end

            path = isabspath(file) ? file : joinpath(data_dir, file)
            times, values = read_time_series_csv(path)

            # Shift and rescale the file's time column into model time. This is how
            # a CSV recorded in hours drives a model running in seconds.
            time_scale = Float64(get(entry, "time_scale", 1.0))
            time_offset = Float64(get(entry, "time_offset", 0.0))
            time_scale > 0.0 ||
                error("time function \"$(name)\": time_scale must be positive, got $(time_scale)")
            if time_scale != 1.0 || time_offset != 0.0
                times = times .* time_scale .+ time_offset
            end

            if (before == :repeat || after == :repeat) && length(times) < 2
                error("time function \"$(name)\": repeat needs at least two samples " *
                      "to define a period")
            end

            file = path
        end

        functions[id] = TimeFunction(id, name, kind, mode,
                                     Float64(get(entry, "value", 0.0)),
                                     t_start, t_end, v_start, v_end,
                                     wave, mean, amplitude, period, phase,
                                     times, values, interp, before, after, average, file)
    end

    return functions
end

#------------------------------------------------------------------------------
# Evaluation
#------------------------------------------------------------------------------
"""
    wave_value(wave::Symbol, θ::Float64) -> Float64

Unit-amplitude, zero-mean wave form of phase angle θ in radians. All four forms
peak at +1 and trough at -1 so `amplitude` means the same thing whichever is
chosen.
"""
@inline function wave_value(wave::Symbol, θ::Float64)
    if wave == :sine
        return sin(θ)
    end

    # Fraction of the way through the cycle, in [0, 1)
    u = mod(θ / (2π), 1.0)

    if wave == :square
        return u < 0.5 ? 1.0 : -1.0
    elseif wave == :triangle
        # Rises over the first half cycle, falls over the second
        return u < 0.5 ? (4.0 * u - 1.0) : (3.0 - 4.0 * u)
    elseif wave == :sawtooth
        return 2.0 * u - 1.0
    end

    return 0.0
end

"""
    table_value(tf::TimeFunction, t::Float64) -> Float64

Point sample of a table curve.

The bracketing sample is found with `searchsortedlast`, which is O(log n) and
carries no state between calls. That matters twice over: the solver shrinks dt
to land exactly on output times, so evaluation times are not evenly spaced, and
a checkpoint restart resumes at an arbitrary time with no history to replay.
"""
function table_value(tf::TimeFunction, t::Float64)
    times = tf.times
    values = tf.values
    n = length(times)
    n == 1 && return values[1]

    t_first = times[1]
    t_last = times[n]

    # Outside the tabulated range
    if t < t_first
        tf.before == :zero && return 0.0
        tf.before == :error &&
            error("time function \"$(tf.name)\": time $(t) is before the first sample " *
                  "$(t_first) in $(tf.file)")
        if tf.before == :repeat
            t = t_first + mod(t - t_first, t_last - t_first)
        else
            return values[1]        # :hold
        end
    elseif t > t_last
        tf.after == :zero && return 0.0
        tf.after == :error &&
            error("time function \"$(tf.name)\": time $(t) is past the last sample " *
                  "$(t_last) in $(tf.file)")
        if tf.after == :repeat
            t = t_first + mod(t - t_first, t_last - t_first)
        else
            return values[n]        # :hold
        end
    end

    i = searchsortedlast(times, t)
    i <= 0 && return values[1]
    i >= n && return values[n]

    tf.interp == :step && return values[i]

    dt_seg = times[i + 1] - times[i]
    dt_seg <= 0.0 && return values[i]
    θ = (t - times[i]) / dt_seg
    return values[i] * (1.0 - θ) + values[i + 1] * θ
end

"""
    table_mean(tf::TimeFunction, t0::Float64, t1::Float64) -> Float64

Mean of a table curve over `[t0, t1]`, by exact integration of the piecewise
linear (or piecewise constant, for `:step`) interpolant.

Point sampling a *flux* loses mass whenever the step over-steps the file's
sample spacing: the samples in between are never visited. Averaging over the
step makes the injected mass correct regardless of how the two spacings relate.
Only meaningful for a Neumann boundary; averaging a prescribed value has no
such justification.
"""
function table_mean(tf::TimeFunction, t0::Float64, t1::Float64)
    span = t1 - t0
    span <= 0.0 && return table_value(tf, t0)

    # Integrate by walking the sample points that fall strictly inside the step,
    # so every segment is integrated on the interval where its own linear form
    # holds. Endpoints are evaluated through table_value so the out-of-range and
    # repeat rules apply here exactly as they do to a point sample.
    integral = 0.0
    t_left = t0
    v_left = table_value(tf, t0)

    i_first = searchsortedfirst(tf.times, t0)
    i_last = searchsortedlast(tf.times, t1)

    for i in i_first:i_last
        t_node = tf.times[i]
        (t_node <= t_left || t_node >= t1) && continue
        v_node = table_value(tf, t_node)
        integral += tf.interp == :step ? v_left * (t_node - t_left) :
                                         0.5 * (v_left + v_node) * (t_node - t_left)
        t_left = t_node
        v_left = v_node
    end

    v_right = table_value(tf, t1)
    integral += tf.interp == :step ? v_left * (t1 - t_left) :
                                     0.5 * (v_left + v_right) * (t1 - t_left)

    return integral / span
end

"""
    evaluate(tf::TimeFunction, t::Float64, dt::Float64 = 0.0) -> Float64

The raw curve value at time `t`, before the mode is applied.

`dt` is only consulted by a table curve carrying `average = true`, which then
returns the mean of the curve over `[t, t + dt]` instead of the point sample.
"""
function evaluate(tf::TimeFunction, t::Float64, dt::Float64 = 0.0)
    if tf.kind == :constant
        return tf.value

    elseif tf.kind == :ramp
        t <= tf.t_start && return tf.v_start
        t >= tf.t_end && return tf.v_end
        θ = (t - tf.t_start) / (tf.t_end - tf.t_start)
        return tf.v_start * (1.0 - θ) + tf.v_end * θ

    elseif tf.kind == :periodic
        θ = 2π * (t - tf.t_start) / tf.period + tf.phase
        return tf.mean + tf.amplitude * wave_value(tf.wave, θ)

    elseif tf.kind == :table
        return (tf.average && dt > 0.0) ? table_mean(tf, t, t + dt) :
                                          table_value(tf, t)
    end

    return 0.0
end

"""
    bc_value(time_functions, tf_id::Int, base_value::Float64, t::Float64, dt::Float64 = 0.0)

The boundary value to impose at time `t`. This is the single entry point every
call site in the solver uses.

Curve id 0 means no curve, and the constant `base_value` read from the mesh file
is returned untouched. Otherwise the curve is evaluated and combined with
`base_value` according to its mode: an `:absolute` curve replaces it, a
`:multiplier` curve scales it.
"""
@inline function bc_value(time_functions, tf_id::Int, base_value::Float64,
                          t::Float64, dt::Float64 = 0.0)
    tf_id == 0 && return base_value
    tf = time_functions[tf_id]
    f = evaluate(tf, t, dt)
    return tf.mode === :multiplier ? base_value * f : f
end

"""
    validate_time_function_ids(time_functions, mesh)

Check every curve id referenced by the mesh file against the curves actually
defined. A dangling reference means the calculation file and the mesh file came
from different runs of the problemtype, which is worth stopping for: silently
treating it as a constant would run a wrong model to completion.
"""
function validate_time_function_ids(time_functions, mesh)
    referenced = Set{Int}()

    for ids in values(mesh.concentration_bc_tf); union!(referenced, ids); end
    for ids in values(mesh.uniform_flow_bc_tf); union!(referenced, ids); end
    union!(referenced, values(mesh.absolute_pressure_tf))
    union!(referenced, values(mesh.temperature_bc_tf))

    missing_ids = sort([id for id in referenced if id != 0 && !haskey(time_functions, id)])
    isempty(missing_ids) && return nothing

    error("the mesh file references time function id(s) $(join(missing_ids, ", ")) " *
          "that are not defined in the [[time_function]] section of the calculation " *
          "file. Regenerate both input files from the same GiD project.")
end

"""
    warn_if_aliased(time_functions, dt::Float64, log_print)

Report any table curve whose samples are closer together than the time step.

The lookup itself does not care how the two spacings relate, but a step larger
than the sample spacing steps straight over samples the file went to the trouble
of recording, and the boundary the solver sees is an aliased version of the
data. That is worth saying out loud once at start-up rather than leaving the
user to discover it in the results.
"""
function warn_if_aliased(time_functions, dt::Float64, log_print)
    dt > 0.0 || return nothing

    for id in sort(collect(keys(time_functions)))
        tf = time_functions[id]
        (tf.kind == :table && length(tf.times) > 1) || continue

        spacing = minimum(diff(tf.times))
        spacing >= dt && continue

        log_print("   ⚠ Warning: time function \"$(tf.name)\" samples every " *
                  "$(spacing) but the time step is $(dt). Samples between steps " *
                  "will not be seen and the boundary condition is aliased.")
        if !tf.average
            log_print("     Set average = true on this curve to integrate it over " *
                      "the step instead (recommended for flux boundaries), or " *
                      "reduce the time step.")
        end
    end

    return nothing
end

export TimeFunction, read_time_functions, read_time_series_csv
export evaluate, bc_value, validate_time_function_ids, warn_if_aliased
