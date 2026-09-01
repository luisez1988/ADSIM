#______________________________________________________
# verify.jl
#
# Core verification driver. For one case it:
#   1. copies the ADSIM input files into src/data/
#   2. runs the solver (buildscripts/run_cli.jl)
#   3. reads the resulting .vtk snapshots directly
#   4. compares each snapshot to the analytical solution
#   5. gates PASS/FAIL on the relative L2 error vs the tolerance
#______________________________________________________

module Verify

using TOML

const VERIF_DIR = normpath(joinpath(@__DIR__, ".."))
const REPO_ROOT = normpath(joinpath(VERIF_DIR, ".."))
const SRC_DIR = joinpath(REPO_ROOT, "src")
const DATA_DIR = joinpath(SRC_DIR, "data")
const OUTPUT_DIR = joinpath(SRC_DIR, "output")

include("vtk_reader.jl")
include("probe_reader.jl")
include("analytical.jl")
include("metrics.jl")
using .VTKReader
using .ProbeReader
using .Analytical
using .Metrics

export CaseResult, run_case, list_cases

struct StepResult
    step::Int
    time::Float64
    rel_l2::Float64
    rmse::Float64
    max_abs::Float64
    rel_linf::Float64
    tv::Float64
end

struct CaseResult
    name::String
    passed::Bool
    tolerance::Float64
    worst_rel_l2::Float64
    steps::Vector{StepResult}
    message::String
    tv_tolerance::Float64
    worst_tv::Float64
end

# The error paths (solver crash, no output, bad config) have no metrics to report, so
# they keep the original six-argument form and leave the TV fields unset.
CaseResult(name, passed, tol, worst, steps, message) =
    CaseResult(name, passed, tol, worst, steps, message, Inf, NaN)

const CASES_DIR = joinpath(VERIF_DIR, "cases")

"""List available cases (subfolders of cases/ that contain a case.toml)."""
function list_cases()
    isdir(CASES_DIR) || return String[]
    return sort([d for d in readdir(CASES_DIR)
                 if isfile(joinpath(CASES_DIR, d, "case.toml"))])
end

# Copy this case's ADSIM input files into src/data so the solver can find them.
function _stage_inputs(input_dir, project)
    isdir(input_dir) || error("case folder not found: $input_dir")
    for ext in (".mesh", "_calc.toml", "_mat.toml")
        src = joinpath(input_dir, "$(project)$(ext)")
        isfile(src) || error("missing input file: $src")
        cp(src, joinpath(DATA_DIR, "$(project)$(ext)"); force=true)
    end
end

# Remove stale outputs (vtk + checkpoints + logs) for a clean run.
function _clean_outputs(project)
    isdir(OUTPUT_DIR) || return
    for f in readdir(OUTPUT_DIR)
        startswith(f, project) && rm(joinpath(OUTPUT_DIR, f); force=true)
    end
end

# Run the solver as a subprocess, exactly the way a user runs it.
function _run_solver(project; verbose=false)
    runner = joinpath(REPO_ROOT, "buildscripts", "run_cli.jl")
    cmd = Cmd(`$(Base.julia_cmd()) $runner $project`; dir=REPO_ROOT)
    if verbose
        run(cmd)
    else
        out = IOBuffer()
        run(pipeline(cmd; stdout=out, stderr=out))
    end
    return nothing
end

function _collect_vtks(project)
    files = filter(readdir(OUTPUT_DIR)) do f
        startswith(f, project * "_") && endswith(f, ".vtk")
    end
    sort!(files)
    return [joinpath(OUTPUT_DIR, f) for f in files]
end

# Probe files are named <project>_node<id>_stage<n>.csv. A case names the node it
# wants; the stage is whatever the run produced, so it is matched rather than assumed.
function _collect_probe(project, node_id)
    isdir(OUTPUT_DIR) || return nothing
    pat = Regex("^" * project * "_node" * string(node_id) * "_stage(\\d+)\\.csv\$")
    files = filter(f -> match(pat, f) !== nothing, readdir(OUTPUT_DIR))
    isempty(files) && return nothing
    sort!(files)
    return joinpath(OUTPUT_DIR, files[end])   # latest stage
end

# Score a probe time series against the analytical solution.
#
# The metrics and result structs are shared with the VTK path unchanged: they take
# two plain equal-length vectors, so a time series substitutes for a spatial profile
# without touching them. A 0-D case has no position, so x is passed as zero and the
# analytical type is expected to ignore it.
function _score_probe(name, project, node_id, field, spec, tol, tv_tol)
    path = _collect_probe(project, node_id)
    path === nothing &&
        return CaseResult(name, false, tol, Inf, StepResult[],
                          "no probe output for node $node_id (is [probing] set in the calc file?)")

    probe = read_probe(path)
    times, num_all = probe_series(probe, field)

    # Skip t = 0: the initial condition, as in the VTK path.
    keep = findall(t -> t > 0, times)
    isempty(keep) &&
        return CaseResult(name, false, tol, Inf, StepResult[], "probe has no samples beyond t=0")

    t_kept = times[keep]
    num = num_all[keep]
    exact = [Analytical.evaluate(spec, 0.0, t) for t in t_kept]

    r = relative_l2(num, exact)
    # On a 0-D case the "profile" is a time history, so total variation measures temporal
    # oscillation rather than spatial - the same defect seen along the other axis.
    tv = tv_ratio(num, exact)
    steps = [StepResult(length(t_kept), t_kept[end], r, rmse(num, exact),
                        max_abs_error(num, exact), relative_linf(num, exact), tv)]

    passed = r <= tol && tv <= tv_tol
    msg = if r > tol
        "relative L2 $(round(r*100; digits=3))% exceeds tol $(round(tol*100; digits=3))%"
    elseif tv > tv_tol
        "total variation $(round(tv; digits=3))x the exact profile exceeds tol $(tv_tol)x"
    else
        "ok"
    end
    return CaseResult(name, passed, tol, r, steps, msg, tv_tol, tv)
end

"""
    run_case(name; verbose=false) -> CaseResult

Verify a single case identified by the config file `cases/<name>.toml`.
"""
function run_case(name::AbstractString; verbose::Bool=false)
    case_dir = joinpath(CASES_DIR, name)
    cfg = TOML.parsefile(joinpath(case_dir, "case.toml"))
    project = cfg["project"]
    field = cfg["field"]
    axis = Symbol(get(cfg, "axis", "y"))
    origin = Float64(get(cfg, "boundary_coord", 0.0))
    tol = Float64(get(cfg, "tolerance", 0.02))
    # Total-variation gate. The default admits 15% more variation than the exact profile,
    # which covers ordinary discretisation wiggle on a coarse mesh while still failing on
    # a node-to-node sawtooth, whose TV runs to several times the exact value. Cases
    # whose reference profile is genuinely near-flat should raise it in case.toml, since
    # the ratio is then a small number divided by a smaller one.
    tv_tol = Float64(get(cfg, "tv_tolerance", 1.15))
    check_steps = Int.(get(cfg, "check_steps", Int[]))
    component = Symbol(get(cfg, "component", "magnitude"))  # for vector fields
    # Centre of the radial coordinate; used only when axis = "radius".
    center = (Float64(get(cfg, "center_x", 0.0)), Float64(get(cfg, "center_y", 0.0)))
    snapshot = get(cfg, "snapshot", "all")                  # "all" or "last" (steady state)
    source = get(cfg, "source", "vtk")                      # "vtk" or "probe" (node time series)
    probe_node = Int(get(cfg, "probe_node", 0))
    spec = cfg["analytical"]

    try
        _stage_inputs(case_dir, project)
        _clean_outputs(project)
        _run_solver(project; verbose=verbose)
    catch err
        return CaseResult(name, false, tol, Inf, StepResult[], "solver error: $(err)")
    end

    # A 0-D case is a time history at one node, not a profile at each instant, so it
    # is scored off the probe CSV instead of the VTK series.
    if source == "probe"
        probe_node > 0 ||
            return CaseResult(name, false, tol, Inf, StepResult[],
                              "source = \"probe\" requires probe_node in case.toml")
        return _score_probe(name, project, probe_node, field, spec, tol, tv_tol)
    elseif source != "vtk"
        return CaseResult(name, false, tol, Inf, StepResult[],
                          "unknown source '$source' (use \"vtk\" or \"probe\")")
    end

    vtks = _collect_vtks(project)
    isempty(vtks) && return CaseResult(name, false, tol, Inf, StepResult[], "no VTK output produced")

    # "last" compares only the final (steady-state) snapshot; "all" every step.
    snapshot == "last" && (vtks = vtks[end:end])

    steps = StepResult[]
    worst = 0.0
    worst_tv = 0.0
    for path in vtks
        vtk = read_vtk(path)
        # Skip t = 0: the initial condition, and singular for advection solutions.
        vtk.time <= 0 && continue
        !isempty(check_steps) && !(vtk.time_step in check_steps) && continue

        pos, num = line_profile(vtk, field; axis=axis, origin=origin,
                                component=component, center=center)
        exact = [Analytical.evaluate(spec, x, vtk.time) for x in pos]

        # Total variation is only meaningful on an ordered profile; sort by position so
        # the metric does not depend on the order line_profile happens to return.
        ord = sortperm(pos)
        tv = tv_ratio(num[ord], exact[ord])

        r = relative_l2(num, exact)
        push!(steps, StepResult(vtk.time_step, vtk.time, r, rmse(num, exact),
                                max_abs_error(num, exact), relative_linf(num, exact), tv))
        worst = max(worst, r)
        worst_tv = max(worst_tv, tv)
    end

    isempty(steps) && return CaseResult(name, false, tol, Inf, steps, "no comparable snapshots (all t=0?)")

    passed = worst <= tol && worst_tv <= tv_tol
    msg = if worst > tol
        "worst relative L2 $(round(worst*100; digits=2))% exceeds tol $(round(tol*100; digits=2))%"
    elseif worst_tv > tv_tol
        "worst total variation $(round(worst_tv; digits=3))x the exact profile exceeds tol $(tv_tol)x " *
        "- the profile oscillates node to node even though its L2 error is within tolerance"
    else
        "ok"
    end
    return CaseResult(name, passed, tol, worst, steps, msg, tv_tol, worst_tv)
end

end # module
