#!/usr/bin/env julia
#
# Precompile execution file for PackageCompiler's create_app.
#
# Running a real analysis here is what makes the compiled app start fast: every
# method the run touches is traced and baked into the system image.
#
# With no arguments (how create_app invokes this) the trace runs a bundled case
# inside a temporary workspace rather than in src/. Running it in src/ would write
# VTK files, probe series and a checkpoint into the working tree, and the
# checkpoint would silently push the next real run into a new stage. The copied
# calculation file is shortened so tracing costs seconds while still exercising
# diffusion, advection, reaction kinetics, node probing and checkpoint writing.
#
# With arguments it behaves as a plain launcher, running the named project out of
# src/data as before.

include(joinpath(@__DIR__, "..", "src", "ADSIM.jl"))
using .ADSIM

# Case used to trace compilation. Must exist in src/data and should touch as many
# solver paths as possible.
const TRACE_PROJECT = "Base_carbonation_tests"

# Simulated time for the traced run. Only needs to be long enough to reach the
# first output step and write a checkpoint.
const TRACE_SIMULATION_TIME = 2.0
const TRACE_TIME_PER_STEP = 1.0

# Probe interval for the traced run. Pinned rather than inherited so that tracing
# covers the sampling path regardless of how the case file happens to be tuned.
const TRACE_PROBE_INTERVAL = 0.5

"""
    prepare_trace_workspace(src_dir) -> String

Copy the trace case into a temporary directory and shorten its calculation file.

Returns the workspace path, which contains a `data/` folder holding the case.
"""
function prepare_trace_workspace(src_dir)
    workspace = mktempdir(; prefix = "adsim_precompile_")
    data_dir = joinpath(workspace, "data")
    mkdir(data_dir)

    for ext in (".mesh", "_mat.toml")
        cp(joinpath(src_dir, "data", "$(TRACE_PROJECT)$(ext)"),
           joinpath(data_dir, "$(TRACE_PROJECT)$(ext)"))
    end

    # Shorten the run and switch probing on, so the probe writer is traced too
    calc = read(joinpath(src_dir, "data", "$(TRACE_PROJECT)_calc.toml"), String)
    calc = replace(calc, r"(?m)^total_simulation_time = .*$" =>
                          "total_simulation_time = $(TRACE_SIMULATION_TIME)")
    calc = replace(calc, r"(?m)^time_per_step = .*$" =>
                          "time_per_step = $(TRACE_TIME_PER_STEP)")
    calc = replace(calc, r"(?m)^data_saving_interval = .*$" =>
                          "data_saving_interval = $(TRACE_PROBE_INTERVAL)")
    calc = replace(calc, r"(?m)^number_of_nodes = .*$" => "number_of_nodes = 1")
    calc = replace(calc, r"(?m)^nodes_to_probe = .*$" => "nodes_to_probe = [1]")
    write(joinpath(data_dir, "$(TRACE_PROJECT)_calc.toml"), calc)

    return workspace
end

function run_cli()
    src_dir = abspath(joinpath(@__DIR__, "..", "src"))

    if isempty(ARGS)
        # Build-time tracing: run the bundled case in a throwaway workspace
        workspace = prepare_trace_workspace(src_dir)
        push!(ARGS, TRACE_PROJECT)
        cd(workspace)
        try
            ADSIM.main()
        finally
            cd(src_dir)
            rm(workspace; recursive = true, force = true)
        end
    else
        # Plain launcher: run the requested project out of src/data
        cd(src_dir)
        ADSIM.main()
    end
end

run_cli()
