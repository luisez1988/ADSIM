#!/usr/bin/env julia
#
# Precompile execution file for PackageCompiler's create_app.
#
# A real analysis traces every method it uses into the system image (fast app start).
#
# Without arguments (create_app call): bundled case run in a temporary workspace
# instead of src/. In src/ it would write VTK files, probe series and a checkpoint to
# the working tree, and the checkpoint would move the next real run to a new stage.
# Shortened calculation file: seconds of tracing, still covering diffusion, advection,
# reaction kinetics, node probing and checkpoint writing.
#
# With arguments: plain launcher for the named project in src/data.

include(joinpath(@__DIR__, "..", "src", "ADSIM.jl"))
using .ADSIM

# Trace case. Must exist in src/data, covering as many solver paths as possible
const TRACE_PROJECT = "Base_carbonation_tests"

# Simulated time of the traced run: first output step and one checkpoint
const TRACE_SIMULATION_TIME = 2.0
const TRACE_TIME_PER_STEP = 1.0

# Probe interval of the traced run. Fixed here such that the sampling path is traced
# for any case-file settings
const TRACE_PROBE_INTERVAL = 0.5

"""
    prepare_trace_workspace(src_dir) -> String

Copy the trace case into a temporary directory and shorten its calculation file.

Returns the workspace path (case in its `data/` folder).
"""
function prepare_trace_workspace(src_dir)
    workspace = mktempdir(; prefix = "adsim_precompile_")
    data_dir = joinpath(workspace, "data")
    mkdir(data_dir)

    for ext in (".mesh", "_mat.toml")
        cp(joinpath(src_dir, "data", "$(TRACE_PROJECT)$(ext)"),
           joinpath(data_dir, "$(TRACE_PROJECT)$(ext)"))
    end

    # Shorter run, probing on (probe writer traced)
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
