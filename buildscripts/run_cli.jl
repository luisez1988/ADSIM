#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "ADSIM.jl"))
using .ADSIM

function run_cli()
    if isempty(ARGS)
        println("Usage: run_cli.jl <project_name>")
        println("Example: run_cli.jl Reaction_test")
        exit(1)
    end

    # Pass command-line arguments straight through to ADSIM.main()
    ADSIM.main()
end

run_cli()
