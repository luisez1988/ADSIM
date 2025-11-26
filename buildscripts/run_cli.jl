#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "ADSIM.jl"))
using .ADSIM

function run_cli()
    if isempty(ARGS)
        # When PackageCompiler runs this for precompilation without arguments,
        # use a default test case to trace compilation paths
        push!(ARGS, "Reaction_test")
    end

    # Change to src/ directory where data files are located
    src_dir = joinpath(@__DIR__, "..", "src")
    cd(src_dir)

    # Pass command-line arguments straight through to ADSIM.main()
    ADSIM.main()
end

run_cli()
