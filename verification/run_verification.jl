#!/usr/bin/env julia
#______________________________________________________
# run_verification.jl
#
# Run ADSIM's verification suite: execute each case's
# solver run, compare against the analytical solution, and
# report PASS/FAIL on the relative L2 error.
#
# Usage:
#   julia verification/run_verification.jl              # run all cases
#   julia verification/run_verification.jl diffusion_1d # run selected cases
#   julia verification/run_verification.jl --verbose    # show solver output
#   julia verification/run_verification.jl --list       # list available cases
#
# Exit code is 0 when every case passes, 1 otherwise
# (so it can be wired into CI or a pre-commit check).
#______________________________________________________

include(joinpath(@__DIR__, "src", "verify.jl"))
using .Verify
using Printf

function main()
    args = copy(ARGS)
    verbose = false
    if "--verbose" in args
        verbose = true
        filter!(a -> a != "--verbose", args)
    end
    if "--list" in args
        println("Available verification cases:")
        for c in list_cases()
            println("  - ", c)
        end
        return 0
    end

    cases = isempty(args) ? list_cases() : args
    if isempty(cases)
        println("No verification cases found in verification/cases/.")
        return 1
    end

    println("="^70)
    println("ADSIM verification suite  (",length(cases)," case(s))")
    println("="^70)

    results = Verify.CaseResult[]
    for name in cases
        @printf("\n▶ %-16s running solver ...\n", name)
        r = run_case(name; verbose=verbose)
        push!(results, r)

        status = r.passed ? "PASS ✅" : "FAIL ❌"
        @printf("  %-16s %s\n", name, status)
        if isempty(r.steps)
            @printf("      %s\n", r.message)
        else
            @printf("      worst rel-L2 = %.3f%%   (tolerance %.3f%%)\n",
                    r.worst_rel_l2 * 100, r.tolerance * 100)
            @printf("      %-6s %-12s %-10s %-10s %-10s\n", "step", "time[s]", "relL2[%]", "rmse", "maxAbs")
            for s in r.steps
                @printf("      %-6d %-12.1f %-10.3f %-10.4g %-10.4g\n",
                        s.step, s.time, s.rel_l2 * 100, s.rmse, s.max_abs)
            end
        end
    end

    println("\n", "="^70)
    println("Summary")
    println("="^70)
    npass = count(r -> r.passed, results)
    for r in results
        @printf("  %-16s %s   worst rel-L2 %.3f%% / tol %.3f%%\n",
                r.name, r.passed ? "PASS" : "FAIL",
                isfinite(r.worst_rel_l2) ? r.worst_rel_l2 * 100 : NaN, r.tolerance * 100)
    end
    @printf("\n%d / %d cases passed\n", npass, length(results))

    return npass == length(results) ? 0 : 1
end

exit(main())
