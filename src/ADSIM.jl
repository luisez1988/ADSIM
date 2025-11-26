module ADSIM

"""
ADSIM main entrypoint. Loads `kernel.jl`, which defines `main()` and the
supporting routines for running a project from the command line.
"""
include("kernel.jl")

export main

end
