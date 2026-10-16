module ADSIM

# Load the main driver and supporting routines defined in kernel.jl.
#
# kernel.jl includes version.jl itself, therefore it also works when run directly as
# a script. Including it here as well would define a second ADSIMVersion module
# under the same name, leaving a conflicting binding behind ("ignoring conflicting
# import of ADSIMVersion.get_version"). On Julia 1.12 that state breaks
# precompilation: loading the cached image segfaults in jl_validate_binding_partition.
include("kernel.jl")

using .ADSIMVersion: get_version, get_version_string

"""
    julia_main()::Cint

Entry point of the compiled application.

PackageCompiler's `create_app` calls this module, not `main`, and needs the `Cint`
return type. `main` handles its own errors and exits directly on failure, so this
only has to catch whatever is not covered there (e.g., a failure before the log file is open, which
would otherwise surface as an unhandled exception with no explanation).

# Returns
- `Cint`: 0 on success, 1 if an error escaped `main`
"""
function julia_main()::Cint
    try
        main()
    catch e
        println(stderr, "\nADSIM terminated with an unhandled error:")
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        return Cint(1)
    end
    return Cint(0)
end

export main, julia_main, get_version, get_version_string

end
