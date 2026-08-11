module ADSIM

# Load version information
include("version.jl")
using .ADSIMVersion: get_version, get_version_string

# Load the main driver and supporting routines defined in kernel.jl.
include("kernel.jl")

export main, get_version, get_version_string

end
