#______________________________________________________
# ADSIM Version Module
# Centralized version management
# Version: 0.1.0
#______________________________________________________

module ADSIMVersion

"""
    get_version()

Return the version string (e.g., "0.1.0") from the VERSION file at the project root;
"0.0.0-dev" if missing.
"""
function get_version()
    version_file = joinpath(@__DIR__, "..", "VERSION")
    if isfile(version_file)
        return strip(read(version_file, String))
    else
        @warn "VERSION file not found, returning default version"
        return "0.0.0-dev"
    end
end

"""
    get_version_string()

Return the display version string, e.g. "ADSIM v0.1.0".
"""
function get_version_string()
    return "ADSIM v$(get_version())"
end

export get_version, get_version_string

end # module
