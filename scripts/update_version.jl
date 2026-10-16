#!/usr/bin/env julia
#______________________________________________________
# ADSIM Version Update Script
# Version string update across all project files
#______________________________________________________

using Pkg

"""
    update_version(new_version::String)

Update the version string (e.g., "0.1.0") in all ADSIM project files.
The script entry point reads it from VERSION when no argument is given.

# Updated Files
1. VERSION (source of truth)
2. Project.toml
3. src/ADSIM.jl header
4. src/kernel.jl header
5. Problemtype/ADSIM_2025.gid/ADSIM_2025.tcl
6. Problemtype/ADSIM_2025.gid/ADSIM_2025.xml
7-9. Problemtype/ADSIM_2025.gid/scripts/*.tcl (3 files)
"""
function update_version(new_version::String)
    # Validate version format (semantic versioning)
    if !occursin(r"^\d+\.\d+\.\d+(-[a-zA-Z0-9\.]+)?$", new_version)
        error("Invalid version format. Use semantic versioning (e.g., 0.1.0, 1.0.0-alpha)")
    end

    project_root = dirname(@__DIR__)
    println("Updating ADSIM version to: $new_version")
    println("Project root: $project_root")
    println("="^60)

    files_updated = 0
    
    # 1. Update VERSION file
    version_file = joinpath(project_root, "VERSION")
    write(version_file, new_version)
    println("✓ Updated: VERSION")
    files_updated += 1

    # 2. Update Project.toml
    project_toml = joinpath(project_root, "Project.toml")
    if isfile(project_toml)
        content = read(project_toml, String)
        content = replace(content, r"version = \"[^\"]+\"" => "version = \"$new_version\"")
        write(project_toml, content)
        println("✓ Updated: Project.toml")
        files_updated += 1
    end

    # 3. Update src/ADSIM.jl header
    adsim_jl = joinpath(project_root, "src", "ADSIM.jl")
    if isfile(adsim_jl)
        content = read(adsim_jl, String)
        content = replace(content, r"# v\d+\.\d+\.\d+(-[a-zA-Z0-9\.]+)?" => "# v$new_version")
        content = replace(content, r"# v\d+\.x\.x" => "# v$new_version")
        write(adsim_jl, content)
        println("✓ Updated: src/ADSIM.jl")
        files_updated += 1
    end

    # 4. Update src/kernel.jl header
    kernel_jl = joinpath(project_root, "src", "kernel.jl")
    if isfile(kernel_jl)
        content = read(kernel_jl, String)
        content = replace(content, r"# v\d+\.\d+\.\d+(-[a-zA-Z0-9\.]+)?" => "# v$new_version")
        content = replace(content, r"# v\d+\.x\.x" => "# v$new_version")
        write(kernel_jl, content)
        println("✓ Updated: src/kernel.jl")
        files_updated += 1
    end

    # 5. Update ADSIM_2025.tcl
    tcl_file = joinpath(project_root, "Problemtype", "ADSIM_2025.gid", "ADSIM_2025.tcl")
    if isfile(tcl_file)
        content = read(tcl_file, String)
        content = replace(content, r"set adsim_version \"ADSIM v[^\"]+\"" => "set adsim_version \"ADSIM v$new_version\"")
        write(tcl_file, content)
        println("✓ Updated: Problemtype/ADSIM_2025.gid/ADSIM_2025.tcl")
        files_updated += 1
    end

    # 6. Update ADSIM_2025.xml
    xml_file = joinpath(project_root, "Problemtype", "ADSIM_2025.gid", "ADSIM_2025.xml")
    if isfile(xml_file)
        content = read(xml_file, String)
        content = replace(content, r"<Version>[^<]+</Version>" => "<Version>$new_version</Version>")
        write(xml_file, content)
        println("✓ Updated: Problemtype/ADSIM_2025.gid/ADSIM_2025.xml")
        files_updated += 1
    end

    # 7-9. Update script files in Problemtype/ADSIM_2025.gid/scripts/
    scripts_dir = joinpath(project_root, "Problemtype", "ADSIM_2025.gid", "scripts")
    if isdir(scripts_dir)
        for script_file in ["WriteMaterialData.tcl", "WriteCalculationData.tcl", "WriteMeshFile.tcl"]
            script_path = joinpath(scripts_dir, script_file)
            if isfile(script_path)
                content = read(script_path, String)
                # Update "2025 v0.x.x" format
                content = replace(content, r"2025 v\d+\.x\.x" => "2025 v$new_version")
                content = replace(content, r"2025 v\d+\.\d+\.\d+(-[a-zA-Z0-9\.]+)?" => "2025 v$new_version")
                # Update "ADSIM_2025" in mesh file
                if script_file == "WriteMeshFile.tcl"
                    content = replace(content, "Version: ADSIM_2025" => "Version: ADSIM v$new_version")
                end
                write(script_path, content)
                println("✓ Updated: Problemtype/ADSIM_2025.gid/scripts/$script_file")
                files_updated += 1
            end
        end
    end

    println("="^60)
    println("Version update complete!")
    println("Files updated: $files_updated")
    println("\nNext steps:")
    println("1. Review changes with: git diff")
    println("2. Commit changes: git commit -am \"Bump version to v$new_version\"")
    println("3. Create tag: git tag -a v$new_version -m \"Release v$new_version\"")
    println("4. Push tag: git push origin v$new_version")
    
    return files_updated
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        # Read from VERSION file
        version_file = joinpath(dirname(@__DIR__), "VERSION")
        if isfile(version_file)
            new_version = strip(read(version_file, String))
            println("No version specified, using VERSION file: $new_version")
        else
            println("Usage: julia update_version.jl <version>")
            println("Example: julia update_version.jl 0.1.0")
            exit(1)
        end
    else
        new_version = ARGS[1]
    end
    
    update_version(new_version)
end
