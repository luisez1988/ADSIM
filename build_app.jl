using PackageCompiler
create_app(pwd(), joinpath(pwd(), "ADSIM_app");
    precompile_execution_file="buildscripts/run_cli.jl",
    force=true)
