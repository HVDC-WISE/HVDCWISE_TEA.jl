export run_tea

function run_tea(path2grid::String, path2data::String, model_type::Type, solver; kwargs...)

    # Get the path to micro-scenario data
    pathlist = Vector{String}[]
    for folder in readdir(path2data)
        path = joinpath(path2data, folder)
        push!(pathlist, joinpath.(path, readdir(path)))
    end
    path2scenario = vec(collect(Base.product(pathlist...)))
    # Generate list of scenarios to be analysed
    paths = [[path2grid, collect(h)] for h in path2scenario]
    # Make directory where to save results
    mkdir(joinpath(path2data, "results"))
    # Instantiate results database
    database = init_database(first(paths), model_type, solver; kwargs...)
    # Solve OPF for multiple scenarios in parallel
    _DC.pmap(x -> run_opf(x, database, model_type, solver; kwargs...), paths)
    _DC.rmprocs(_DC.workers())
end


function run_opf(paths::Vector{Any}, database, model_type::Type, solver; kwargs...)

    # Derive micro-scenario name and create results folder
    scenario = [join([basename(dirname(h)), basename(h)], "-") for h in last(paths)]
    scenario = join(scenario, "_")
    pathbase = dirname(dirname(first(last(paths))))
    path2out = mkdir(joinpath(pathbase, "results", scenario))
    # Run hybrid AC/DC OPF
    results = solve_mc_acdcopf(paths, model_type, solver; kwargs...)
    # Record OPF results and save them to CSV file
    export_results(results["solution"], database, path2out)
end