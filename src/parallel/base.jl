export run_tea

function run_tea(path2grid::String, path2data::String, model_type::Type, solver; kwargs...)

    # Generate list of scenarios to be analysed
    paths = [[path2grid, joinpath(path2data, h)] for h in readdir(path2data)]
    # Make directory where to save results
    mkdir(joinpath(path2data, "results"))
    # Solve OPF for multiple scenarios in parallel
    _DC.pmap(x -> run_opf(x, model_type, solver; kwargs...), paths)
    _DC.rmprocs(_DC.workers())
end


function run_opf(paths::Vector{String}, model_type::Type, solver; kwargs...)

    # Create scenario folder
    path2out = mkdir(joinpath(dirname(last(paths)), "results", basename(last(paths))))
    # Instantiate results init_database
    database = init_database(first(paths), model_type, solver; kwargs...)
    # Run hybrid AC/DC OPF
    results = solve_mc_acdcopf(paths, model_type, solver; kwargs...)
    # Record OPF results and save them to CSV file
    export_results(results["solution"], database, path2out)
end