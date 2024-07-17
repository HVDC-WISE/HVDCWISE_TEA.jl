export run_tea


function run_opf(paths::Vector{Any}, hours::Int, database::Dict{String, <:Any}, model_type::Type, solver; kwargs...)

    # Derive micro-scenario name and create scenario results folder
    scenario = scenario_name(paths)
    # Instantiate power system model
    network = parse_data(paths...)
    energy = Dict{String, Any}()
    # Run hybrid AC/DC OPF
    for nw_ids in Iterators.partition(_FP.nw_ids(network), hours)

        slice = slice_multinetwork(network, nw_ids)
        set_storage_energy_initial!(slice, energy)
        solution = solve_mc_acdcopf(slice, model_type, solver; kwargs...)

        if is_feasible(solution)
            energy = get_storage_energy_final(solution["solution"])
            export_results(solution["solution"], deepcopy(database), scenario)
        else
            msg = "OPF between hours $(first(nw_ids)) and $(last(nw_ids)) did not converged. Exiting iteration for scenario $(basename(scenario))."
            Memento.error(_LOGGER, msg)
            break
        end
    end
end


function run_tea(path2grid::String, path2data::String, hours::Int, model_type::Type, solver; kwargs...)

    # Generate list of paths to scenarios to be analysed
    scenarios = scenario_product(path2grid, path2data)
    # Instantiate results database
    database = init_database(first(scenarios), model_type, solver; kwargs...)
    # Solve OPF for multiple scenarios in parallel
    _DC.pmap(x -> run_opf(x, hours, database, model_type, solver; kwargs...), scenarios)
    _DC.rmprocs(_DC.workers())
end