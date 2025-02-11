
function scenario_name(paths::Vector{Any})

    grid = chop(basename(first(paths))[:], tail=2)
    scenario = [join([basename(dirname(h)), basename(h)], "-") for h in last(paths)]
    scenario = join(scenario, "_")
    return mkdir(joinpath(dirname(first(paths)), grid, scenario))
end

function scenario_product(path2grid::String, path2data::String)

    # Make directory where to save results
    grid = chop(basename(path2grid)[:], tail=2)
    mkdir(joinpath(dirname(path2grid), grid))

    pathlist = Vector{String}[]
    for folder in readdir(path2data)
        path = joinpath(path2data, folder)
        push!(pathlist, joinpath.(path, readdir(path)))
    end
    path2scenario = vec(collect(Base.product(pathlist...)))

    return [[path2grid, collect(h)] for h in path2scenario]
end