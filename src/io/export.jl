
function init_database(paths::Vector{Any}, model_type::Type, solver; kwargs...)

    names = ["qg", "qf", "qt", "qpr_fr", "qtf_to", "qconv", "qgrid", "iconv", "pdcg_shunt", "phi", "vmfilt", "vmconv", "vm"]

    # Solve exemplary model
    grid = parse_data(paths..., 10)
    results = solve_mc_acdcopf(grid, model_type, solver; kwargs...)["solution"]
    # Keep only results for one timestep
    grid = grid["nw"]["1"]
    results = results["nw"]["1"]
    filter!(x -> isa(x.second, Dict) & !isempty(x.second), results)

    db = Dict{String,Any}()
    for (key, value) in results
        db[key] = Dict{String, Matrix{Float64}}()
        for var in keys(last(first(value)))
            if key != "busdc" && var in names
                continue
            end
            dim = maximum(length.(_PM.component_table(results, key, var)[:, 2]))
            db[key][var] = fill(Inf, length(grid[key]), dim)
            if !contains(key, "bus") && dim > 1
                for (row, info) in grid[key]
                    ref = first(filter(x -> contains(x, "confi"), collect(keys(info))))
                    # Find conductor or converter pole to which the value belongs to
                    if info[ref] == 1
                        if info["connect_at"] == 0
                            cols = [1, 2]
                        elseif info["connect_at"] == 1
                            cols = [1, 3]
                        elseif info["connect_at"] == 2
                            cols = [2, 3]
                        end
                    elseif info[ref] == 2
                        cols = [1, 2, 3]
                    end
                    if dim == 2
                        pop!(cols)
                    end
                    db[key][var][parse(Int, row), cols] .= 0
                end
            elseif contains(key, "bus") && dim > 1
                replace!(db[key][var], Inf=>0.0)
            end
            replace!(db[key][var], Inf=>NaN)
        end
    end
    return db
end


function export_results(solution::Dict{String, Any}, database, path::String)

    con = ["p", "n", "r"]
    nws = string.(sort!(parse.(Int, collect(keys(solution["nw"])))))

    for (comp, vars) in database
        path2comp = mkdir(joinpath(path, comp))
        for var in keys(vars)
            path2var = joinpath(path2comp, var * ".csv")
            ls = Vector{Float64}[]
            data = _PM.component_table(solution, comp, var)
            names = reduce(vcat, [string(h) .* con[axes(vars[var], 2)] for h in axes(vars[var], 1)])
            for nw in nws
                ref = copy(vars[var])
                if size(ref, 2) > 1
                    for v in eachrow(data[nw])
                        cols = findall(iszero, ref[first(v), :])
                        ref[first(v), cols] = last(v)
                    end
                    ref = reduce(vcat, eachrow(ref))
                else
                    ref[data[nw][:, 1]] = data[nw][:, 2]
                    ref = vec(ref)
                end
                push!(ls, ref)
            end
            res = Float32.(permutedims(reduce(hcat, ls)))
            # Save file with name of variable
            CSV.write(path2var, _DF.DataFrame(collect(eachcol(res)), names))
        end
    end
end
