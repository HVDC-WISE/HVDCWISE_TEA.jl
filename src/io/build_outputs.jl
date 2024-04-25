function build_outputs(case_name::String, raw_results::Dict)
    results_structure = Dict(
        "bus" => Dict("va" => "°?", "vm" => "pu"),  # FIXME verify the angle units and convert magnitude into kV
        "branch" => Dict("pf" => "MW", "pt" => "MW", "qf" => "MVAr", "qt" => "MVAr"),
        "gen" => Dict("pg" => "MW", "pgcurt" => "MW", "qg" => "MVAr"),
        "load" => Dict("ered" => "MWh", "eshift_down" => "MWh", "eshift_up" => "MWh", "pflex" => "MW", "pred" => "MW", "pshift_down" => "MW", "pshift_up" => "MW"),
        "storage" => Dict("ps" => "MW", "qs" => "MVAr", "qsc" => "MVAr", "sc" => "MW", "sd" => "MW", "se" => "MWh")
    )

    work_dir = joinpath(_HWTEA_dir, "test\\data\\$case_name")
    excel_path = joinpath(work_dir, "$case_name"*"_results.xlsx")
    nw = raw_results["solution"]["nw"]
    hours = [h for h in keys(nw)]
    @assert length(hours) > 0 "The length of the raw results is 0"
    h1 = hours[1]  # h1 should be "1"
    base_mva = nw[h1]["baseMVA"]

    full_grid = Dict{String, Any}()  # grid model with all components, including the ones which are sometimes in default
    for h in hours
        # Loop on all hours is needed because components with status=0 are removed from the results
        for component_type in keys(nw[h])
            if typeof(nw[h][component_type]) == Dict{String,Any}
                if !haskey(full_grid, component_type)
                    full_grid[component_type] = Vector{Any}()
                end
                for id in keys(nw[h][component_type])
                    push!(full_grid[component_type], id)
                end
            end
        end
    end

    XLSX.openxlsx(excel_path, mode="w") do xf
        for component_type in keys(full_grid)
            comp_results = Dict{String, Any}()
            comp_ids = full_grid[component_type]
            @assert length(comp_ids) > 0 "There is no component of type $component_type in the raw results"
            id1 = comp_ids[1]  # id1 should be "1"
            comp_attributes = keys(nw[h1][component_type][id1])
            for attribute in comp_attributes
                if attribute in keys(results_structure[component_type])
                    if occursin("MW", results_structure[component_type][attribute])
                        base_value = base_mva
                    else
                        base_value = 1
                    end
                    for id in comp_ids
                        name = "$id-$attribute"
                        comp_results[name] = Vector{Any}()
                        for h in hours
                            if component_type in keys(nw[h]) && id in keys(nw[h][component_type])
                                push!(comp_results[name], nw[h][component_type][id][attribute] * base_value)
                            else
                                push!(comp_results[name], nothing)
                            end
                        end
                    end
                end
            end
            comp_table = DataFrame(comp_results)  # dictionary converted into DataFrame
            if "Sheet1" in XLSX.sheetnames(xf)
                comp_sheet = xf["Sheet1"]
                XLSX.rename!(comp_sheet, component_type)
            else
                XLSX.addsheet!(xf, component_type)
                comp_sheet = xf[component_type]
            end
            nrows, ncols = size(comp_table)
            comp_sheet[2:3,1] = ["Time", "h"]
            comp_sheet[4:nrows+3, 1] = [i for i in 1:nrows]
            comp_sheet[1:3,2] = ["Id", "Attribute", "Unit"]
            for (col, name) in enumerate(names(comp_table))
                id, attribute = split(name, "-")
                unit = results_structure[component_type][attribute]
                comp_sheet[1:3,col+2] = [parse(Int, id), attribute, unit]
            end
            for row in 1:nrows
                for col in 1:ncols
                    comp_sheet[row+3, col+2] = comp_table[row, col]
                end
            end
        end
    end
end