import CSV
import DataFrames
import JSON
import XLSX

function build_outputs(work_dir::String, case_name::String, raw_results::Dict)
    if !isdir(work_dir)
        mkpath(work_dir)
    end
    json_path = joinpath(work_dir, "$case_name"*"_results.json")
    stringdata = JSON.json(raw_results, 4)
    open(json_path, "w") do f
        write(f, stringdata)
    end

    build_outputs_from_json(json_path)
end

function build_outputs_from_csv(results_dir::String, base_mva::Int, micro_scenarios::Bool=true)
    # If micro_scenarios: the results directory contains 1 folder per micro-scenario. Else it contains the results of 1 micro-scenario (at the root of results_dir)
    # println("results_dir: $results_dir")
    if micro_scenarios
        results_list = Vector{String}()
        for micro_scenario in readdir(results_dir)
            micro_dir = joinpath(results_dir, micro_scenario)
            #println("micro_dir: $micro_dir")
            if isdir(micro_dir)
                push!(results_list, micro_dir)
            end
        end
    else
        results_list = [results_dir]
    end
    # println("results_list: $results_list")
    # Micro-scenarios results are aggregated in the same Excel file
    excel_path = joinpath(results_dir, "Results.xlsx")
    
    XLSX.openxlsx(excel_path, mode="w") do xf  # Create Excel
        for result_folder in results_list
            micro_scenario_name = basename(result_folder)
            for component in readdir(result_folder)
                # println("component: $component")
                comp_folder = joinpath(result_folder, component)
                # println("comp_folder: $comp_folder")
                if isdir(comp_folder)
                    component = basename(comp_folder)
                    if component in XLSX.sheetnames(xf)
                        # println("Use sheet $component")
                        # println("$component in $(XLSX.sheetnames(xf))")
                        comp_sheet = xf[component]
                        @assert XLSX.get_dimension(comp_sheet).stop.row_number > 1
                    else
                        if "Sheet1" in XLSX.sheetnames(xf)  # Rename by default sheet
                            comp_sheet = xf["Sheet1"]
                            XLSX.rename!(comp_sheet, component)
                        else
                            #println("Add sheet $component")
                            XLSX.addsheet!(xf, component)  # Add sheet
                            comp_sheet = xf[component]
                        end
                        comp_sheet[1, 1:3] = ["Time", "Micro-scenario", "Attribute"]
                        # println("comp_sheet: $comp_sheet")
                    end
                    for file_name in readdir(comp_folder)
                        if occursin(".csv", file_name)
                            # println("attribute: $attribute")
                            attribute_file = joinpath(comp_folder, file_name)
                            # println("attribute_file: $attribute_file")
                            attribute = file_name[1:length(file_name)-4]
                            # Read csv
                            csv_data = CSV.File(attribute_file, delim=',') |> DataFrames.DataFrame
                            n_hours = first(size(csv_data))
                            n_rows = XLSX.get_dimension(comp_sheet).stop.row_number

                            if n_rows in [0, 1]
                                # println(names(csv_data))
                                for (col, component_id) in enumerate(names(csv_data))
                                    comp_sheet[1, col+3] = component_id
                                end
                                n_rows = 1
                            end

                            if first(attribute) == 'p'
                                # pu powers must be converted into MW
                                base_value = base_mva
                            else
                                base_value = 1
                            end
                            for component_id in names(csv_data)
                                #=println("component_id: $component_id")
                                println("comp_sheet: $comp_sheet")
                                println("comp_sheet[1,:]: $(comp_sheet[1,:])")=#
                                @assert component_id in comp_sheet[1,:]  "component_id $component_id is not in the Excel col names for $component: $(comp_sheet[1,:])"
                                col = findall(x->x==component_id, comp_sheet[1,:])
                                #println("col: $col")
                                @assert length(col) == 1 && typeof(col[1]) == CartesianIndex{2}
                                col = col[1][2]
                                # col = findall(x->x==component_id, comp_sheet[1,:])[1]  # col index in the Excel file
                                # println("col: $col")
                                @assert col > 3
                                # Copy csv data in the sheet
                                comp_sheet[n_rows+1:n_rows+n_hours, col] = float(csv_data[:, component_id]) * base_value
                            end
                            comp_sheet[n_rows+1:n_rows+n_hours, 1] = 1:n_hours
                            comp_sheet[n_rows+1:n_rows+n_hours, 2] = [micro_scenario_name for _ in 1:n_hours]
                            comp_sheet[n_rows+1:n_rows+n_hours, 3] = [attribute for _ in 1:n_hours]
                        end
                    end
                end
            end
        end
    end
end

function build_outputs_from_json(json_path::String)
    main_attributes = Dict(
        "branch" => ["pf", "pt"],
        "branchdc" => ["i_from", "i_to"],
        "bus" => ["va"],
        "busdc" => ["vm"],
        "convdc" => ["iconv_dc", "iconv_dcg", "iconv_dcg_shunt", "pconv", "pdc", "pdcg", "pgrid", "ppr_fr", "ptf_to", "vaconv", "vafilt"],
        "gen" => ["pg", "pgcurt"],
        "load" => ["pcurt", "pflex"],
        )
    file_name = basename(json_path)
    case_name = file_name[1:length(file_name)-5]
    raw_results = JSON.parsefile(json_path)
    nw = raw_results["solution"]["nw"]

    # Sorting of the hours (to avoid having 1 10 ... 19 2 ...)
    hours = sort!([parse(Int, h) for h in keys(nw)])
    hours = ["$h" for h in hours]
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
                    if !in(id, full_grid[component_type])
                        push!(full_grid[component_type], id)
                    end
                end
            end
        end
    end

    excel_path = joinpath(dirname(json_path), "$case_name.xlsx")
    XLSX.openxlsx(excel_path, mode="w") do xf  # Create Excel
        for component in keys(main_attributes)
            if component in keys(nw[h1])
                if "Sheet1" in XLSX.sheetnames(xf)  # Rename by default sheet
                    comp_sheet = xf["Sheet1"]
                    XLSX.rename!(comp_sheet, component)
                else
                    XLSX.addsheet!(xf, component)  # Add sheet
                    comp_sheet = xf[component]
                end
                comp_sheet[1, 1:3] = ["Time", "Micro-scenario", "Attribute"]
                # Add column_names
                for (col, component_id) in enumerate(keys(full_grid[component]))
                    comp_sheet[1, col+3] = component_id
                end
                row = 1
                for h in keys(nw)
                    if component in keys(nw[h])
                        if length(keys(nw[h][component])) > 0
                            comp_id1 = first(keys(nw[h][component]))
                            for attribute in main_attributes[component]
                                if nw[h][component][comp_id1][attribute] isa Vector
                                    n_values = length(nw[h][component][comp_id1][attribute])
                                    poles = ["-p", "-n", "-r"][1:n_values]
                                else
                                    n_values = 1
                                    poles = [""]
                                end
                                for k in 1:n_values
                                    pole = poles[k]
                                    comp_sheet[row+k, 1] = h
                                    comp_sheet[row+k, 2] = case_name
                                    comp_sheet[row+k, 3] = "$attribute$pole"
                                end

                                for component_id in keys(nw[h][component])
                                    #=
                                    if !in(component_id, comp_sheet[1,:])
                                        println("component_id: $(typeof(component_id))")
                                        println("comp_sheet[1,4]: $(comp_sheet[1,4])")
                                        println("comp_sheet[1,4]: $(typeof(comp_sheet[1,4]))")
                                    end
                                    =#
                                    column_titles = ["$id" for id in comp_sheet[1,:]]
                                    @assert component_id in column_titles  "component_id $component_id is not in the Excel col names for $component: $(comp_sheet[1,:])"
                                    col = findall(x->x==component_id, column_titles)
                                    @assert length(col) == 1 && typeof(col[1]) == CartesianIndex{2}
                                    col = col[1][2]
                                    @assert attribute in keys(nw[h][component][component_id])

                                    value = nw[h][component][component_id][attribute]
                                    if n_values == 1
                                        values = [value]
                                    else
                                        values = value
                                    end
                                    for value in values
                                        row = row + 1
                                        if first(attribute) == 'p'  # pu powers must be converted into MW
                                            comp_sheet[row, col] = float(value) * base_mva
                                        else
                                            comp_sheet[row, col] = value
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end