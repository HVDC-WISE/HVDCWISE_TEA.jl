import CSV
import DataFrames
import JSON
import XLSX


function build_user_results(work_dir::String, base_mva::Int, matlab_octave_path::String)
    # Identify simulation results directory
    simulation_dir = joinpath(work_dir, "simulation_interface")
    macro_scenario = ""
    for fifo_name in readdir(simulation_dir)
        dir = joinpath(simulation_dir, fifo_name)
        if isdir(dir) && fifo_name != "Input_series"
            @assert (macro_scenario == "") "Impossible to identify the simulation results directory.\n$(joinpath(simulation_dir, macro_scenario))\n$dir"
            macro_scenario = fifo_name
        end
    end

    # Write simulation results folder name in a file of the matlab tool folder
    matlab_tool_path = joinpath(pwd(), "src", "matlab_tools")
    write(joinpath(matlab_tool_path,"simulation_results_path.txt"), joinpath(simulation_dir, macro_scenario))

    # Build user results (Julia code)
    user_results_dir = joinpath(work_dir, "user_interface", "results")
    if isdir(user_results_dir)
        rm(user_results_dir, recursive=true)
    end
    mkpath(user_results_dir)
    opf_results = gather_opf_results(work_dir, macro_scenario, base_mva)
    opex = opex_summary(opf_results)
    capex = capex_summary(work_dir)
    totex = totex_summary(opex, capex)

    # Build user results (Matlab code)
    # println("Run compute_KPIs.m in $matlab_tool_path. Then write 'y' and press twice ENTER.")
    # a = readline();  # TODO automatically run src/matlab_tools/compute_KPIs.m or recode it in Julia
    println("Computing KPIs")
    run_matlab_script(joinpath(matlab_tool_path, "compute_KPIs.m"), matlab_octave_path)

    # Delete study files in the matlab tool folder
    rm(joinpath(matlab_tool_path, "grid_model.m"))
    rm(joinpath(matlab_tool_path,"simulation_results_path.txt"))
end


function gather_opf_results(work_dir::String, macro_scenario::String, base_mva::Int)
    simulation_results_dir = joinpath(work_dir, "simulation_interface", macro_scenario)

    result_folders = Vector{String}()
    for micro_scenario in readdir(simulation_results_dir)
        micro_dir = joinpath(simulation_results_dir, micro_scenario)
        if isdir(micro_dir)
            push!(result_folders, micro_dir)
        end
    end

    opf_results = Dict()

    for result_folder in result_folders
        micro_scenario_name = basename(result_folder)
        mkpath(joinpath(work_dir, "user_interface", "results", micro_scenario_name))
        excel_path = joinpath(work_dir, "user_interface", "results", micro_scenario_name, "OPF_results.xlsx")
        XLSX.openxlsx(excel_path, mode="w") do xf  # Create Excel
            opf_results[micro_scenario_name] = Dict()
            for component in readdir(result_folder)
                comp_folder = joinpath(result_folder, component)
                if isdir(comp_folder)
                    # component = basename(comp_folder)
                    opf_results[micro_scenario_name][component] = Dict()

                    if component in XLSX.sheetnames(xf)
                        comp_sheet = xf[component]
                        @assert XLSX.get_dimension(comp_sheet).stop.row_number > 1
                    else
                        if "Sheet1" in XLSX.sheetnames(xf)  # Rename by default sheet
                            comp_sheet = xf["Sheet1"]
                            XLSX.rename!(comp_sheet, component)
                        else
                            XLSX.addsheet!(xf, component)  # Add sheet
                            comp_sheet = xf[component]
                        end
                        comp_sheet[1, 1:4] = ["Time", "Micro-scenario", "Attribute", "Unit"]
                    end
                    for file_name in readdir(comp_folder)
                        if occursin(".csv", file_name)
                            attribute_file = joinpath(comp_folder, file_name)
                            attribute = file_name[1:length(file_name)-4]
                            opf_results[micro_scenario_name][component][attribute] = Dict()
                            # Read csv
                            csv_data = CSV.File(attribute_file, delim=',') |> DataFrames.DataFrame
                            n_hours = first(size(csv_data))
                            n_rows = XLSX.get_dimension(comp_sheet).stop.row_number

                            if n_rows in [0, 1]
                                for (col, component_id) in enumerate(names(csv_data))
                                    comp_sheet[1, col+4] = component_id
                                end
                                n_rows = 1
                            end

                            if first(attribute) == 'p'
                                # pu powers must be converted into MW
                                base_value = base_mva
                                unit = "MW"
                            else
                                base_value = 1
                                unit = "pu"
                            end
                            for component_id in names(csv_data)
                                @assert component_id in comp_sheet[1,:]  "component_id $component_id is not in the Excel col names for $component: $(comp_sheet[1,:])"
                                col = findall(x->x==component_id, comp_sheet[1,:])
                                @assert length(col) == 1 && typeof(col[1]) == CartesianIndex{2}
                                col = col[1][2]
                                @assert col > 3
                                values = float(csv_data[:, component_id]) * base_value
                                comp_sheet[n_rows+1:n_rows+n_hours, col] = values
                                opf_results[micro_scenario_name][component][attribute][component_id] = values
                            end
                            comp_sheet[n_rows+1:n_rows+n_hours, 1] = 1:n_hours
                            comp_sheet[n_rows+1:n_rows+n_hours, 2] = [micro_scenario_name for _ in 1:n_hours]
                            comp_sheet[n_rows+1:n_rows+n_hours, 3] = [attribute for _ in 1:n_hours]
                            comp_sheet[n_rows+1:n_rows+n_hours, 4] = [unit for _ in 1:n_hours]
                        end
                    end
                end
            end
        end
    end
    return opf_results
end


function opex_summary(opf_results::Dict)  # TODO
    opex = Dict()

    return opex
end

function capex_summary(work_dir::String)  # TODO
    capex = Dict()

    return capex
end

function totex_summary(opex::Dict, capex::Dict)  # TODO
    totex = Dict()

    return totex
end
