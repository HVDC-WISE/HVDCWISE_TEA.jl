import CSV
import DataFrames
import JSON
import XLSX


function build_user_results(work_dir::String, base_mva::Int)
    # Identify simulation results directory
    simulation_dir = joinpath(work_dir, "simulation_interface")
    macro_scenario = ""
    for fifo_name in readdir(simulation_dir)
        dir = joinpath(simulation_dir, fifo_name)
        if isdir(dir) && fifo_name != "Inputs_series"
            @assert (macro_scenario == "") "Impossible to identify the simulation results directory.\n$(joinpath(simulation_dir, macro_scenario))\n$dir"
            macro_scenario = fifo_name
        end
    end

    # Build user results
    user_results_dir = joinpath(work_dir, "user_interface", "results")
    mkpath(user_results_dir)  # Create the folder if it does not exist yet
    gather_opf_results(work_dir, macro_scenario, base_mva)

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

    # Micro-scenarios results will be aggregated in the same Excel file
    excel_path = joinpath(work_dir, "user_interface", "results", "OPF_results.xlsx")
    
    XLSX.openxlsx(excel_path, mode="w") do xf  # Create Excel
        for result_folder in result_folders
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
                                @assert component_id in comp_sheet[1,:]  "component_id $component_id is not in the Excel col names for $component: $(comp_sheet[1,:])"
                                col = findall(x->x==component_id, comp_sheet[1,:])
                                @assert length(col) == 1 && typeof(col[1]) == CartesianIndex{2}
                                col = col[1][2]
                                @assert col > 3
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
