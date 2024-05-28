using CSV
import DataFrames
using XLSX

function build_raw_inputs(input_dir::String, output_dir::String, macro_scenario::String, base_mva::Int, micro_scenarios::Bool=false)
    build_m(input_dir, output_dir, macro_scenario, base_mva)
    if micro_scenarios  # in run_study
        work_dir = joinpath(output_dir, macro_scenario)
        if isdir(work_dir)
            rm(work_dir; force=true, recursive=true)
        end
        mkpath(work_dir)
        for file_name in readdir(input_dir)
            if occursin("_series.xlsx", file_name)
                micro_scenario = file_name[1:length(file_name)-12]
                micro_dir = joinpath(work_dir, micro_scenario)
                mkpath(micro_dir)
                build_csv(input_dir, micro_dir, micro_scenario, base_mva, true)
            end
        end
    else  # in test_load_case
        micro_scenario = macro_scenario
        build_csv(input_dir, output_dir, micro_scenario, base_mva, false)
    end
end
function build_m(input_dir::String, output_dir::String, macro_scenario::String, base_mva::Int)  # .m file building
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    useful_model_sheets = Dict(
    "bus" => "AC buses",
    "busdc" => "DC buses",
    "branch" => "AC branches",
    "branch_currents" => "AC branch rating",
    "branchdc" => "DC branches",
    "convdc" => "converters",
    "pst" => "phase-shift transformers",
    "gen" => "generators",
    "gencost" => "dispatchable generators operating cost model",
    "ndgen" => "non-dispatchable generators",
    "load_extra" => "load additional data",
    "storage" => "storages",
    "storage_extra" => "storage additional data"
    )

    excel_path = joinpath(input_dir, "$macro_scenario"*"_model.xlsx")
    @assert isfile(excel_path) "$excel_path is not a file"
    excel_file = XLSX.readxlsx(excel_path)
    matpower_path = joinpath(output_dir, "$macro_scenario.m")
    if isfile(matpower_path)
        rm(matpower_path; force=true)
    end
    open(matpower_path, "w") do f
        write(f, "function mpc = $macro_scenario()\n\n")
        write(f, "%% MATPOWER Case Format : Version 2\nmpc.version = '2';\n\n")
        write(f, "%% system MVA base\nmpc.baseMVA = $base_mva;\n\n")
        write(f, "mpc.time_elapsed = 1.0;\n")
        sheet_names = XLSX.sheetnames(excel_file)
        for sheet_name in sheet_names
            if haskey(useful_model_sheets, sheet_name)
                sheet = excel_file[sheet_name]
                n_rows = XLSX.get_dimension(sheet).stop.row_number
                n_cols = XLSX.get_dimension(sheet).stop.column_number
                if n_rows > 3 && n_cols > 1  # 3 first rows are for Attribute, Unit and Description. First col is for row title and comments
                    @assert [sheet[i,1] for i in 1:3] == ["Attribute", "Unit", "Description"] "$excel_path sheet $sheet_name A1:A3 is $([sheet[i,1] for i in 1:3]) instead of [Attribute, Unit, Description]"
                    if sum([string(sheet[row,col]) != "missing" for row=4:n_rows for col=2:n_cols]) != 0  # sheet not empty
                        write(f, "\n% $(useful_model_sheets[sheet_name])\n")
                        write(f, "%column_names% $(join(sheet[1,2:end], " "))\n")
                        write(f, "mpc.$(sheet_name) = [\n")
                        for row in 4:n_rows
                            row_values = [string(sheet[row,col]) for col in 2:n_cols]
                            if sum([value != "missing" for value in row_values]) > 0  # line not empty
                                n_missing_values = sum([value == "missing" for value in row_values])
                                @assert n_missing_values == 0 "$n_missing_values values are missing in sheet $sheet_name of $excel_path"
                                write(f, join(row_values, " ") * ";\n")
                            end
                        end
                        write(f, "];\n")
                    end
                end
            end
        end
        write(f, "end\n")
    end
end

function build_csv(input_dir::String, output_dir::String, micro_scenario::String, base_mva::Int, pu::Bool=false)  # .csv files building
    useful_series_sheets = [
        "branch|status",
        "branchdc|status_n",
        "branchdc|status_p",
        "branchdc|status_r",
        "convdc|status_n",
        "convdc|status_p",
        "gen|Pmax",
        "gen|status",
        "load|Pd",
        "storage|energy_inflow",  # stationary_energy_inflow
        "storage|energy_outflow",  # stationary_energy_outflow
        "storage|status"
        ]
    
        attribute_conversion = Dict(
            "branch" => Dict("status" => "br_status"),
            "branchdc" => Dict("status_n" => "status_n", "status_p" => "status_p", "status_r" => "status_r"),
            "convdc" => Dict("status_n" => "status_n", "status_p" => "status_p"),
            "gen" => Dict("Pmax" => "pmax", "status" => "gen_status"),
            "load"=> Dict("Pd" => "pd"),
            "storage" => Dict("status" => "status", "energy_inflow" => "stationary_energy_inflow", "energy_outflow" => "stationary_energy_outflow")
        )
    
        for component_type in ["branch", "branchdc", "convdc", "gen", "load", "storage"]
            comp_dir_path = joinpath(output_dir, component_type)
            if isdir(comp_dir_path)
                rm(comp_dir_path; force=true, recursive=true)
            end
        end
        excel_path = joinpath(input_dir, "$micro_scenario"*"_series.xlsx")
        excel_file = XLSX.readxlsx(excel_path)
        sheet_names = XLSX.sheetnames(excel_file)
        for sheet_name in sheet_names
            if sheet_name in useful_series_sheets
                component_type, attribute = split(sheet_name, "|")
                attribute = attribute_conversion[component_type][attribute]
                sheet = excel_file[sheet_name]
                n_rows = XLSX.get_dimension(sheet).stop.row_number
                n_cols = XLSX.get_dimension(sheet).stop.column_number
                if n_rows > 2 && n_cols > 1  # Second row is for component id. First column (for time) is not conserved.
                    @assert sheet[2,1] == "Time\\Id" "$excel_path $sheet_name A2 should be Time\\Id"
                    if sum([string(sheet[row,col]) == "missing" for row=3:n_rows for col=2:n_cols]) == 0  # sheet not empty
                        #=
                        real_n_rows = 0
                        row = 1
                        while row < n_rows && real_n_rows == 0
                            if sum([string(sheet[row,col]) == "missing" for col=2:n_cols]) > 0
                                real_n_rows = row
                            else
                                row = row + 1
                            end
                        end
                        =#
                        comp_dir_path = joinpath(output_dir, component_type)
                        if isdir(comp_dir_path) == false
                            mkdir(comp_dir_path)
                        end
                        csv_path = joinpath(comp_dir_path, "$attribute.csv")
                        if pu && attribute in ["pmax", "pd", "stationary_energy_inflow", "stationary_energy_outflow"]  # power or energy data
                        # if pu && "$(first(attribute))" in ["p", "e"]  # power or energy data
                            df = DataFrames.DataFrame(Dict{String, Any}(
                                "$(sheet[2,j])" => [float(sheet[i,j]) /base_mva for i in 3:n_rows] for j in 2:n_cols)
                                )

                            #=
                            data = [[sheet[2,j] for j in 2:n_cols]]
                            for i in 3:n_rows
                                push!(data, [sheet[i,j] for j in 2:n_cols])
                            end
                            =#

                            #data = [[float(sheet[i,j]) / base_mva for j in 2:n_cols] for i in 3:n_rows]
                            #pushfirst!(data, [sheet[2,j] for j in 2:n_cols]) # adding row for component id

                            # println("n_rows: $n_rows. n_cols: $n_cols. attribute: $attribute")
                            # println([sheet[2,j] for j in 2:n_cols])
                            # pushfirst!(sheet[2,2:n_cols], data)
                            # data = float(DataFrames.DataFrame(sheet[2:end, 2:end],:auto)) / base_mva
                            # values = parse(Float64, sheet[2:end, 2:end])  # float(sheet[2:end, 2:end]) / base_mva
                        else
                            df = DataFrames.DataFrame(Dict{String, Any}(
                                "$(sheet[2,j])" => [sheet[i,j] for i in 3:n_rows] for j in 2:n_cols)
                                )
                            # data = sheet[2:n_rows, 2:n_cols]

                            # data = DataFrames.DataFrame(sheet[2:end, 2:end],:auto)
                            # println(first(attribute))
                            # values = sheet[2:end, 2:end]
                        end
                        #=
                        if pu && first(attribute) in ["p", "e"]  # power or energy data
                            values = float(sheet[2:end, 2:end]) / base_mva
                        else
                            values = sheet[2:end, 2:end]
                        end
                        =#
                        # CSV.write(csv_path, DataFrames.DataFrame(data,:auto), writeheader=false)  # auto is needed to convert matrix into table
                        CSV.write(csv_path, df, writeheader=true)  # auto is needed to convert matrix into table
                    end
                end
            end
        end
end

#=
function build_raw_inputs(input_dir::String, output_dir::String, case_name::String, base_mva::Int, pu::Bool=false)
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    ## .m file building
    
    useful_model_sheets = Dict(
    "bus" => "AC buses",
    "busdc" => "DC buses",
    "branch" => "AC branches",
    "branch_currents" => "AC branch rating",
    "branchdc" => "DC branches",
    "convdc" => "converters",
    "pst" => "phase-shift transformers",
    "gen" => "generators",
    "gencost" => "dispatchable generators operating cost model",
    "ndgen" => "non-dispatchable generators",
    "load_extra" => "load additional data",
    "storage" => "storages",
    "storage_extra" => "storage additional data"
    )

    excel_path = joinpath(input_dir, "$case_name"*"_model.xlsx")
    @assert isfile(excel_path) "$excel_path is not a file"
    excel_file = XLSX.readxlsx(excel_path)
    matpower_path = joinpath(output_dir, "$case_name.m")
    if isfile(matpower_path)
        rm(matpower_path; force=true)
    end
    open(matpower_path, "w") do f
        write(f, "function mpc = $case_name()\n\n")
        write(f, "%% MATPOWER Case Format : Version 2\nmpc.version = '2';\n\n")
        write(f, "%% system MVA base\nmpc.baseMVA = $base_mva;\n\n")
        write(f, "mpc.time_elapsed = 1.0;\n")
        sheet_names = XLSX.sheetnames(excel_file)
        for sheet_name in sheet_names
            if haskey(useful_model_sheets, sheet_name)
                sheet = excel_file[sheet_name]
                n_rows = XLSX.get_dimension(sheet).stop.row_number
                n_cols = XLSX.get_dimension(sheet).stop.column_number
                if n_rows > 3 && n_cols > 1  # 3 first rows are for Attribute, Unit and Description. First col is for row title and comments
                    @assert [sheet[i,1] for i in 1:3] == ["Attribute", "Unit", "Description"] "$excel_path sheet $sheet_name A1:A3 is $([sheet[i,1] for i in 1:3]) instead of [Attribute, Unit, Description]"
                    if sum([string(sheet[row,col]) != "missing" for row=4:n_rows for col=2:n_cols]) != 0  # sheet not empty
                        write(f, "\n% $(useful_model_sheets[sheet_name])\n")
                        write(f, "%column_names% $(join(sheet[1,2:end], " "))\n")
                        write(f, "mpc.$(sheet_name) = [\n")
                        for row in 4:n_rows
                            row_values = [string(sheet[row,col]) for col in 2:n_cols]
                            if sum([value != "missing" for value in row_values]) > 0  # line not empty
                                n_missing_values = sum([value == "missing" for value in row_values])
                                @assert n_missing_values == 0 "$n_missing_values values are missing in sheet $sheet_name of $excel_path"
                                write(f, join(row_values, " ") * ";\n")
                            end
                        end
                        write(f, "];\n")
                    end
                end
            end
        end
        write(f, "end\n")
    end

    ## .csv files building

    useful_series_sheets = [
    "branch|status",
    "branchdc|status_n",
    "branchdc|status_p",
    "branchdc|status_r",
    "convdc|status_n",
    "convdc|status_p",
    "gen|Pmax",
    "gen|status",
    "load|Pd",
    "storage|energy_inflow",  # stationary_energy_inflow
    "storage|energy_outflow",  # stationary_energy_outflow
    "storage|status"
    ]

    attribute_conversion = Dict(
        "branch" => Dict("status" => "br_status"),
        "branchdc" => Dict("status_n" => "status_n", "status_p" => "status_p", "status_r" => "status_r"),
        "convdc" => Dict("status_n" => "status_n", "status_p" => "status_p"),
        "gen" => Dict("Pmax" => "pmax", "status" => "gen_status"),
        "load"=> Dict("Pd" => "pd"),
        "storage" => Dict("status" => "status", "energy_inflow" => "stationary_energy_inflow", "energy_outflow" => "stationary_energy_outflow")
    )

    for component_type in ["branch", "branchdc", "convdc", "gen", "load", "storage"]
        comp_dir_path = joinpath(output_dir, component_type)
        if isdir(comp_dir_path)
            rm(comp_dir_path; force=true, recursive=true)
        end
    end
    excel_path = joinpath(input_dir, "$case_name"*"_series.xlsx")
    excel_file = XLSX.readxlsx(excel_path)
    sheet_names = XLSX.sheetnames(excel_file)
    for sheet_name in sheet_names
        if sheet_name in useful_series_sheets
            component_type, attribute = split(sheet_name, "|")
            attribute = attribute_conversion[component_type][attribute]
            sheet = excel_file[sheet_name]
            n_rows = XLSX.get_dimension(sheet).stop.row_number
            n_cols = XLSX.get_dimension(sheet).stop.column_number
            if n_rows > 2 && n_cols > 1  # Second row is for component id. First column (for time) is not conserved.
                @assert sheet[2,1] == "Time\\Id" "$excel_path $sheet_name A2 should be Time\\Id"
                if sum([string(sheet[row,col]) == "missing" for row=3:n_rows for col=2:n_cols]) == 0  # sheet not empty
                    comp_dir_path = joinpath(output_dir, component_type)
                    if isdir(comp_dir_path) == false
                        mkdir(comp_dir_path)
                    end
                    csv_path = joinpath(comp_dir_path, "$attribute.csv")
                    if pu && first(attribute) in ["p", "e"]  # power or energy data
                        values = float(sheet[2:end, 2:end]) / base_mva
                    else
                        values = sheet[2:end, 2:end]
                    end
                    CSV.write(csv_path, DataFrame(sheet[2:end, 2:end],:auto),writeheader=false)  # auto is needed to convert matrix into table
                end
            end
        end
    end
end
=#