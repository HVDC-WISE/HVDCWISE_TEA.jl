using CSV
import DataFrames
using XLSX

function build_raw_inputs(input_dir::String, output_dir::String, case_name::String, base_mva::Int)
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
                    CSV.write(csv_path, DataFrame(sheet[2:end, 2:end],:auto),writeheader=false)  # auto is needed to convert matrix into table
                end
            end
        end
    end
end
