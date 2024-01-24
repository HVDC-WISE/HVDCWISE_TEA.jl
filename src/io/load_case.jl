# load_case.jl

## Import relevant packages
import FlexPlan as _FP
import PowerModels as _PM
import CSV
using  DataFrames
import HiGHS
import HVDCWISE_TEA as _HWTEA

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

function load_case(test_case_name)
    # define optimizer linear solver
    optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    ## Load test case
    file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.m_files\\$test_case_name.m")

    ## Read CSV files
    # Non-dispatchable generation related files
    ndgen_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\ndgenerators_MW.csv")

    # Load related files
    fix_loads_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\loads_fix_MW.csv")
    flex_loads_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\loads_flex_MW.csv")

    # Storage related file
    # TODO process those data into time_series["storage"]["id"]["stationary_energy_inflow\outflow"]["timestamp"] (inflow = X & outflow = 0 if X>=0, else outflow = -X & inflow = 0 )
    hydro_dam_storage_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\hydro_dam_MW.csv")

    # Statuses file
    generation_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\generators_statuses.csv")
    ac_lines_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\AC_lines_statuses.csv")
    dc_lines_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\DC_lines_statuses.csv")
    converters_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\converters_statuses.csv")
    pst_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\PST_statuses.csv")
    storage_statuses_path = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\storage_statuses.csv")
    
    global time_series = Dict{String, Any}()
        
    load_files = [fix_loads_path,flex_loads_path]
 
    # Modifying the generation status with the timeseries data
    gen_statuses_data = read_csv_file(generation_statuses_path)[:,2:end]
    if !isempty(gen_statuses_data)
        for (_, col) in enumerate(names(gen_statuses_data))
            for (row_idx, value) in enumerate(gen_statuses_data[!, col])
                timestamp= Int64(row_idx)
                gen_id = col
                gen_status = value
                if !haskey(time_series, "gen")
                    time_series["gen"] = Dict{String, Any}()
                end
                if !haskey(time_series["gen"], "$gen_id")
                    time_series["gen"]["$gen_id"]= Dict{String, Any}()
                end
                if !haskey(time_series["gen"]["$gen_id"], "gen_status")
                    time_series["gen"]["$gen_id"]["gen_status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["gen"][gen_id]["gen_status"]) < timestamp
                    push!(time_series["gen"][gen_id]["gen_status"], 1)
                end
                time_series["gen"]["$gen_id"]["gen_status"][timestamp] = gen_status 
            end
        end
    end

    # Setting generation and load values
    set_generation_values!(time_series, ndgen_path, gen_statuses_data)
    set_loads_values!(time_series, load_files) 

    # Modifying the AC lines status with the timeseries data
    ac_lines_statuses_data = read_csv_file(ac_lines_statuses_path)[:,2:end]
    if !isempty(ac_lines_statuses_data)
        for (_, col) in enumerate(names(ac_lines_statuses_data))
            for (row_idx, value) in enumerate(ac_lines_statuses_data[!, col])
                timestamp = Int64(row_idx)
                branch_id = col
                br_status = value
                if !haskey(time_series, "branch")
                    time_series["branch"] = Dict{String, Any}()
                end
                if !haskey(time_series["branch"], "$branch_id")
                    time_series["branch"]["$branch_id"]= Dict{String, Any}()
                end
                if !haskey(time_series["branch"]["$branch_id"], "br_status")
                    time_series["branch"]["$branch_id"]["br_status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["branch"]["$branch_id"]["br_status"]) < timestamp
                    push!(time_series["branch"]["$branch_id"]["br_status"], 1)
                end
                time_series["branch"]["$branch_id"]["br_status"][timestamp] = br_status
            end
        end
    end

    # Modifying the storage status with the timeseries data
    storage_statuses_data = read_csv_file(storage_statuses_path)[:,2:end]
    if !isempty(storage_statuses_data)
        for (_, col) in enumerate(names(storage_statuses_data))
            for (row_idx, value) in enumerate(storage_statuses_data[!, col])
                timestamp = Int64(row_idx)
                storage_id = col
                storage_status = value
                if !haskey(time_series, "storage")
                    time_series["storage"] = Dict{String, Any}()
                end
                if !haskey(time_series["storage"], "$storage_id")
                    time_series["storage"]["$storage_id"]= Dict{String, Any}()
                end
                if !haskey(time_series["storage"]["$storage_id"], "status")
                    time_series["storage"]["$storage_id"]["status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["storage"]["$storage_id"]["status"]) < timestamp
                    push!(time_series["storage"]["$storage_id"]["status"], 1)
                end
                time_series["storage"]["$storage_id"]["status"][timestamp] = storage_status
            end
        end
    end

    # Modifying the pst status with the timeseries data
    pst_statuses_data = read_csv_file(pst_statuses_path)[:,2:end]
    if !isempty(pst_statuses_data)
        for (_, col) in enumerate(names(pst_statuses_data))
            for (row_idx, value) in enumerate(pst_statuses_data[!, col])
                timestamp = Int64(row_idx)
                pst_id = col
                pst_status = value
                if !haskey(time_series, "pst")
                    time_series["pst"] = Dict()
                end
                if !haskey(time_series["pst"], "$pst_id")
                    time_series["pst"]["$pst_id"]= Dict{String, Any}()
                end
                if !haskey(time_series["pst"]["$pst_id"], "br_status")
                    time_series["pst"]["$pst_id"]["br_status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["pst"]["$pst_id"]["br_status"]) < timestamp
                    push!(time_series["pst"]["$pst_id"]["br_status"], 1)
                end
                time_series["pst"]["$pst_id"]["br_status"][timestamp] = pst_status
            end
        end
    end
    
    # Modifying the DC lines status with the timeseries data
    dc_line_statuses_data = read_csv_file(dc_lines_statuses_path)
    if !isempty(dc_line_statuses_data)
        # Extract header information to check the phase order in the CSV
        header_phases_dc = [split(name, "_")[2] for name in names(dc_line_statuses_data)[2:end]]
        expected_phases_dc = [i % 3 == 1 ? "+" : i % 3 == 2 ? "-" : "MR" for i in 1:length(header_phases_dc)]

        if expected_phases_dc != header_phases_dc
            error("Unexpected phase order in DC line status file header.\n Order X_+;X_-;X_MR for component X expected")
        end

        for (row_idx, timestamp_str) in enumerate(dc_line_statuses_data[!, 1])
            timestamp = Int64(timestamp_str)
            for col_group_idx in 1:3:size(dc_line_statuses_data, 2)-2  # iterating by group of 3 columns (+, - , MR) for each DC lines
                # Branchdc_id is inferred from the column index
                branchdc_id = div(col_group_idx, 3) + 1

                # Separate variables for each status type
                status_p = dc_line_statuses_data[row_idx, col_group_idx + 1] # positive phase is the first column of the group
                status_n = dc_line_statuses_data[row_idx, col_group_idx + 2] # negative phase is the second column of the group
                status_r = dc_line_statuses_data[row_idx, col_group_idx + 3] # metallic return phase is the third column of the group
                
                # Assign values to the corresponding fields
                if !haskey(time_series, "branchdc")
                    time_series["branchdc"] = Dict{String, Any}()
                end
                if !haskey(time_series["branchdc"], "$branchdc_id")
                    time_series["branchdc"]["$branchdc_id"] = Dict{String, Any}()
                end
                if !haskey(time_series["branchdc"]["$branchdc_id"], "status")
                    time_series["branchdc"]["$branchdc_id"]["status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["branchdc"]["$branchdc_id"]["status"]) < timestamp
                    push!(time_series["branchdc"]["$branchdc_id"]["status"], 1)
                end
                time_series["branchdc"]["$branchdc_id"]["status"][timestamp] = [status_p, status_n, status_r]
            end
        end
    end

    # Modifying the ACDC converters status with the timeseries data
    converters_statuses_data = read_csv_file(converters_statuses_path)
    if !isempty(converters_statuses_data)
        
        # Extract header information to check the phase order in the CSV
        header_phases_conv = [split(name, "_")[2] for name in names(converters_statuses_data)[2:end]]
        expected_phases_conv = [i % 2 == 1 ? "+" : "-" for i in 1:length(header_phases_conv)]

        if expected_phases_conv != header_phases_conv
            error("Unexpected phase order in converter status file header.\n Order X_+;X_-; for component X expected")
        end

        for (row_idx, timestamp_str) in enumerate(converters_statuses_data[!, 1])
            timestamp = Int64(timestamp_str)
            for col_group_idx in 1:2:size(converters_statuses_data, 2)-1 # iterating by group of 2 columns (+, -) for each AC/DC conv
                # Conv_id is inferred from the column index
                conv_id = div(col_group_idx, 2) + 1

                # Separate variables for each status type
                status_p = converters_statuses_data[row_idx, col_group_idx + 1] # positive converter is the first column of the group
                status_n = converters_statuses_data[row_idx, col_group_idx + 2] # negative converter is the second column of the group

                # Assign values to the corresponding fields
                if !haskey(time_series, "convdc")
                    time_series["convdc"] = Dict{String, Any}()
                end
                if !haskey(time_series["convdc"], "$conv_id")
                    time_series["convdc"]["$conv_id"] = Dict{String, Any}()
                end
                if !haskey(time_series["convdc"]["$conv_id"], "status")
                    time_series["convdc"]["$conv_id"]["status"] = Vector{Any}()
                end
                # Make sure the vector is large enough to accommodate the timestamp
                while length(time_series["convdc"]["$conv_id"]["status"]) < timestamp
                    push!(time_series["convdc"]["$conv_id"]["status"], 1)
                end
                time_series["convdc"]["$conv_id"]["status"][timestamp] = [status_p, status_n]
            end
        end
    end

    global data = _HWTEA.parse_data(file, time_series) # read the data and add multinetwork dimension

    # HVDCWiseTEA settings
    s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => false)

    # TEST
    _PM.propagate_topology_status!(data)  # TODO check if usefull
    #data["per_unit"] = false
    #_PM.make_mixed_units!(data) # TODO check if usefull
    #_PM.make_per_unit!(data) # TODO check if usefull

   return _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = s)
end

# Function to read CSV file if it exists
function read_csv_file(file_path)
    if isfile(file_path)
        return CSV.File(file_path, delim=';') |> DataFrames.DataFrame
    else
        return DataFrames.DataFrame()  # Return an empty DataFrame if the file doesn't exist
    end
end

# Function to set the generation values accordingt to generation files
function set_generation_values!(time_series, gen_file, gen_statuses_data)    
    gen_data = read_csv_file(gen_file)
    if !isempty(gen_data)
        # Extract generator IDs (excluding the first column 'Time')
        gen_ids = string.(names(gen_data[:, 2:end]))
        # Extract timestamps
        timestamps = gen_data.Time

        # Iterate over timestamps
        for (timestamp_idx, timestamp) in enumerate(timestamps)
            # Iterate over generator IDs
            for (gen_id_idx, gen_id) in enumerate(gen_ids)
                gen_value = gen_data[!, gen_id_idx + 1][timestamp_idx]  # Extract gen_value for the current gen_id

                # Update time_series with aggregated generation values using haskey
                if !haskey(time_series, "gen")
                    time_series["gen"] = Dict{String, Any}()
                end
                if !haskey(time_series["gen"], gen_id)
                    time_series["gen"][gen_id] = Dict{String, Any}()
                end
                if !haskey(time_series["gen"][gen_id], "pmax")
                    time_series["gen"][gen_id]["pmax"] = Vector{Any}()
                end
                if !haskey(time_series["gen"][gen_id], "gen_status")
                    time_series["gen"][gen_id]["gen_status"] = Vector{Any}()
                end

                # Make sure the vectors are large enough to accommodate the timestamp
                while length(time_series["gen"][gen_id]["pmax"]) < timestamp
                    push!(time_series["gen"][gen_id]["pmax"], 0)
                end
                while length(time_series["gen"][gen_id]["gen_status"]) < timestamp
                    push!(time_series["gen"][gen_id]["gen_status"], 1)
                end
                
                if !isempty(gen_statuses_data)
                    # Check gen_status at the current timestamp
                    gen_status = gen_statuses_data[!, gen_id][timestamp_idx]

                    if gen_status == 1
                        # Set pmax if gen_status is 1
                        time_series["gen"][gen_id]["pmax"][timestamp] += gen_value
                    else
                        # Set gen_status if gen_status is not 1
                        time_series["gen"][gen_id]["gen_status"][timestamp] = 0
                    end
                else
                    # Set pmax
                    time_series["gen"][gen_id]["pmax"][timestamp] += gen_value
                end        
            end
        end
    end
end



# Function to set the loads values of the several load files
function set_loads_values!(time_series, load_files)
    for load_file in load_files
        load_data = read_csv_file(load_file)
        if !isempty(load_data)
            # Extract generator IDs (excluding the first column 'Time')
            load_ids = string.(names(load_data[:, 2:end]))
            # Extract timestamps
            timestamps = load_data.Time

            # Iterate over timestamps
            for (timestamp_idx, timestamp) in enumerate(timestamps)
                # Iterate over generator IDs
                for (load_id_idx, load_id) in enumerate(load_ids)
                    load_value = load_data[!, load_id_idx + 1][timestamp_idx]  # Extract load_value for the current gen_id

                    # Update time_series with aggregated generation values using haskey
                    if !haskey(time_series, "load")
                        time_series["load"] = Dict{String, Any}()
                    end
                    if !haskey(time_series["load"], load_id)
                        time_series["load"][load_id] = Dict{String, Any}()
                    end
                    if !haskey(time_series["load"][load_id], "pd")
                        time_series["load"][load_id]["pd"] = Vector{Any}()
                    end
                    # Make sure the vector is large enough to accommodate the timestamp
                    while length(time_series["load"][load_id]["pd"]) < timestamp
                        push!(time_series["load"][load_id]["pd"], 0)
                    end
                    time_series["load"][load_id]["pd"][timestamp] += load_value
                end
            end
        end
    end
end

