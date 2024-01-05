# load_case.jl

## Import relevant packages
import FlexPlan as _FP
import PowerModels as _PM
import CSV
import DataFrames
import HiGHS
import HVDCWISE_TEA as _HWTEA

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package


function load_case(test_case_name)
    optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => true)

    ## Load test case
    file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.m_files\\$test_case_name.m")
    data = _FP.parse_file(file) # Parse input file to obtain data dictionary

    # Read CSV files
    demand_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\load_fix_MW.csv")
    generation_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_generators.csv")
    AC_line_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_lines_AC.csv")
    DC_line_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_lines_DC.csv")
    converters_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_converters.csv")

    demand = CSV.read(demand_file, DataFrames.DataFrame)[:, 2:end]
    demand_pu = demand[:, 1]
    
    for (l,load) in data["load"]
        demand_pu = demand_pu[:, parse(Int, l)] ./ load["pd"]   

    end 
    
    # get all rows of the 2nd column, the first row is the timestamp & the index of the generators in the mpc.generator table in the .m file
    generation_status = CSV.read(generation_status_file,DataFrames.DataFrame)[:,2:end]

    # get all rows of the 2nd column, the first row is the timestamp & the index of the AC line in the mpc.generator table in the .m file
    AC_line_status = CSV.read(AC_line_status_file,DataFrames.DataFrame)[:,2:end]
        
    # Get the number of time points
    number_of_hours = size(CSV.read(joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\load_fix_MW.csv"), DataFrames.DataFrame), 1)

    # Initialize genprofile arrays
    genprofile = ones(length(data["gen"]), number_of_hours)
    # Populate genprofile based on generation status CSV file
    genprofile[:, :] .= generation_status[:, 1]'

    # Initialize loadprofile arrays
    loadprofile = ones(length(data["load"]), number_of_hours)
    # Populate loadprofile based on demand CSV file
    loadprofile[:, :] .= demand_pu'

    # Create time series data to be passed to the data dictionary
    time_series = _FP.make_time_series(data, number_of_hours; loadprofile = permutedims(loadprofile), genprofile = permutedims(genprofile))
    data = _HWTEA.parse_data(file, time_series) # reread the data and add multinetwork dimension

    # Modifying the AC line status with the timeseries data
    for (col_idx, col) in enumerate(names(AC_line_status))
        for (row_idx, value) in enumerate(AC_line_status[!, col])
            timestamp = row_idx
            branch_id = col
            br_status = value
            data["nw"]["$timestamp"]["branch"]["$branch_id"]["br_status"] = br_status
        end
    end

    # Modifying the DC lines status with the timeseries data
    if isfile(DC_line_status_file)
        DC_line_status = CSV.read(DC_line_status_file, DataFrames.DataFrame)
        
        # Extract header information to check the phase order in the CSV
        header_phases_dc = [split(name, "_")[2] for name in names(DC_line_status)[2:end]]
        expected_phases_dc = [i % 3 == 1 ? "+" : i % 3 == 2 ? "-" : "MR" for i in 1:length(header_phases_dc)]

        if expected_phases_dc != header_phases_dc
            error("Unexpected phase order in DC line status file header.\n Order X_+;X_-;X_MR for component X expected")
        end

        for (row_idx, timestamp) in enumerate(DC_line_status[!, 1])
            for col_group_idx in 1:3:size(DC_line_status, 2)-2  # iterating by group of 3 columns (+, - , MR) for each DC lines
                # Branchdc_id is inferred from the column index
                branchdc_id = div(col_group_idx, 3) + 1

                # Separate variables for each status type
                status_p = DC_line_status[row_idx, col_group_idx + 1] # positive phase is the first column of the group
                status_n = DC_line_status[row_idx, col_group_idx + 2] # negative phase is the second column of the group
                status_r = DC_line_status[row_idx, col_group_idx + 3] # metallic return phase is the third column of the group

                # Assign values to the corresponding fields
                data["nw"]["$timestamp"]["branchdc"]["$branchdc_id"]["status_p"] = status_p
                data["nw"]["$timestamp"]["branchdc"]["$branchdc_id"]["status_n"] = status_n
                data["nw"]["$timestamp"]["branchdc"]["$branchdc_id"]["status_r"] = status_r
            end
        end
    end

    # Modifying the ACDC converters status with the timeseries data
    if isfile(converters_status_file)
        converters_status = CSV.read(converters_status_file, DataFrames.DataFrame)
        
        # Extract header information to check the phase order in the CSV
        header_phases_conv = [split(name, "_")[2] for name in names(converters_status)[2:end]]
        expected_phases_conv = [i % 2 == 1 ? "+" : "-" for i in 1:length(header_phases_conv)]

        if expected_phases_conv != header_phases_conv
            error("Unexpected phase order in converter status file header.\n Order X_+;X_-; for component X expected")
        end

        for (row_idx, timestamp) in enumerate(converters_status[!, 1])
            for col_group_idx in 1:2:size(converters_status, 2)-1 # iterating by group of 2 columns (+, -) for each AC/DC conv
                # Conv_id is inferred from the column index
                conv_id = div(col_group_idx, 2) + 1

                # Separate variables for each status type
                status_p = converters_status[row_idx, col_group_idx + 1] # positive converter is the first column of the group
                status_n = converters_status[row_idx, col_group_idx + 2] # negative converter is the second column of the group

                # Assign values to the corresponding fields
                data["nw"]["$timestamp"]["convdc"]["$conv_id"]["status_p"] = status_p
                data["nw"]["$timestamp"]["convdc"]["$conv_id"]["status_n"] = status_n
            end
        end
    end

    # HVDCWiseTEA settings
    s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => false)

    # TEST
    _PM.propagate_topology_status!(data)  # TODO check if usefull
    #data["per_unit"] = false
    #_PM.make_mixed_units!(data) # TODO check if usefull
    #_PM.make_per_unit!(data) # TODO check if usefull

   return _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = s)
end
