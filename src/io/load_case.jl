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
    DC_line_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_lines_DC_temp.csv")
    converters_status_file = joinpath(_HWTEA_dir, "test\\data\\$test_case_name\\.csv_files\\status_converters_temp.csv")

    demand = CSV.read(demand_file, DataFrames.DataFrame)[:, 2:end]
    demand_pu = demand[:, 1] ./ maximum(demand[:, 1])

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
        DC_line_status = CSV.read(DC_line_status_file, DataFrames.DataFrame)[:, 2:end]
        for (col_idx, col) in enumerate(names(DC_line_status))
            for (row_idx, value) in enumerate(DC_line_status[!, col])
                timestamp = row_idx
                branchdc_id = col
                status = value
                data["nw"]["$timestamp"]["branchdc"]["$branchdc_id"]["status"] = status
            end
        end
    end

    # Modifying the ACDC converters  status with the timeseries data
    if isfile(converters_status_file)
        converters_status = CSV.read(converters_status_file, DataFrames.DataFrame)[:, 2:end]
        for (col_idx, col) in enumerate(names(converters_status))
            for (row_idx, value) in enumerate(converters_status[!, col])
                timestamp = row_idx
                conv_id = col
                status = value
                data["nw"]["$timestamp"]["convdc"]["$conv_id"]["status"] = status
            end
        end
    end

    # HVDCWiseTEA settings
    s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => false)
   return _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = s)
end
