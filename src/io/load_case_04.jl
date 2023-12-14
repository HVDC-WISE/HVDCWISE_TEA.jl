## Import relevant packages
import FlexPlan as _FP
import PowerModels as _PM
import CSV
import DataFrames
import HiGHS


nlp_optimizer = _FP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)

## Load test case
file = normpath(@__DIR__,"../../test/data/Test04/.m_files/test04.m")
data = _FP.parse_file(file) # Parse input file to obtain data dictionary

# Read CSV files
# get all rows of the 2nd column, the first row is the timestamp & the index of the buses to which the load is connected
demand = CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/load_fix_MW.csv"),DataFrames.DataFrame)[:,2:end]
#TODO a vérifier
demand_pu = demand[:, 1] ./maximum(demand[:, 1])

# get all rows of the 2nd column, the first row is the timestamp & the index of the generators in the mpc.generator table in the .m file
generation_status = CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/status_generators.csv"),DataFrames.DataFrame)[:,2:end]

# get all rows of the 2nd column, the first row is the timestamp & the index of the AC line in the mpc.generator table in the .m file
AC_line_status = CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/status_lines_AC.csv"),DataFrames.DataFrame)[:,2:end]

# get all rows of the 2nd column, the first row is the timestamp & the index of the generators in the mpc.generator table in the .m file
DC_line_status = CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/status_lines_DC_temp.csv"),DataFrames.DataFrame)[:,2:end]

# get all rows of the 2nd column, the first row is the timestamp & the index of the generators in the mpc.generator table in the .m file
converters_status = CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/status_converters_temp.csv"),DataFrames.DataFrame)[:,2:end]

#Get the number of time points
number_of_hours = size(CSV.read(normpath(@__DIR__,"../../test/data/Test04/.csv_files/load_fix_MW.csv"),DataFrames.DataFrame),1)

#=
Following dimensions:
    hour: the finest time granularity that can be represented in a model. During an hour, each continuous variable has a constant value.
    year: an investment period. Different investment decisions can be made in different years.
    scenario: one of the different possible sets of values related to renewable generation and consumption data.
=#
planning_horizon = 10 # Years to scale generation costs
_FP.add_dimension!(data, :hour, number_of_hours) # Add dimension, e.g. hours
_FP.add_dimension!(data, :scenario, Dict(1 => Dict{String,Any}("probability"=>1)))
_FP.add_dimension!(data, :year, 1; metadata = Dict{String,Any}("scale_factor"=>planning_horizon))

_FP.scale_data!(data) # Scale investment & operational cost data based on planning years & hours


# Initialize genprofile arrays
genprofile = ones(length(data["gen"]), number_of_hours)
# Populate genprofile based on generation status CSV file
genprofile[:, :] .= generation_status[:, 1]'

# Initialize loadprofile arrays
loadprofile = ones(length(data["load"]), number_of_hours)
# Populate loadprofile based on demand CSV file
loadprofile[:, :] .= demand_pu'

# Create time series data to be passed to the data dictionay
time_series = _FP.make_time_series(data, number_of_hours; loadprofile = permutedims(loadprofile), genprofile = permutedims(genprofile))
mn_data = _FP.make_multinetwork(data, time_series) # Create the multinetwork data dictionary

#Modifying the AC line status with the timeseries data
for (col_idx, col) in enumerate(names(AC_line_status))
    for (row_idx, value) in enumerate(AC_line_status[!, col])
        timestamp = row_idx
        branch_id = col
        br_status = value
        mn_data["nw"]["$timestamp"]["branch"]["$branch_id"]["br_status"] = br_status
    end
end

# TODO replace mn_data["nw"]["timestamp"]["branchdc"]["branch_id"]["status"] : mn_data["nw"]["1"]["branchdc"]["1"]["status"] = 0 or 1
# TODO once PR https://github.com/Electa-Git/PowerModelsMCDC.jl/pull/17/ is validated, include all line status (depending on the DC network topology)
# TODO pole + : mn_data["nw"]["timestamp"]["branchdc"]["branch_id"]["positive_pole_status"] : mn_data["nw"]["1"]["branchdc"]["1"]["status_p"] = 0 or 1
# TODO pole - : mn_data["nw"]["timestamp"]["branchdc"]["branch_id"]["negative_pole_status"] : mn_data["nw"]["1"]["branchdc"]["1"]["status_n"] = 0 or 1
# TODO pole DMR : mn_data["nw"]["timestamp"]["branchdc"]["branch_id"]["MR_pole_status"] : mn_data["nw"]["1"]["branchdc"]["1"]["status_r"] = 0 or 1
#Modifying the DC line status with the timeseries data
for (col_idx, col) in enumerate(names(DC_line_status))
    for (row_idx, value) in enumerate(DC_line_status[!, col])
        timestamp = row_idx
        branchdc_id = col
        status = value
        mn_data["nw"]["$timestamp"]["branchdc"]["$branchdc_id"]["status"] = status
    end
end

# TODO replace mn_data["nw"]["timestamp"]["convdc"]["conv_id"]["conv_pole_status"] : mn_data["nw"]["1"]["convdc"]["1"]["status"] = 0 or 1
# TODO once PR https://github.com/Electa-Git/PowerModelsMCDC.jl/pull/17/ is validated, include all converter status (depending on the DC network topology)
# TODO replace mn_data["nw"]["timestamp"]["convdc"]["conv_id"]["positive_pole_status"] : mn_data["nw"]["1"]["convdc"]["1"]["status_p"] = 0 or 1
# TODO replace mn_data["nw"]["timestamp"]["convdc"]["conv_id"]["negative_pole_status"] : mn_data["nw"]["1"]["convdc"]["1"]["status_n"] = 0 or 1
#Modifying the converter status with the timeseries data
for (col_idx, col) in enumerate(names(converters_status))
    for (row_idx, value) in enumerate(converters_status[!, col])
        timestamp = row_idx
        conv_id = col
        status = value
        mn_data["nw"]["$timestamp"]["convdc"]["$conv_id"]["status"] = status
    end
end


# FlexPlan settings
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "add_co2_cost" => false)
results = _FP.stoch_flex_tnep(mn_data, _PM.DCPPowerModel, nlp_optimizer; setting = s)
results_2 = _FP.simple_stoch_flex_tnep(mn_data,_PM.DCPPowerModel, nlp_optimizer; setting = s)

println("----------------STOCHASTIC FLEX TNEP----------------")
_PM.print_summary(results["solution"]["nw"]["1"])
_PM.print_summary(results["solution"]["nw"]["2"])
_PM.print_summary(results["solution"]["nw"]["3"])
_PM.print_summary(results["solution"]["nw"]["4"])
_PM.print_summary(results["solution"]["nw"]["5"])
_PM.print_summary(results["solution"]["nw"]["6"])
println("----------------SIMPLE STOCHASTIC FLEX TNEP----------------")
_PM.print_summary(results_2["solution"]["nw"]["1"])
_PM.print_summary(results_2["solution"]["nw"]["2"])
_PM.print_summary(results_2["solution"]["nw"]["3"])
_PM.print_summary(results_2["solution"]["nw"]["4"])
_PM.print_summary(results_2["solution"]["nw"]["5"])
_PM.print_summary(results_2["solution"]["nw"]["6"])