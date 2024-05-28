import HVDCWISE_TEA as _HWTEA
import HiGHS
import Ipopt
#=
import FlexPlan as _FP
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC
import CSV
using  DataFrames

import HVDCWISE_TEA as _HWTEA
import JSON

include("build_raw_inputs.jl")
include("build_outputs.jl")
include("parse.jl")
=#

function load_case(case_name::String, base_mva::Int=100, save_results::Bool=false)
    _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package
    input_dir = joinpath(_HWTEA_dir, "test\\data\\$case_name")
    work_dir = joinpath(_HWTEA_dir, "output\\data\\$case_name")
    raw_results = load_case(input_dir, work_dir, case_name, base_mva, save_results)
    return raw_results
end

function load_case(input_dir::String, work_dir::String, case_name::String, base_mva::Int=100, save_results::Bool=false)

    ## Buil .m & .csv files
    build_raw_inputs(input_dir, work_dir, case_name, base_mva)

    ## Define the path of the .m file
    m_file = joinpath(work_dir, "$case_name.m")

    ## Read CSV files and save data in the time_series dictionary
    time_series = build_time_series(work_dir, base_mva)

    # parse data
    sn_data = parse_file(m_file)  # Parse initial single-network structure
    hours = length(find_value(time_series))  # Find number of simulated hours (a single dimension is used)
    _FP.add_dimension!(sn_data, :hour, hours)  # Add hour dimension to single-network data
    mn_data = _FP.make_multinetwork(sn_data, time_series; share_data=false)  # Create the multinetwork data dictionary
    _PMMCDC.make_multiconductor!(mn_data)  # Convert the DC grid (if any) of each single-network model to multi-conductor

    # global data = _HWTEA.parse_data(m_file, time_series) # read the data and add multinetwork dimension # deprecated function

    if save_results
        json_path = joinpath(work_dir, "$case_name"*"_inputs.json")
        stringdata = JSON.json(mn_data, 4)
        open(json_path, "w") do f
            write(f, stringdata)
        end
    end

    # define optimizer linear solver
    # optimizer = _HWTEA.optimizer_with_attributes(Ipopt.Optimizer)
    optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    # HVDCWiseTEA settings
    s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => false)

    # TEST
    _PM.propagate_topology_status!(mn_data)  # TODO check if usefull
    #mn_data["per_unit"] = false
    #_PM.make_mixed_units!(mn_data) # TODO check if usefull
    #_PM.make_per_unit!(mn_data) # TODO check if usefull

   raw_results = _HWTEA.solve_mc_acdcopf(mn_data, _PM.DCPPowerModel, optimizer; setting = s)
   if save_results
    build_outputs(work_dir, case_name, raw_results)
   end
   return raw_results
end

function build_time_series(work_dir::String, base_mva::Int=100)

    ## Read CSV files and save data in the time_series dictionary
    time_series = Dict{String, Any}()
    for component_type in ["branch", "branchdc", "convdc", "gen", "load", "storage"]
        comp_dir = joinpath(work_dir, "$component_type")
        if isdir(comp_dir) && isempty(readdir(comp_dir)) == false
            time_series["$component_type"] = Dict{String, Any}()
            for file_name in readdir(comp_dir)
                file_path = joinpath(comp_dir, "$file_name")
                @assert isfile(file_path) "$file_path is not a file"
                file_data = CSV.File(file_path, delim=',') |> DataFrames.DataFrame
                attribute_name = file_name[1:length(file_name)-4]
                if occursin("status", attribute_name)
                    base_value = 1
                else  # All attributes except status are provided in per unit
                    base_value = base_mva
                end
                if length(time_series["$component_type"]) == 0
                    for component_id in names(file_data)
                        time_series["$component_type"]["$component_id"] = Dict{String, Any}()
                    end
                end
                for component_id in names(file_data)
                    time_series["$component_type"]["$component_id"]["$attribute_name"] = float(file_data[:, component_id]) / base_value
                end
            end
        end
    end

   return time_series
end
