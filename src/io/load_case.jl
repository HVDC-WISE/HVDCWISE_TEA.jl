import PowerModels as _PM
import CSV
using  DataFrames
import HiGHS
import HVDCWISE_TEA as _HWTEA
# import InfrastructureModels as _IM

include("build_raw_inputs.jl")
include("build_outputs.jl")

function load_case(case_name::String, base_mva::Int=100, save_results::Bool=false)
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

    global data = _HWTEA.parse_data(m_file, time_series) # read the data and add multinetwork dimension

    # define optimizer linear solver
    optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    # HVDCWiseTEA settings
    s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => false)

    # TEST
    _PM.propagate_topology_status!(data)  # TODO check if usefull
    #data["per_unit"] = false
    #_PM.make_mixed_units!(data) # TODO check if usefull
    #_PM.make_per_unit!(data) # TODO check if usefull

   raw_results = _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = s)
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
