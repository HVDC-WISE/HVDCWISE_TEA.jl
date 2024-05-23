import HVDCWISE_TEA as _HWTEA

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package
include("load_case.jl")

function run_study(dir_name::String, case_name::String)
    work_dir = joinpath(_HWTEA_dir, "studies\\$case_name\\$dir_name")
    baseMVA = 100
    load_case(work_dir, work_dir, case_name, baseMVA, true)
end

function run_study_2(dir_name::String, case_name::String)
    work_dir = joinpath(_HWTEA_dir, "studies\\$case_name\\$dir_name")
    base_mva = 100

    m_file = joinpath(work_dir, "$case_name.m")
    time_series = build_time_series(work_dir, base_mva)
    data = _HWTEA.parse_data(m_file, time_series)

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
    build_outputs(work_dir, case_name, raw_results)
end

@time run_study_2("draft", "IEEE39")
# @time run_study("3_hours_without_storage", "IEEE39")
# run_study("1_day_without_storage", "IEEE39")
# run_study("1_day_with_storage", "IEEE39")
# @time run_study("1_day_with_storage", "IEEE39")
# @time run_study("1_week_with_storage", "IEEE39")
# @time run_study("1_year_with_storage", "IEEE39")