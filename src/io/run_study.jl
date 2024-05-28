using Distributed
using HVDCWISE_TEA
import HVDCWISE_TEA as _HWTEA
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA))) # Root directory of HVDCWISE_TEA package


_PM.silence()

function run_study(work_dir::String, case_name::String)
    base_mva = 100
    # build_raw_inputs(work_dir) -> the .m at the root and 1 folder per micro with the .csv
    build_raw_inputs(work_dir, work_dir, case_name, base_mva, true)
    #
    # run_tea

    path2grid = joinpath(work_dir, "$case_name.m")
    path2data = joinpath(work_dir, case_name)

    optimizer = _HWTEA.optimizer_with_attributes(Ipopt.Optimizer)
    setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

    run_tea(path2grid, path2data, _PM.DCPPowerModel, optimizer; setting = setting)  # FIXME csv in pu or in MW ? # FIXME uncomment

    # build_outputs_from_csv
    results_dir = joinpath(work_dir, case_name, "results")
    build_outputs_from_csv(results_dir, base_mva, true) # FIXME uncomment
    #
end

# work_dir = joinpath(_HWTEA_dir, "studies\\IEEE39\\3_hours_with_storage")
# case_name = "IEEE39"
# @time run_study(work_dir, case_name)

# work_dir = joinpath(_HWTEA_dir, "studies\\IEEE39\\draft")
# case_name = "test05b"
# @time load_case(work_dir, work_dir, case_name, 100, true)

@time run_study(joinpath(_HWTEA_dir, "studies\\IEEE39\\3_hours_without_storage"), "IEEE39")

# @time run_study_2("draft", "IEEE39")
# @time run_study("3_hours_without_storage", "IEEE39")
# run_study("1_day_without_storage", "IEEE39")
# run_study("1_day_with_storage", "IEEE39")
# @time run_study("1_day_with_storage", "IEEE39")
# @time run_study("1_week_with_storage", "IEEE39")
# @time run_study("1_year_with_storage", "IEEE39")