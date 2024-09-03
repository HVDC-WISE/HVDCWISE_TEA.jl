using Distributed
using HVDCWISE_TEA
import PowerModels as _PM
import HiGHS

const _HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA)));

## Workers setup

nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using HiGHS
    HVDCWISE_TEA.silence()
end

## Solver parameters

setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);
setting_opt = Dict("presolve" => "on", "solver" => "ipm", "run_crossover" => "off", "ipm_iteration_limit" => 3000, "time_limit" => 3600.0, "output_flag" => false)
optimizer = HVDCWISE_TEA.optimizer_with_attributes(HiGHS.Optimizer, setting_opt...)


## Input files

hours_per_subsimulation = 168  # The yearly problem is split into subproblems of this size (1 week is 168h)
# path2grid = joinpath(_HWTEA_dir, "test/data", "grids/adcdc", "case39_mcdc.m")
# path2data = joinpath(_HWTEA_dir, "test/data", "timeseries", "example_mc")
work_dir = joinpath(_HWTEA_dir, "studies\\simple_use_case")
# work_dir = joinpath(_HWTEA_dir, "studies\\2024-08-23 case39")

## Solve the multiperiod OPF problem

# run_tea(path2grid, path2data, hours_per_subsimulation, _PM.DCPPowerModel, optimizer; setting = setting)

# build_grid_model(work_dir, 100)
# gather_opf_results(work_dir, "Macro", 100)
# run_simulation(work_dir, hours_per_subsimulation, 100)

@time run_study(work_dir, hours_per_subsimulation)