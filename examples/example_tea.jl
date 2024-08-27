using Distributed
using HVDCWISE_TEA
import PowerModels as _PM
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA)));

## Workers setup

nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using Ipopt
    HVDCWISE_TEA.silence()
end

## Solver parameters
optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

## Input files

hours_per_subsimulation = 96  # The yearly problem is split into subproblems of this size (1 week is 168h)
# path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case39_mcdc.m")
# path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example_mc")
work_dir = joinpath(_HWTEA_dir, "studies\\simple_use_case")
# work_dir = joinpath(_HWTEA_dir, "studies\\2024-08-23 case39")

## Solve the multiperiod OPF problem

# run_tea(path2grid, path2data, hours_per_subsimulation, _PM.DCPPowerModel, optimizer; setting = setting)
@time run_study(work_dir, hours_per_subsimulation)
