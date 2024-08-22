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

## Path and solver parameters

hours = 168
path2grid = joinpath(_HWTEA_dir, "test/data", "grids/adcdc", "case39_mcdc.m")
path2data = joinpath(_HWTEA_dir, "test/data", "timeseries", "example_mc")
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);
setting_opt = Dict("solver" => "ipm", "ipm_iteration_limit" => 3000, "time_limit" => 3600.0, "output_flag" => false)
optimizer = HVDCWISE_TEA.optimizer_with_attributes(HiGHS.Optimizer, setting_opt)


## Solve the multiperiod OPF problem

run_tea(path2grid, path2data, hours, _PM.DCPPowerModel, optimizer; setting = setting)

