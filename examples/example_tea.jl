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

## Path and solver parameters

hours = 96
path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case39_mcdc.m")
path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example_mc")
optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

## Solve the multiperiod OPF problem

run_tea(path2grid, path2data, hours, _PM.DCPPowerModel, optimizer; setting = setting)
