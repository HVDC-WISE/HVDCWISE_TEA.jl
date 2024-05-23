using Distributed
using HVDCWISE_TEA
import PowerModels as _PM
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA))) # Root directory of HVDCWISE_TEA package

## Workers setup

nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere using HVDCWISE_TEA
@everywhere using Ipopt

## Path and solver parameters

path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case5_3grids_MC.m")
path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example_mc")
optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

## Solve the multiperiod OPF problem

run_tea(path2grid, path2data, _PM.DCPPowerModel, optimizer; setting = setting)
