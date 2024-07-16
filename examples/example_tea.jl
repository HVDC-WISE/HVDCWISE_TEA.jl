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

path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case39_mcdc.m")
path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example_mc")
# path2grid = joinpath(_HWTEA_dir, "studies/39bus/case39_mcdc.m")
# path2data = joinpath(_HWTEA_dir, "studies/39bus/1_week_without_storage")

optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

## Solve the multiperiod OPF problem

run_tea(path2grid, path2data, _PM.DCPPowerModel, optimizer; setting = setting)

# Gather results in an Excel file

macro_results_dir = joinpath(path2data, "results")
build_outputs_from_csv(macro_results_dir, 100, true)  # base_mva=100 in the .m file