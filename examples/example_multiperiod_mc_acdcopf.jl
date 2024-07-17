import PowerModels as _PM
import HVDCWISE_TEA as _HWTEA
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

## Path and solver parameters

path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case39_mcdc.m")
path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example_39")
optimizer = _HWTEA.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false)

## Solve the multiperiod OPF problem

results = _HWTEA.solve_mc_acdcopf([path2grid, [path2data]], _PM.DCPPowerModel, optimizer; setting = setting)
