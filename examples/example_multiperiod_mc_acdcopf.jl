import PowerModels as _PM
import HVDCWISE_TEA as _HWTEA
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

## Input parameters

path2grid = joinpath(_HWTEA_dir, "test/data/grids/acdc/case5_3grids_MC.m")
path2data = joinpath(_HWTEA_dir, "test/data/timeseries/example")
data = _HWTEA.parse_data(path2grid, path2data)

## Solve the multiperiod OPF problem

optimizer = _HWTEA.optimizer_with_attributes(Ipopt.Optimizer)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>true), "conv_losses_mp" => false)
results = _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = setting)
