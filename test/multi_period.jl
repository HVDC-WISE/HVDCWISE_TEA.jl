
import FlexPlan as _FP
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC
import HVDCWISE_TEA as _HWTEA
import HiGHS

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

## Time series
hours = 24
time_series = Dict{String, Any}(
    "gen" => Dict(
        "1" => Dict("gen_status" => fill(1, hours)),
        "2" => Dict("gen_status" => fill(1, hours)),
        "3" => Dict("gen_status" => fill(1, hours)),
        "4" => Dict("gen_status" => fill(1, hours)),
        "5" => Dict("gen_status" => fill(1, hours))
    ),
    "load" => Dict(
        "1" => Dict("pd" => fill(20.0, hours) ./ 100),
        "2" => Dict("pd" => fill(45.0, hours) ./ 100),
        "3" => Dict("pd" => fill(40.0, hours) ./ 100),
        "4" => Dict("pd" => fill(60.0, hours) ./ 100),
        "5" => Dict("pd" => fill(20.0, hours) ./ 100),
        "6" => Dict("pd" => fill(45.0, hours) ./ 100),
        "7" => Dict("pd" => fill(40.0, hours) ./ 100),
        "8" => Dict("pd" => fill(60.0, hours) ./ 100),
        "9" => Dict("pd" => reduce(vcat, [fill(3.0, 5), fill(15, 2), fill(70, 3), fill(50, 6), fill(70, 2), fill(400, 4), fill(25, 2)]) ./ 100)
    )
)

## Input parameters

path = joinpath(_HWTEA_dir, "test/data/case5/case5_3grids_MC.m")
data = _HWTEA.parse_data(path, time_series)

## Solve the multiperiod OPF problem

optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
s = Dict("output" => Dict("branch_flows" => true, "duals" =>true), "conv_losses_mp" => false)

results = _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, optimizer; setting = s)


