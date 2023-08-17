import PowerModels as _PM
import HVDCWISE_TEA as _HWTEA
import HiGHS
import Memento
using Test


# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(_PM), "error")
Memento.setlevel!(Memento.getlogger(_HWTEA), "error")

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

##

lp_optimizer = _HWTEA.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
s = Dict("output" => Dict("branch_flows" => true, "duals" =>true), "conv_losses_mp" => false)


@testset "HVDCWISE_TEA.jl" begin

    # Network component flexibility
    include("flexibility.jl")

    # Problem specification
    include("prob.jl")

    # Exported symbols
    include("export.jl")

end;
