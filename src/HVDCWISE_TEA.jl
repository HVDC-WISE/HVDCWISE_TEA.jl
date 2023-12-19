module HVDCWISE_TEA


## Imports

import Memento
import JuMP
import InfrastructureModels as _IM
import PowerModels as _PM
import PowerModelsMCDC as _PMMCDC
import FlexPlan as _FP
import CbaOPF as _CBA


## Memento settings

# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(HVDCWISE_TEA)`
# NOTE: If this line is not included then the precompiled `HVDCWISE_TEA._LOGGER` won't be registered at runtime.
__init__() = Memento.register(_LOGGER)


## Includes

include("core/constraint_template.jl")
include("core/data.jl")
include("core/objective.jl")
include("core/ref_extension.jl")
include("core/variable.jl")

include("form/acp.jl")
include("form/dcp.jl")

include("io/parse.jl")

include("prob/opf_mc_acdc.jl")


## Exports

# HVDCWISE_TEA exports everything except internal symbols, which are defined as those whose name
# starts with an underscore. If you don't want all of these symbols in your environment,
# then use `import HVDCWISE_TEA` instead of `using HVDCWISE_TEA`.

# Do not add HVDCWISE_TEA-defined symbols to this exclude list. Instead, rename them with an
# underscore.
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all=true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_") || startswith(sym_string, "@_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
         Base.isidentifier(sym_string[2:end])))
       continue
    end
    #println("$(sym)")
    @eval export $sym
end

# The following items are also exported for user-friendlyness when calling `using HVDCWISE_TEA`,
# so that users do not need to import JuMP to use a solver with HVDCWISE_TEA.
import JuMP: optimizer_with_attributes
export optimizer_with_attributes

import JuMP: TerminationStatusCode
export TerminationStatusCode

import JuMP: ResultStatusCode
export ResultStatusCode

for status_code_enum in [TerminationStatusCode, ResultStatusCode]
    for status_code in instances(status_code_enum)
        @eval import JuMP: $(Symbol(status_code))
        @eval export $(Symbol(status_code))
    end
end


end
