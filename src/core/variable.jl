# Additional variables modified or not used in PowerModels, PowerModelsMCDC, FlexPlan

## Variables for flexible demand

function variable_flexible_demand(pm::_PM.AbstractPowerModel; kwargs...)

    _FP.variable_total_flex_demand(pm; kwargs...)
    _FP.variable_demand_curtailment(pm; kwargs...)

    if haskey(pm.data, "dim")
        _FP.variable_demand_reduction(pm; kwargs...)
        _FP.variable_demand_shifting_upwards(pm; kwargs...)
        _FP.variable_demand_shifting_downwards(pm; kwargs...)

        _FP.variable_energy_not_consumed(pm; kwargs...)
        _FP.variable_total_demand_shifting_upwards(pm; kwargs...)
        _FP.variable_total_demand_shifting_downwards(pm; kwargs...)
    end
end
