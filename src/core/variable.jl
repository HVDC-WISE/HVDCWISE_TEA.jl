
"variables for flexible load"
function variable_flexible_demand(pm::_PM.AbstractPowerModel; kwargs...)

    _FP.variable_total_flex_demand(pm; kwargs...)
    _FP.variable_demand_curtailment(pm; kwargs...)

    if haskey(pm.data, "dim")
        _FP.variable_demand_reduction(pm; kwargs...)
        _FP.variable_demand_shifting_upwards(pm; kwargs...)
        _FP.variable_demand_shifting_downwards(pm; kwargs...)
    end
end

"variable for bus active power slack"
function variable_slack_power_real(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    var1 = _PM.var(pm, nw)[:p_slack_up] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_p_slack_up",
        lower_bound = 0,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "p_slack_up_start", 0.0)
    )
    var2 = _PM.var(pm, nw)[:p_slack_down] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_p_slack_down",
        lower_bound = 0,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "p_slack_down_start", 0.0)
    )
    report && _PM.sol_component_value(pm, nw, :bus, :p_slack_up, _PM.ids(pm, nw, :bus), var1)
    report && _PM.sol_component_value(pm, nw, :bus, :p_slack_down, _PM.ids(pm, nw, :bus), var2)
end

"variable for bus reactive power slack"
function variable_slack_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    var1 = _PM.var(pm, nw)[:q_slack_up] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_q_slack_up",
        lower_bound = 0,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "q_slack_up_start", 0.0)
    )
    var2 = _PM.var(pm, nw)[:q_slack_down] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_q_slack_down",
        lower_bound = 0,
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "q_slack_down_start", 0.0)
    )
    report && _PM.sol_component_value(pm, nw, :bus, :q_slack_up, _PM.ids(pm, nw, :bus), var1)
    report && _PM.sol_component_value(pm, nw, :bus, :q_slack_down, _PM.ids(pm, nw, :bus), var2)
end