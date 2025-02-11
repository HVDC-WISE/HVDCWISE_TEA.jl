
"Variable for the actual active load demand at each load accounting for market flexibility"
function variable_total_demand(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    var = _PM.var(pm, nw)[:pflex] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_pflex",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "pd")
    )

    if bounded
        limit = 1.0E-3 / _PM.ref(pm, nw, :baseMVA)
        for (i, load) in _PM.ref(pm, nw, :load)
            if load["pd"] <= limit || in(i, _PM.ids(pm, nw, :fixed_load))
                JuMP.fix(var[i], load["pd"])
            else
                JuMP.set_lower_bound(var[i], 0.0)
                JuMP.set_upper_bound(var[i], load["pd"] * (1 + load["pshift_up_rel_max"]))
            end
        end
    end
    report && _PM.sol_component_value(pm, nw, :load, :pflex, _PM.ids(pm, nw, :load), var)
end

"Variable for power not consumed (voluntary load reduction) at each flexible load"
function variable_demand_reduction(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    var = _PM.var(pm, nw)[:pred] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :flex_load)], base_name="$(nw)_pred", start = 0.0
    )

    if bounded
        limit = 1.0E-3 / _PM.ref(pm, nw, :baseMVA)
        for (i, load) in _PM.ref(pm, nw, :flex_load)
            if load["pd"] <= limit
                JuMP.fix(var[i], 0.0)
            else
                JuMP.set_lower_bound(var[i], 0.0)
                JuMP.set_upper_bound(var[i], load["pd"] * load["pred_rel_max"])
            end
        end
    end

    if report
        _PM.sol_component_value(pm, nw, :load, :pred, _PM.ids(pm, nw, :flex_load), var)
        _PM.sol_component_fixed(pm, nw, :load, :pred, _PM.ids(pm, nw, :fixed_load), 0.0)
    end
end

"Variable for upward demand shifting at each flexible load"
function variable_demand_shifting_upwards(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    var = _PM.var(pm, nw)[:pshift_up] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :flex_load)], base_name="$(nw)_pshift_up", start = 0
    )

    if bounded
        limit = 1.0E-3 / _PM.ref(pm, nw, :baseMVA)
        for (i, load) in _PM.ref(pm, nw, :flex_load)
            if load["pd"] <= limit
                JuMP.fix(var[i], 0.0)
            else
                JuMP.set_lower_bound(var[i], 0.0)
                JuMP.set_upper_bound(var[i], load["pd"] * load["pshift_up_rel_max"])
            end
        end
    end

    if report
        _PM.sol_component_value(pm, nw, :load, :pshift_up, _PM.ids(pm, nw, :flex_load), var)
        _PM.sol_component_fixed(pm, nw, :load, :pshift_up, _PM.ids(pm, nw, :fixed_load), 0.0)
    end
end

"Variable for downward demand shifting at each flexible load"
function variable_demand_shifting_downwards(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    var = _PM.var(pm, nw)[:pshift_down] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :flex_load)], base_name="$(nw)_pshift_down", start = 0
    )

    if bounded
        limit = 1.0E-3 / _PM.ref(pm, nw, :baseMVA)
        for (i, load) in _PM.ref(pm, nw, :flex_load)
            if load["pd"] <= limit
                JuMP.fix(var[i], 0.0)
            else
                JuMP.set_lower_bound(var[i], 0.0)
                JuMP.set_upper_bound(var[i], load["pd"] * load["pshift_down_rel_max"])
            end
        end
    end
    if report
        _PM.sol_component_value(pm, nw, :load, :pshift_down, _PM.ids(pm, nw, :flex_load), var)
        _PM.sol_component_fixed(pm, nw, :load, :pshift_down, _PM.ids(pm, nw, :fixed_load), 0.0)
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