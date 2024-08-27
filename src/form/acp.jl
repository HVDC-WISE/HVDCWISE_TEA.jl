
## variables

"variables for active and reactive bus slack"
function variable_slack_power(pm::_PM.AbstractACPModel; kwargs...)
    variable_slack_power_real(pm; kwargs...)
    variable_slack_power_imaginary(pm; kwargs...)
end

"variables for AC/DC converters"
function variable_mcdc_converter(pm::_PM.AbstractACPModel; kwargs...)
    _PMMCDC.variable_conv_tranformer_flow(pm; kwargs...)
    _PMMCDC.variable_conv_reactor_flow(pm; kwargs...)

    _PMMCDC.variable_converter_active_power(pm; kwargs...)
    _PMMCDC.variable_converter_reactive_power(pm; kwargs...)
    _PMMCDC.variable_acside_current(pm; kwargs...)
    _PMMCDC.variable_dcside_current(pm; kwargs...)
    _PMMCDC.variable_dcside_current_ground(pm; kwargs...)
    _PMMCDC.variable_dcside_current_grounding_shunt(pm; kwargs...)
    _PMMCDC.variable_dcside_power(pm; kwargs...)
    _PMMCDC.variable_dcside_ground_power(pm; kwargs...)
    _PMMCDC.variable_dcside_grounding_shunt_power(pm; kwargs...)
    _PMMCDC.variable_converter_firing_angle(pm; kwargs...)

    _PMMCDC.variable_converter_filter_voltage(pm; kwargs...)
    _PMMCDC.variable_converter_internal_voltage(pm; kwargs...)

    _PMMCDC.variable_converter_to_grid_active_power(pm; kwargs...)
    _PMMCDC.variable_converter_to_grid_reactive_power(pm; kwargs...)
end

## Constraints

# Power balance of hybrid AC/DC multi-conductor network including storage & flexible demand
function constraint_power_balance_ac(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_pst, bus_gens, bus_convs_ac, bus_loads, bus_shunts, bus_storage, gs, bs)
    vm            = _PM.var(pm, n,  :vm, i)
    p_slack_up    = _PM.var(pm, n, :p_slack_up, i)
    q_slack_up    = _PM.var(pm, n, :q_slack_up, i)
    p_slack_down  = _PM.var(pm, n, :p_slack_down, i)
    q_slack_down  = _PM.var(pm, n, :q_slack_down, i)

    p             = get(_PM.var(pm, n), :p, Dict())
    q             = get(_PM.var(pm, n), :q, Dict())
    pg            = get(_PM.var(pm, n), :pg, Dict())
    qg            = get(_PM.var(pm, n), :qg, Dict())
    ps            = get(_PM.var(pm, n), :ps, Dict())
    qs            = get(_PM.var(pm, n), :qs, Dict())
    ppst          = get(_PM.var(pm, n), :ppst, Dict())
    qpst          = get(_PM.var(pm, n), :qpst, Dict())
    pflex         = get(_PM.var(pm, n), :pflex, Dict())
    qflex         = get(_PM.var(pm, n), :qflex, Dict())
    pconv_grid_ac = get(_PM.var(pm, n), :pconv_tf_fr, Dict())
    qconv_grid_ac = get(_PM.var(pm, n), :qconv_tf_fr, Dict())

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(ppst[a] for a in bus_arcs_pst)
        + sum(sum(pconv_grid_ac[c][d] for d in 1:length(_PM.var(pm, n, :pconv_tf_fr, c))) for c in bus_convs_ac)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pflex[d] for d in bus_loads)
        - sum(gs[s] for s in bus_shunts)*vm^2
        + p_slack_up
        - p_slack_down
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(qpst[a] for a in bus_arcs_pst)
        + sum(sum(qconv_grid_ac[c][d] for d in 1:length(_PM.var(pm, n, :qconv_tf_fr, c))) for c in bus_convs_ac)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qflex[d] for d in bus_loads)
        - sum(bs[s] for s in bus_shunts)*vm^2
        + q_slack_up
        - q_slack_down
    )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

## Expressions

"expression for nodal slack power cost"
function add_slack_operation_cost!(cost, pm::_PM.AbstractACPModel, n::Int)

    baseMVA = _PM.ref(pm, n, :baseMVA)
    for (i, bus) in _PM.ref(pm, n, :bus)
        p_slack_up   = _PM.var(pm, n, :p_slack_up, i)
        q_slack_up   = _PM.var(pm, n, :q_slack_up, i)
        p_slack_down = _PM.var(pm, n, :p_slack_down, i)
        q_slack_down = _PM.var(pm, n, :q_slack_down, i)
        
        JuMP.add_to_expression!(cost, 50000 * baseMVA, (p_slack_up + p_slack_down + q_slack_up + q_slack_down))
    end
    return cost
end