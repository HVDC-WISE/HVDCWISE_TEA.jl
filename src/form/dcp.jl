
## Variables

"variables for bus power slack"
function variable_slack_power(pm::_PM.AbstractDCPModel; kwargs...)
    variable_slack_power_real(pm; kwargs...)
end

"variables for AC/DC converters"
function variable_mcdc_converter(pm::_PM.AbstractDCPModel; kwargs...)
    variable_conv_tranformer_flow(pm; kwargs...)
    variable_conv_reactor_flow(pm; kwargs...)

    _PMMCDC.variable_converter_active_power(pm; kwargs...)
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
end

function variable_conv_tranformer_flow(pm::_PM.AbstractDCPModel; kwargs...)
    _PMMCDC.variable_conv_transformer_active_power_to(pm; kwargs...)
end

function variable_conv_reactor_flow(pm::_PM.AbstractDCPModel; kwargs...)
    _PMMCDC.variable_conv_reactor_active_power_from(pm; kwargs...)
end

## Constraints

# Power balance of hybrid AC/DC multi-conductor network including storage & flexible demand
function constraint_power_balance_ac(pm::_PM.AbstractDCPModel, n::Int, i::Int, bus_arcs, bus_gens, bus_convs_ac, bus_loads, bus_shunts, bus_storage, gs, bs)
    p_slack_up    = get(_PM.var(pm, n), :p_slack_up, Dict(i => 0.0))[i]
    p_slack_down  = get(_PM.var(pm, n), :p_slack_down, Dict(i => 0.0))[i]

    p             = get(_PM.var(pm, n), :p, Dict())
    pg            = get(_PM.var(pm, n), :pg, Dict())
    ps            = get(_PM.var(pm, n), :ps, Dict())
    pflex         = get(_PM.var(pm, n), :pflex, Dict())
    pconv_grid_ac = get(_PM.var(pm, n), :pconv_tf_fr, Dict())

    vm = 1

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(sum(pconv_grid_ac[c][d] for d in first(axes(_PM.var(pm, n, :pconv_tf_fr, c)))) for c in bus_convs_ac)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pflex[d] for d in bus_loads)
        - sum(gs[s] for s in bus_shunts)*vm^2
        + p_slack_up
        - p_slack_down
    )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
    end
end

## Expressions


## Expressions

"expression for nodal slack power cost"
function add_slack_operation_cost!(cost, pm::_PM.AbstractDCPModel, n::Int)

    value = _PM.ref(pm, n, :baseMVA) * 50000
    for (i, bus) in _PM.ref(pm, n, :bus)
        p_slack_up   = get(_PM.var(pm, n), :p_slack_up, Dict(i => 0.0))[i]
        p_slack_down = get(_PM.var(pm, n), :p_slack_down, Dict(i => 0.0))[i]

        JuMP.add_to_expression!(cost, value, (p_slack_up + p_slack_down))
    end
    return cost
end