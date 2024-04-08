
## Power balance

# Power balance of hybrid AC/DC multi-conductor network including storage & flexible demand
function constraint_power_balance_ac(pm::_PM.AbstractDCPModel, n::Int, i::Int, bus_arcs, bus_arcs_pst, bus_gens, bus_convs_ac, bus_loads, bus_shunts, bus_storage, gs, bs)
    p             = get(_PM.var(pm, n), :p, Dict())
    pg            = get(_PM.var(pm, n), :pg, Dict())
    ps            = get(_PM.var(pm, n), :ps, Dict())
    ppst          = get(_PM.var(pm, n), :ppst, Dict())
    pflex         = get(_PM.var(pm, n), :pflex, Dict())
    pconv_grid_ac = get(_PM.var(pm, n), :pconv_tf_fr, Dict())

    vm = 1

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(ppst[a] for a in bus_arcs_pst)
        + sum(sum(pconv_grid_ac[c][d] for d in first(axes(_PM.var(pm, n, :pconv_tf_fr, c)))) for c in bus_convs_ac)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pflex[d] for d in bus_loads)
        - sum(gs[s] for s in bus_shunts)*vm^2
    )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
    end
end
