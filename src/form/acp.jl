
## Power balance

# Power balance of hybrid AC/DC multi-conductor network including storage & flexible demand
function constraint_power_balance_ac(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_pst, bus_gens, bus_convs_ac, bus_loads, bus_shunts, bus_storage, gs, bs)
    vm            = _PM.var(pm, n,  :vm, i)
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
    )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
