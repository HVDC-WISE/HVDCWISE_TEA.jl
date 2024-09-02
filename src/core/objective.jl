
function objective_min_fuel_cost(pm::_PM.AbstractPowerModel)

    order = _PM.calc_max_cost_index(pm.data) - 1
    cost = order <= 1 ? JuMP.AffExpr() : JuMP.QuadExpr()

    for n in _PM.nw_ids(pm)
        add_gen_operation_cost!(cost, pm, n)
        add_load_operation_cost!(cost, pm, n)
        add_slack_operation_cost!(cost, pm, n)
    end
    JuMP.@objective(pm.model, Min, cost)
end

## Auxiliary functions


function add_gen_operation_cost!(cost, pm::_PM.AbstractPowerModel, n::Int)

    for (i, gen) in _PM.ref(pm, n, :gen)
        pg = _PM.var(pm, n, :pg, i)
        if length(gen["cost"]) == 0
            cost_expr = 0
        elseif length(gen["cost"]) == 1
            cost_expr = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            cost_expr = gen["cost"][1]*pg + gen["cost"][2]
        elseif length(gen["cost"]) == 3
            cost_expr = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
        end
        JuMP.add_to_expression!(cost, cost_expr)
    end
    if get(pm.setting, "add_co2_cost", false)
        co2_emission_cost = pm.ref[:it][_PM.pm_it_sym][:co2_emission_cost]
        for (i, g) in _PM.ref(pm, n, :dgen)
            JuMP.add_to_expression!(cost, g["emission_factor"] * co2_emission_cost, _PM.var(pm, n, :pg, i))
        end
    end
    for (i, g) in _PM.ref(pm, n, :ndgen)
        JuMP.add_to_expression!(cost, g["cost_curt"], _PM.var(pm, n, :pgcurt, i))
    end
    return cost
end

function add_load_operation_cost!(cost, pm::_PM.AbstractPowerModel, n::Int)

    for (i, l) in get(_PM.ref(pm, n), :flex_load, Dict())
        # Splitting into half and half allows for better cost attribution when running single-period problems or problems with no integral constraints
        JuMP.add_to_expression!(cost, 0.5*l["cost_shift"], _PM.var(pm, n, :pshift_up, i))
        JuMP.add_to_expression!(cost, 0.5*l["cost_shift"], _PM.var(pm, n, :pshift_down, i))
        JuMP.add_to_expression!(cost, l["cost_red"], _PM.var(pm, n, :pred, i))
    end
    return cost
end
