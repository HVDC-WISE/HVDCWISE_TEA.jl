export solve_mc_acdcopf

function solve_mc_acdcopf(paths::Vector{Any}, model_type::Type, solver; kwargs...)
    data = parse_data(paths...)
    return solve_mc_acdcopf(data, model_type, solver; kwargs...)
end

function solve_mc_acdcopf(data::Dict{String,Any}, model_type::Type, solver; kwargs...)
    return _PM.solve_model(
        data, model_type, solver, build_mc_acdcopf;
        ref_extensions=[_FP.ref_add_gen!, _FP.ref_add_storage!, _CBA.ref_add_pst!, _PMMCDC.add_ref_dcgrid!, ref_add_flex_load!],
        solution_processors = [_PM.sol_data_model!],
        multinetwork = haskey(data, "dim") ? true : false,
        kwargs...
    )
end


function build_mc_acdcopf(pm::_PM.AbstractPowerModel; objective::Bool=true)

    for n in _FP.nw_ids(pm)

        # AC grid variables
        _PM.variable_bus_voltage(pm; nw = n)
        _PM.variable_branch_power(pm; nw = n)

        # AC generators
        _PM.variable_gen_power(pm; nw = n)
        _FP.expression_gen_curtailment(pm; nw = n)

        # AC phase-shift transformer
        _CBA.variable_pst(pm; nw = n)

        # AC bus slack variables
        variable_slack_power(pm; nw = n)

        # AC load
        variable_total_demand(pm; nw = n)

        if haskey(pm.data, "dim")
            # Flexible demand
            variable_flexible_demand(pm; nw = n)
            # Storage
            _PM.variable_storage_power(pm; nw = n)
            _FP.variable_absorbed_energy(pm; nw = n)
        end

        # AC/DC converter variables
        variable_mcdc_converter(pm; nw = n)

        # DC grid variables
        _PMMCDC.variable_mcdcgrid_voltage_magnitude(pm; nw = n)
        _PMMCDC.variable_mc_dcbranch_current(pm; nw = n)

    end

    # Objective
    if objective
        objective_min_fuel_cost(pm)
    end

    for n in _FP.nw_ids(pm)

        ## Steady state constraints

        # AC grid constraints
        _PM.constraint_model_voltage(pm; nw = n)

        for i in _PM.ids(pm, n, :ref_buses)
            _PM.constraint_theta_ref(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :bus)
            constraint_power_balance_ac(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branch)
            _PM.constraint_ohms_yt_from(pm, i; nw = n)
            _PM.constraint_ohms_yt_to(pm, i; nw = n)
            _PM.constraint_voltage_angle_difference(pm, i; nw = n)
            _PM.constraint_thermal_limit_from(pm, i; nw = n)
            _PM.constraint_thermal_limit_to(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :pst)
            _CBA.constraint_ohms_y_from_pst(pm, i; nw = n)
            _CBA.constraint_ohms_y_to_pst(pm, i; nw = n)
            _CBA.constraint_limits_pst(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :flex_load)
            constraint_flexible_demand(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :storage)
            _FP.constraint_storage_excl_slack(pm, i; nw = n)
            _PM.constraint_storage_thermal_limit(pm, i; nw = n)
            _PM.constraint_storage_losses(pm, i; nw = n)
        end

        # DC grid constraints
        _PMMCDC.constraint_voltage_dc(pm)
        # _PMMCDC.constraint_kcl_ground_dcgrid(pm, n)
        _PMMCDC.constraint_converter_dc_ground_shunt_ohm(pm; nw = n)

        for i in _PM.ids(pm, n, :busdc)
            _PMMCDC.constraint_kcl_shunt_dcgrid(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :branchdc)
            _PMMCDC.constraint_ohms_dc_branch(pm, i; nw = n)
        end

        for i in _PM.ids(pm, n, :convdc)
            _PMMCDC.constraint_converter_losses(pm, i; nw = n)
            _PMMCDC.constraint_converter_current(pm, i; nw = n)
            _PMMCDC.constraint_converter_dc_current(pm, i; nw = n)
            _PMMCDC.constraint_conv_transformer(pm, i; nw = n)
            _PMMCDC.constraint_conv_reactor(pm, i; nw = n)
            _PMMCDC.constraint_conv_filter(pm, i; nw = n)
            if pm.ref[:it][_PM.pm_it_sym][:nw][n][:convdc][i]["islcc"] == 1
                _PMMCDC.constraint_conv_firing_angle(pm, i; nw = n)
            end
        end

        ## Intertemporal constraints

        if haskey(pm.data, "dim")
            if _FP.is_first_id(pm, n, :hour)
                for i in _PM.ids(pm, :storage, nw = n)
                    _FP.constraint_storage_state(pm, i, nw = n)
                end

                for i in _PM.ids(pm, :storage_bounded_absorption, nw = n)
                    _FP.constraint_maximum_absorption(pm, i, nw = n)
                end

            else
                if _FP.is_last_id(pm, n, :hour)
                    for i in _PM.ids(pm, :storage, nw = n)
                        _FP.constraint_storage_state_final(pm, i, nw = n)
                    end

                    for i in _PM.ids(pm, :flex_load, nw = n)
                        _FP.constraint_shift_balance_periodic(pm, i, get(pm.setting, "demand_shifting_balance_period", 24), nw = n)
                    end
                end

                # From second hour to last hour
                prev_n = _FP.prev_id(pm, n, :hour)
                for i in _PM.ids(pm, :storage, nw = n)
                    _FP.constraint_storage_state(pm, i, prev_n, n)
                end

                for i in _PM.ids(pm, :storage_bounded_absorption, nw = n)
                    _FP.constraint_maximum_absorption(pm, i, prev_n, n)
                end
            end
        end
    end
end