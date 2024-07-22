"""
    process_additional_data!(data::Dict{String,Any})

Organise, process and perform checks on data tables that are not used/recognised in PowerModels.
These may include non-dispatchable generators, flexible loads, storages, phase-shift transformers and
DC components (i.e., AC/DC converters and DC branches).

Optional tables are:
 - `gencost` and `storage`, which are processed in PowerModels;
 - `ndgen`, `storage_extra`, `load_extra`, which are processed in FlexPlan;
 - `busdc`, `convdc` and `branchdc`, which are processed in PowerModelsMCDC;
 - `pst`, which is processed in CbaOPF.

Other tables may be added as well: they will be made available in the returned object, without any check.

## Extended help

The available `load_extra` table parameters are:
 - `flex` : boolean thta defines if load is flexible or not. If not flexible, only involutary curtailment is active.
 - `pf_angle` : fixed power factor angle in radians, giving the reactive power as `Q = P ⨉ tan(θ)`
 - `cost_curt` : compensation for involuntary load curtailment (i.e. involuntary demand reduction) [€/MWh];
 - `cost_red` : compensation for consuming less (i.e. voluntary demand reduction) [€/MWh];
 - `cost_shift` : Compensation for demand shifting, applied half to the power shifted upward and half to the power shifted downward [€/MWh]
 - `co2_cost` : CO2 costs for enabling flexible demand [€]
 - `pred_rel_max` : superior bound on voluntary load reduction (not consumed power) as a fraction of the total reference demand [p.u.];
 - `pshift_up_rel_max` : superior bound on upward demand shifted as a fraction of the total reference demand [p.u.];
 - `pshift_down_rel_max`: superior bound on downward demand shifted as a fraction of the total reference demand [p.u.];
 - `ered_rel_max` : superior bound on voluntary energy reduction as a fraction of the total reference demand [p.u.];
 - `eshift_rel_max` : superior bound on shifted energy as a fraction of the total reference demand [p.u.];
 - `tshift_up` : Recovery period for upward demand shifting [h];
 - `tshift_down` : Recovery period for downward demand shifting [h];

"""
function process_additional_data!(data::Dict{String, Any})
    
    _FP._add_gen_data!(data)
    _add_dcgrid_data!(data)
    _add_storage_data!(data)
    _add_flexible_demand_data!(data)
    _add_pst_data!(data)
end


## Auxiliary function


function _add_dcgrid_data!(data::Dict{String, Any})
    if haskey(data, "busdc")
        _PMMCDC.process_additional_data!(data)
    end
end

function _add_flexible_demand_data!(data::Dict{String, Any})

    rescale_cost = x -> x*data["baseMVA"]

    if !haskey(data, "load_extra")
        for (i, load) in data["load"]
            # Whether load is flexible (boolean)
            load["flex"] = 0
            # Compensation for load curtailment (i.e. involuntary demand reduction) (€/MWh)
            load["cost_curt"] = 3000
            _PM._apply_func!(load, "cost_curt", rescale_cost)
        end
    else
        for (i, load_extra) in data["load_extra"]
            # ID of load point
            idx = load_extra["load_id"]
            delete!.(Ref(load_extra), ["source_id", "load_id", "index"])
            merge!(data["load"]["$idx"], load_extra)

            # Rescale cost and power input values to the p.u. values used internally in the model
            _PM._apply_func!(data["load"]["$idx"], "cost_curt", rescale_cost)
            _PM._apply_func!(data["load"]["$idx"], "cost_red", rescale_cost)
            _PM._apply_func!(data["load"]["$idx"], "cost_shift", rescale_cost)
            delete!(data, "load_extra")
        end
    end
end

function _add_storage_data!(data::Dict{String, Any})
    if haskey(data, "storage")
        for (s, storage) in data["storage"]
            rescale_power = x -> x/data["baseMVA"]
            _PM._apply_func!(storage, "max_energy_absorption", rescale_power)
            _PM._apply_func!(storage, "stationary_energy_outflow", rescale_power)
            _PM._apply_func!(storage, "stationary_energy_inflow", rescale_power)
        end
    else
        data["storage"] = Dict{String, Any}()
    end
end

function _add_pst_data!(data::Dict{String, Any})
    if haskey(data, "pst")
        if !haskey(data, "multinetwork") || data["multinetwork"] == false
            _CBA.to_pu_single_network_pst!(data)
            _CBA.fix_data_single_network_pst!(data)
        else
            _CBA.to_pu_multi_network_pst!(data)
            _CBA.fix_data_multi_network_pst!(data)
        end
    else
        data["pst"] = Dict{String, Any}()
    end
end


## Optimisation helper functions


function get_storage_energy_final(solution::Dict{String, Any})

    nw = maximum(x -> parse(Int, x), keys(solution["nw"]))
    data = Dict{String, Any}()
    for (key, value) in get(solution["nw"][string(nw)], "storage", Dict())
        data[key] = Dict{String, Any}("energy" => value["se"])
    end
    return data
end


function set_storage_energy_initial!(network::Dict{String, Any}, energy::Dict{String, Any})

    nw = minimum(x -> parse(Int, x), keys(network["nw"]))
    if nw > 1
        for (key, value) in network["nw"][string(nw)]["storage"]
            merge!(value, energy[key])
        end
    end
    return nothing
end


function is_feasible(solution::Dict{String, Any})
    return in(string(solution["termination_status"]), ["LOCALLY_SOLVED", "OPTIMAL"])
end
