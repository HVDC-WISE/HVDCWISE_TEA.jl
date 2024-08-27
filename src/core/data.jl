
"list of all `PowerModels`` unused properties in `HVDCWISE_TEA` optimisation model"
const _UNUSED_VARS = Dict{String, Vector{String}}(
    "gen" => ["apf", "pc1", "pc2", "qc1max", "qc1min", "qc2max", "qc2min", "ramp_10", "ramp_30",  "ramp_agc", "ramp_q", "shutdown", "startup"],
    "load" => ["cost_inv", "ered_rel_max", "eshift_rel_max", "lifetime", "pf_angle", "tshift_down", "tshift_up"],
    "convdc" => ["droop", "P_g", "Q_g", "Vdcset"],
    "branchdc" => ["l", "c"]
)

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

By default, element properties that are not used in `HVDCWISE_TEA` mathematical model are removed.
"""
function process_additional_data!(data::Dict{String, Any})
    
    _FP._add_gen_data!(data)
    _add_dcgrid_data!(data)
    _add_storage_data!(data)
    _add_flexible_demand_data!(data)
    _add_pst_data!(data)
    _delete_unused_properties!(data)
end


## Auxiliary function


function _add_dcgrid_data!(data::Dict{String, Any})
    if haskey(data, "busdc")
        _PMMCDC.process_additional_data!(data)
    end
end

function _add_flexible_demand_data!(data::Dict{String, Any})

    rescale_cost = x -> x*data["baseMVA"]

    if haskey(data, "load_extra")
        cost_curt = maximum(getindex.(values(data["load_extra"]), "cost_curt"))
    else
        cost_curt = 3000
    end

    for (i, load) in data["load"]
        # Whether load is flexible (boolean)
        load["flex"] = 0
        # Compensation for load curtailment (i.e. involuntary demand reduction) (€/MWh)
        load["cost_curt"] = cost_curt
        _PM._apply_func!(load, "cost_curt", rescale_cost)
    end

    if haskey(data, "load_extra")
        for (i, load_extra) in data["load_extra"]
            # ID of load point
            idx = load_extra["load_id"]
            delete!.(Ref(load_extra), ["source_id", "load_id", "index"])
            merge!(data["load"]["$idx"], load_extra)

            # Rescale cost and power input values to the p.u. values used internally in the model
            _PM._apply_func!(data["load"]["$idx"], "cost_curt", rescale_cost)
            _PM._apply_func!(data["load"]["$idx"], "cost_red", rescale_cost)
            _PM._apply_func!(data["load"]["$idx"], "cost_shift", rescale_cost)
        end
        delete!(data, "load_extra")
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

function _delete_unused_properties!(data::Dict{String, Any})

    for (key, names) in _UNUSED_VARS
        delete!.(values(data[key]), permutedims(names))
    end
    return nothing
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
