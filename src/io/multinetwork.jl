
function offset_dimension!(data::Dict{String, Any}, nw_ids::AbstractVector{Int})

    key = "dim"
    data[key][:li] = collect(nw_ids)
    data[key][:ci] = CartesianIndices(data[key][:li])
    data[key][:offset] = first(nw_ids) - 1
    return nothing
end

function slice_multinetwork(data::Dict{String, Any}, nw_ids::AbstractVector{Int})

    slice = Dict{String,Any}()
    for key in keys(data)
        if key == "nw"
            slice[key] = Dict{String,Any}()
        elseif key == "dim"
            slice[key] = deepcopy(data[key])
        else
            slice[key] = data[key]
        end
    end

    key = "dim"
    slice[key][:li] = collect(nw_ids)
    slice[key][:ci] = CartesianIndices(slice[key][:li])
    slice[key][:offset] = first(nw_ids) - 1

    key = "nw"
    for nw in slice["dim"][:li]
        slice["nw"]["$nw"] = data["nw"]["$nw"]
    end

    return slice
end

"""
    make_multinetwork(sn_data, time_series; <keyword arguments>)

Generate a multinetwork data structure from a single network and a time series.

# Arguments
- `sn_data`: single-network data structure to be replicated.
- `time_series`: data structure containing the time series.
- `global_keys`: keys that are stored once per multinetwork (they are not repeated in each
  `nw`).
- `number_of_nws`: number of networks to be created from `sn_data` and `time_series`;
  default: read from `dim`.
- `nw_id_offset`: optional value to be added to `time_series` ids to shift `nw` ids in
  multinetwork data structure; default: read from `dim`.
- `check_dim`: whether to check for `dim` in `sn_data`; default: `true`.
"""
function make_multinetwork(
        sn_data::Dict{String,Any},
        time_series::Dict{String,Any};
        global_keys = ["dim","multinetwork","name","per_unit","source_type","source_version"],
        number_of_nws::Int = length(_FP.dim(sn_data)[:li]),
        nw_id_offset::Int = _FP.dim(sn_data)[:offset],
        check_dim::Bool = true
    )

    if _IM.ismultinetwork(sn_data)
        Memento.error(_LOGGER, "`sn_data` argument must be a single network.")
    end
    if check_dim && !haskey(sn_data, "dim")
        Memento.error(_LOGGER, "Missing `dim` dict in `sn_data` argument. The function `add_dimension!` must be called before `make_multinetwork`.")
    end

    mn_data = Dict{String,Any}("nw"=>Dict{String,Any}())
    _FP._add_mn_global_values!(mn_data, sn_data, global_keys)
    _add_time_series!(mn_data, sn_data, global_keys, time_series, number_of_nws, nw_id_offset)

    return mn_data
end

# Build multinetwork data structure: for each network, replicate the template and replace with data from time_series
function _add_time_series!(mn_data, sn_data, global_keys, time_series, number_of_nws, offset)
    template_nw = _FP._make_template_nw(sn_data, global_keys)
    for time_series_idx in 1:number_of_nws
        n = time_series_idx + offset
        mn_data["nw"]["$n"] = _build_nw(template_nw, time_series, time_series_idx)
    end
end

# Build the nw by copying the template and substituting data from time_series.
function _build_nw(template_nw, time_series, idx)
    nw = copy(template_nw)
    for (key, element) in time_series
        if haskey(nw, key)
            nw[key] = deepcopy(template_nw[key])
            for (l, element) in time_series[key]
                if haskey(nw[key], l)
                    for (m, property) in time_series[key][l]
                        nw[key][l][m] = property[idx]
                    end
                else
                    Memento.warn(_LOGGER, "Key $l not found, will be ignored.")
                end
            end
        else
            Memento.warn(_LOGGER, "Key $key not found, will be ignored.")
        end
    end
    return nw
end
