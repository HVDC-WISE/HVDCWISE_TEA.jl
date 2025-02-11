## Flexible loads

"Add to `ref` the keys for handling flexible demand"
function ref_add_flex_load!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})

    for (n, nw_ref) in ref[:it][_PM.pm_it_sym][:nw]
        # Loads that can be made flexible, depending on investment decision
        nw_ref[:flex_load] = Dict(x for x in nw_ref[:load] if x.second["flex"] == 1)
        # Loads that are not flexible and do not have an associated investment decision
        nw_ref[:fixed_load] = Dict(x for x in nw_ref[:load] if x.second["flex"] == 0)
    end

    # Compute the total energy demand of each flex load and store it in the first hour nw
    if haskey(data, "dim")
        for nw in _FP.nw_ids(data; hour = 1)
            if haskey(ref[:it][_PM.pm_it_sym][:nw][nw], :time_elapsed)
                time_elapsed = ref[:it][_PM.pm_it_sym][:nw][nw][:time_elapsed]
            else
                Memento.warn(_LOGGER, "network data should specify time_elapsed, using 1.0 as a default")
                time_elapsed = 1.0
            end
            timeseries_nw_ids = _FP.similar_ids(data, nw, hour = 1:_FP.dim_length(data,:hour))
            for (l, load) in ref[:it][_PM.pm_it_sym][:nw][nw][:flex_load]
                # `ref` instead of `data` must be used to access loads, since the former has
                # already been filtered to remove inactive loads.
                load["ed"] = time_elapsed * sum(ref[:it][_PM.pm_it_sym][:nw][n][:load][l]["pd"] for n in timeseries_nw_ids)
            end
        end
    end
end