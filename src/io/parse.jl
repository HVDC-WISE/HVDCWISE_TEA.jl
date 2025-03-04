
"""
    parse_file(file; <keyword arguments>)

Parse a Matpower .m `file` into an `HVDCWISE_TEA` single-network data structure.

Keyword arguments, if any, are forwarded to `PowerModels.parse_file`.
"""
function parse_file(file::String; kwargs...)

    data = _PM.parse_file(file; kwargs...)
    if haskey(data, "multinetwork") && data["multinetwork"] == true
        Memento.error(_LOGGER, "A single-network model must be provided.")
    end
    process_additional_data!(data)

    return data
end

"""
    parse_data(path2grid::String, path2data::String; kwargs...)

Wrapper function for parsing, sequentially:
 - A Matpower .m file located at `path2grid` into an `HVDCWISE_TEA` single-network data structure;
 - A nested folder structure located at `path2data` into a dictionary of timeseries;
 - The single-network structure into a multi-network one by merging it with the timeseries dictionary, if provided;
 - The single or multi-network structure into a multi-conductor one (for the HVDC grid part, if provided).

"""
function parse_data(path2grid::String, path2data::Vector{String}; kwargs...)

    # Parse initial single-network structure
    sn_data = parse_file(path2grid; kwargs...)
    # Parse time series data
    time_series, hours = parse_timeseries(path2data, sn_data)
    # Add hour dimension to single-network data
    _FP.add_dimension!(sn_data, :hour, hours)
    # Create the multinetwork data dictionary
    mn_data = _FP.make_multinetwork(sn_data, time_series; share_data=false)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    _PMMCDC.make_multiconductor!(mn_data)

    return mn_data
end

function parse_data(path2data::Vector{String}, data_sn::Dict{String, <:Any}, nw_ids::AbstractVector{Int}; kwargs...)

    pop!(data_sn, "dim", nothing)
    data_sn = deepcopy(data_sn)
    # Parse time series data withing given interval
    time_series, _ = parse_timeseries(path2data, data_sn; nw_ids=nw_ids)
    # Add hour dimension to single-network data
    _FP.add_dimension!(data_sn, :hour, length(nw_ids))
    offset_dimension!(data_sn, nw_ids)
    # Create the multinetwork data dictionary
    data_mn = make_multinetwork(data_sn, time_series)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    if !haskey(time_series, "busdc") && haskey(data_sn, "busdc")
        nw = first(nw_ids)
        make_multiconductor!(data_mn["nw"]["$nw"]; components=["busdc"])
        make_multiconductor!(data_mn; components=["convdc", "branchdc"])
    else
        make_multiconductor!(data_mn)
    end

    return data_mn
end

function parse_data(path2grid::String, path2data::Vector{String}, num::Int; kwargs...)

    # Parse initial single-network structure
    sn_data = parse_file(path2grid; kwargs...)
    # Parse time series data
    time_series, hours = parse_timeseries(path2data, sn_data)
    # Keep only first `num` steps of the time series
    num = hours > num ? num : hours
    for comp in values(time_series)
        for var in values(comp)
            for array in values(var)
                keepat!(array, 1:num)
            end
        end
    end
    # Add hour dimension to single-network data
    _FP.add_dimension!(sn_data, :hour, num)
    # Create the multinetwork data dictionary
    mn_data = _FP.make_multinetwork(sn_data, time_series; share_data=false)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    _PMMCDC.make_multiconductor!(mn_data)

    return mn_data, hours
end

function parse_data(path2grid::String; kwargs...)

    # Parse initial single-network structure
    sn_data = parse_file(path2grid; kwargs...)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    _PMMCDC.make_multiconductor!(sn_data)

    return sn_data
end

function parse_timeseries(paths::Vector{String}, network::Dict{String, <:Any}; nw_ids::Union{AbstractVector{Int}, Nothing} = nothing)

    dim = Int[]
    data = Dict{String, Any}()

    for path in paths
        for component in readdir(path)

            if !haskey(data, component)
                # Initialise nested dictionary for given component
                data[component] = Dict(key => Dict{String, Any}() for key in keys(network[component]))
            end

            for file in readdir(joinpath(path, component))

                var = string(first(split(file, ".")))
                # Check if the variable name exists for the given component
                if !in(var, keys(first(values(network[component]))))
                    msg = "Variable name $var does not exists for power system component $component."
                    error(msg)
                end
                # Initialise empty vector of the correct type
                type = typeof(first(values(network[component]))[var])
                for value in values(data[component])
                    value[var] = Vector{type}()
                end

                # Read timeseries data for the given component variable
                if isnothing(nw_ids)
                    df = CSV.read(joinpath(path, component, file), _DF.DataFrame, header=true)
                else
                    df = CSV.read(joinpath(path, component, file), _DF.DataFrame, header=true,
                        skipto=first(nw_ids)+1, limit=length(nw_ids)
                    )
                end
                # Check for NaN in timeseries data
                if any(isnan.(Matrix(df)))
                    msg = "NaN values exist in timeseries for variable $var of component $component."
                    error(msg)
                end
                # Check if timeseries length is consistent across different files
                append!(dim, size(df, 1))
                if length(unique(dim)) > 1
                    msg = "Length of timeseries for variable $var of component $component is inconsistent with other timeseries."
                    error(msg)
                end
                df = Dict(names(df) .=> eachcol(df))

                # Assign value to dictionary
                for (key, value) in data[component]
                    if in(key, keys(df))
                        timeseries = convert(Vector{type}, pop!(df, key))
                    else
                        timeseries = repeat([network[component][key][var]], dim[end])
                    end
                    append!(value[var], timeseries)
                end
                # Check if there are keys that do not exist for the given component
                if !isempty(df)
                    msg = "Timeseries data for variable $var of component $component is defined for elements that do not exist."
                    error(msg)
                end
            end
        end
    end
    return data, dim[begin]
end