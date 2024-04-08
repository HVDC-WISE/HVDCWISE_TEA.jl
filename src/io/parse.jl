
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
    parse_data(file::String, time_series::Dict{String, <:Any}; kwargs...)

Wrapper function for parsing, sequentially:
 - A Matpower .m `file` into an `HVDCWISE_TEA` single-network data structure;
 - The single-network structure into a multi-network one by merging it with a `time_series` dictionary, if provided;
 - The single or multi-network structure into a multi-conductor one (for the HVDC grid part).

"""
function parse_data(file::String; kwargs...)

    # Parse initial single-network structure
    sn_data = parse_file(file; kwargs...)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    _PMMCDC.make_multiconductor!(sn_data)

    return sn_data
end

function parse_data(file::String, time_series::Dict{String, <:Any}; kwargs...)

    # Parse initial single-network structure
    sn_data = parse_file(file; kwargs...)
    # Find number of simulated hours (a single dimension is used)
    hours = length(find_value(time_series))
    # Add hour dimension to single-network data
    _FP.add_dimension!(sn_data, :hour, hours)
    # Create the multinetwork data dictionary
    mn_data = _FP.make_multinetwork(sn_data, time_series; share_data=false)
    # Convert the DC grid (if any) of each single-network model to multi-conductor
    _PMMCDC.make_multiconductor!(mn_data)

    return mn_data
end
