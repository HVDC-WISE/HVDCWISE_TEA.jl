
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