using PlotlyJS

function create_plotly_figure(work_dir::String, model_data::Dict)
    # model_data : Dict(sheet_name => Dict(id => Dict(col_name => value)))

    # Check that all buses have coordinates set. If not, map generation is skpipped
    has_coord = true
    n_ac_bus = length(model_data["bus"])
    n_dc_bus = length(model_data["busdc"])
    for bus_id in 1:n_ac_bus
        if ismissing(model_data["bus"][bus_id]["latitude"]) || ismissing(model_data["bus"][bus_id]["longitude"])
            has_coord = false
            break
        end
    end
    for busdc_id in 1:n_dc_bus
        if ismissing(model_data["busdc"][busdc_id]["latitude"]) || ismissing(model_data["busdc"][busdc_id]["longitude"])
            has_coord = false
            break
        end
    end

    if !has_coord
        println("Coordinates are missing for some nodes. No map will be generated")
        return
    end

    # Create bus scatter plot data
    bus_latitudes = []
    bus_longitudes = []
    bus_infos = []
    for bus_id in 1:n_ac_bus
        push!(bus_latitudes, model_data["bus"][bus_id]["latitude"])
        push!(bus_longitudes, model_data["bus"][bus_id]["longitude"])
        bus_text = "AC bus<br>Bus id: $bus_id<br>Rated voltage: $(model_data["bus"][bus_id]["base voltage"])"
        push!(bus_infos, bus_text)
    end
    for busdc_id in 1:n_dc_bus
        push!(bus_latitudes, model_data["busdc"][busdc_id]["latitude"])
        push!(bus_longitudes, model_data["busdc"][busdc_id]["longitude"])
        bus_text = "DC bus<br>Bus id: $busdc_id<br>Rated voltage: $(model_data["busdc"][busdc_id]["base voltage"]) kV"
        push!(bus_infos, bus_text)
    end
    bus_scatter_plot = scattergeo(lat=bus_latitudes, lon=bus_longitudes, mode="markers", hovertext=bus_infos, name="Buses")

    # Create branch line plot data
    n_ac_branch = length(model_data["branch"])
    n_dc_branch = length(model_data["branchdc"])
    branches_latitudes = []
    branches_longitudes = []
    branches_infos_latitudes = []
    branches_infos_longitude = []
    branches_infos = []
    for branch_id in 1:n_ac_branch
        from_bus_id = model_data["branch"][branch_id]["from bus id"]
        to_bus_id = model_data["branch"][branch_id]["to bus id"]

        from_bus_latitude = bus_latitudes[from_bus_id]
        from_bus_longitude = bus_longitudes[from_bus_id]
        to_bus_latitude = bus_latitudes[to_bus_id]
        to_bus_longitude = bus_longitudes[to_bus_id]
        append!(branches_latitudes, [from_bus_latitude, to_bus_latitude, nothing])
        append!(branches_longitudes, [from_bus_longitude, to_bus_longitude, nothing])

        edge_text = "AC branch<br>Branch id: $branch_id<br>" *
                    "Power rating: $(model_data["branch"][branch_id]["power rating"]) MW<br>" *
                    "Resistance: $(model_data["branch"][branch_id]["resistance"]) Ω<br>" *
                    "Reactance: $(model_data["branch"][branch_id]["reactance"]) Ω<br>" *
                    "Lenght: $(model_data["branch"][branch_id]["length"]) km<br>" *
                    "Parallel lines: $(model_data["branch"][branch_id]["number of lines"])"
        push!(branches_infos_latitudes, (from_bus_latitude + to_bus_latitude) / 2)
        push!(branches_infos_longitude, (from_bus_longitude + to_bus_longitude) / 2)
        push!(branches_infos, edge_text)
    end

    for branchdc_id in 1:n_dc_branch
        from_busdc_id = model_data["branchdc"][branchdc_id]["from bus id"]
        to_busdc_id = model_data["branchdc"][branchdc_id]["to bus id"]

        from_busdc_latitude = bus_latitudes[n_ac_bus+from_busdc_id]
        from_busdc_longitude = bus_longitudes[n_ac_bus+from_busdc_id]
        to_busdc_latitude = bus_latitudes[n_ac_bus+to_busdc_id]
        to_busdc_longitude = bus_longitudes[n_ac_bus+to_busdc_id]
        append!(branches_latitudes, [from_busdc_latitude, to_busdc_latitude, nothing])
        append!(branches_longitudes, [from_busdc_longitude, to_busdc_longitude, nothing])

        edge_text = "DC branch<br>Branch id: $branchdc_id<br>" *
                    "Power rating: $(model_data["branchdc"][branchdc_id]["power rating"]) MW<br>" *
                    "Resistance: $(model_data["branchdc"][branchdc_id]["resistance"]) Ω<br>" *
                    "Configuration: $(model_data["branchdc"][branchdc_id]["configuration"]==1 ? "Monopolar" : "Bipolar")<br>" *
                    "Lenght: $(model_data["branchdc"][branchdc_id]["length"]) km<br>" *
                    "Parallel lines: 1"
        push!(branches_infos_latitudes, (from_busdc_latitude + to_busdc_latitude) / 2)
        push!(branches_infos_longitude, (from_busdc_longitude + to_busdc_longitude) / 2)
        push!(branches_infos, edge_text)
    end
    branches_lines_plot = scattergeo(lat=branches_latitudes, lon=branches_longitudes, mode="lines", hoverinfo="skip", name="Branches")
    branches_infos_scatter_plot = scattergeo(lat=branches_infos_latitudes, lon=branches_infos_longitude, mode="text", hovertext=branches_infos, showlegend=false, hoverinfo="text")

    # Create plot layout
    layout = Layout(
        Dict(
            "title" => "Grid Map",
            "xaxis" => Dict("title" => "Longitude (°)"),
            "yaxis" => Dict("title" => "Latitude (°)"),
            "hovermode" => "closest",
            "geo" => Dict(
                "showcountries" => true,  # Show borders
                "fitbounds" => "locations",  # Center/scale view on node location
                "projection" => Dict("type" => "conic equidistant"),
            )
        )
    )

    # Create and save plot
    model_plot = plot(
        [bus_scatter_plot, branches_lines_plot, branches_infos_scatter_plot],
        layout
    )
    output_file_path = joinpath(work_dir, "user_interface", "inputs", "model.html")
    open(output_file_path, "w") do io
        PlotlyBase.to_html(io, model_plot.plot)
    end
    println("Model map saved in $work_dir/user_interface/inputs/model.html")
end