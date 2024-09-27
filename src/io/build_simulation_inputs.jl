
using CSV
import DataFrames
using XLSX

function build_simulation_inputs(work_dir::String, base_mva::Int=100)
    model_data = build_grid_model(work_dir, base_mva)
    power_series_info = build_power_series(work_dir, base_mva, model_data)
    n_power_series = power_series_info["n_series"]  # Not used
    n_hours = power_series_info["n_hours"]

    build_availability_series(work_dir, n_hours)
end

function base_to_pu(value::Float64, unit::String, base_mva::Int, base_kv::Int)
    if in(unit, ["MVA", "MW", "MVAr", "MWh"])
        return value / base_mva
    elseif unit == "kV"
        return value / base_kv
    elseif unit == "kA"
        return value *  (sqrt(3) * base_kv) / base_mva
    elseif unit == "Ohm"
        return value * base_mva / base_kv^2
    else
        return value
    end
end

function pu_to_base(value::Float64, unit::String, base_mva::Int, base_kv::Int)
    if unit in ["MVA", "MW", "MVAr", "MWh"]
        return value * base_mva
    elseif unit == "kV"
        return value * base_kv
    elseif unit == "kA"
        return value * base_mva / (sqrt(3) * base_kv)
    elseif unit == "Ohm"
        return value * base_kv^2 / base_mva
    else
        return value
    end
end

function build_grid_model(work_dir::String, base_mva::Int)
    # Identify the macro-scenario name and the paths of the required input files

    user_dir = joinpath(work_dir, "user_interface")
    @assert isdir(user_dir) "The work_dir ($work_dir) should contain a folder 'user_interface'"
    user_inputs_dir = joinpath(user_dir, "inputs")
    @assert isdir(user_inputs_dir) "The user_dir ($user_dir) should contain a folder 'inputs'"

    macro_scenario = ""

    for file_name in readdir(user_inputs_dir)
        if occursin("_model.xlsx", file_name)
            @assert (macro_scenario == "")  "2 files could contain the grid model data: $macro_scenario" * "_model.xlsx and $file_name. Please delete or rename the wrong file."
            macro_scenario = file_name[1:length(file_name)-11]
        end
    end
    @assert (macro_scenario != "") "The file containing the grid model data should end with '_model.xlsx'. This file has not been found in $user_inputs_dir."
    model_path = joinpath(user_inputs_dir, "$macro_scenario"*"_model.xlsx")
    @assert isfile(model_path)  "$model_path is not a file"

    costs_path = joinpath(user_inputs_dir, "costs_data.xlsx")
    @assert isfile(costs_path)  "The file containing the cost data should be 'costs_data.xlsx'. This file has not been found in $user_inputs_dir."
    
    default_path = joinpath(user_inputs_dir, "default_values.xlsx")
    @assert isfile(default_path)  "The file containing the default values should be 'default_values.xlsx'. This file has not been found in $user_inputs_dir."

    # Read default values and build a nested dictionary from it
    
    default_data = Dict("bus" => Dict(), "busdc" => Dict(), "branch" => Dict(), "branch_currents" => Dict(), "branchdc" => Dict(), 
    "convdc" => Dict(), "emission_factors" => Dict(), "gen" => Dict(), "gencost" => Dict(), "ndgen" => Dict(), 
    "load_extra" => Dict(), "storage" => Dict(), "storage_extra" => Dict())
    sorted_default_attributes = Dict(comp_name => Vector{String}() for comp_name in keys(default_data))

    default_file = XLSX.readxlsx(default_path)
    sheet_names = XLSX.sheetnames(default_file)
    @assert "default" in sheet_names  "The default values should be in the sheet 'default' of $default_path. This sheet has not been found."
    default_sheet = default_file["default"]
    n_rows = XLSX.get_dimension(default_sheet).stop.row_number
    n_cols = XLSX.get_dimension(default_sheet).stop.column_number
    @assert n_rows >= 153 && n_cols >= 5  "Sheet 'default' of $default_path should have 153 rows and 5 columns, not $n_rows and $n_cols"
    if n_rows > 153
        @assert (Set([string(default_sheet[i,j]) for i in 153:n_rows for j in 1:5]) == "missing") "Rows > 153 should be empty in sheet 'default' of $default_path. $([default_sheet[i,j] for i in 153:n_rows for j in 1:5])"
    end
    n_rows = 153
    n_cols = 5
    @assert [default_sheet[1,j] for j in 1:5] == ["Component", "Attribute", "Unit", "Description", "Default value"] "$default_path sheet 'default' A1:E1 is $([default_sheet[1,j] for j in 1:5]) instead of [Component, Attribute, Unit, Description, Default value]"
    @assert (Set([default_sheet[i,1] for i in 2:n_rows]) == Set(keys(default_data)))  "components in default values should be $(Set(keys(default_data))), not $(Set([default_sheet[i,1] for i in 2:n_rows]))"
    
    for i in 2:n_rows
        comp_name = default_sheet[i,1]
        attribute = default_sheet[i,2]
        default_data[comp_name][attribute] = Dict("unit" => string(default_sheet[i,3]), "value" => default_sheet[i,5])
        push!(sorted_default_attributes[comp_name], attribute)
    end

    # Read grid model and build a nested dictionary from it

    model_attributes = Dict(
        "bus" => ["base voltage"],
        "busdc" => ["base voltage"], 
        "branch" => ["from bus id",  "to bus id", "type", "power rating", "resistance", "reactance"], 
        "branchdc" => ["from bus id",  "to bus id", "type", "power rating", "resistance"], 
        "convdc" => ["AC bus id",  "DC bus id", "type", "power rating", "configuration", "connection type"],
        "gen" => ["bus id", "type", "power rating"], 
        "load" => ["bus id", "type", "power rating", "max power shift up", "max power shift down", "max voluntary reduced power"], 
        "storage" => ["bus id", "type", "production rating", "consumption rating", "initial energy", "energy rating", "production efficiency", "consumption efficiency", "self discharge rate"])
    model_units = Dict(comp_name => Dict() for comp_name in keys(model_attributes))
    model_data = Dict(comp_name => Dict() for comp_name in keys(model_attributes))

    model_file = XLSX.readxlsx(model_path)
    sheet_names = XLSX.sheetnames(model_file)

    for comp_name in keys(model_data)
        @assert comp_name in sheet_names  "Sheet $comp_name is missing in $model_path"
        comp_sheet = model_file[comp_name]
        n_rows = XLSX.get_dimension(comp_sheet).stop.row_number
        n_cols = XLSX.get_dimension(comp_sheet).stop.column_number
        @assert n_rows >= 3 && n_cols >= length(model_attributes[comp_name])  "Sheet $comp_name of $model_path should have at least 3 rows and $(length(model_attributes[comp_name])) columns, not $n_rows and $n_cols"
        @assert comp_sheet[1,2] == "$comp_name id" "$model_path sheet $comp_name B1 is $(comp_sheet[1,2]) instead of '$comp_name id'"
        @assert [comp_sheet[i,1] for i in 1:3] == ["Attribute", "Unit", "Description"] "$model_path sheet $comp_sheet A1:A3 is $([comp_sheet[i,1] for i in 1:3]) instead of [Attribute, Unit, Description]"
        if n_cols > length(model_attributes[comp_name]) + 2
            @assert (Set([string(comp_sheet[i,j]) for i in 4:n_rows for j in length(model_attributes[comp_name])+3:n_cols]) == Set(["missing"]))  "Columns > $(length(model_attributes[comp_name])+3) in $model_path sheet $comp_sheet should be empty.\n$(Set([string(comp_sheet[i,j]) for i in 4:n_rows for j in length(model_attributes[comp_name])+3:n_cols]))"
            n_cols = length(model_attributes[comp_name]) + 2
        end
        @assert [comp_sheet[1,j] for j in 3:n_cols] == model_attributes[comp_name] "$model_path sheet $comp_name B3:end3 is $([comp_sheet[1,j] for j in 2:4]) instead of $(model_attributes[comp_name])"
        for j in 2:n_cols
            model_units[comp_name][comp_sheet[1,j]] = comp_sheet[2,j]
        end
        for i in 4:n_rows
            if Set([comp_sheet[i,j] for j in 2:n_cols]) != Set(["missing"])  # The row is not empty
                @assert !in("missing", [string(comp_sheet[i,j]) for j in 2:n_cols])  "$model_path sheet $comp_name row $i. Some data is missing: $([comp_sheet[i,j] for j in 2:n_cols])"
                @assert comp_sheet[i,2] == i-3  "$model_path sheet $comp_name row $i. $comp_name id should be $(i-3) instead of $(comp_sheet[i,2])"
                comp_id = i-3
                model_data[comp_name][comp_id] = Dict(comp_sheet[1,j] => comp_sheet[i,j] for j in 3:n_cols)
            end
        end
    end

    # Read grid costs and build a nested dictionary only for load and gen (other data are for KPI computation in post-processing)
    
    costs_attributes = Dict(
        "load" => ["curtailment cost", "reduction cost", "shifting cost"], 
        "gen" => ["production cost", "curtailment cost", "CO2 emission"])
    costs_units = Dict(comp_name => Dict() for comp_name in keys(costs_attributes))
    costs_data = Dict(comp_name => Dict() for comp_name in keys(costs_attributes))

    costs_file = XLSX.readxlsx(costs_path)
    sheet_names = XLSX.sheetnames(costs_file)

    for comp_name in keys(costs_data)
        @assert comp_name in sheet_names  "Sheet $comp_name is missing in $costs_path"
        comp_sheet = costs_file[comp_name]
        n_rows = XLSX.get_dimension(comp_sheet).stop.row_number
        n_cols = XLSX.get_dimension(comp_sheet).stop.column_number
        @assert n_rows >= 3 && n_cols >= length(costs_attributes[comp_name])  "Sheet $comp_name of $costs_path should have at least 3 rows and $(length(costs_attributes[comp_name])) columns, not $n_rows and $n_cols"
        @assert comp_sheet[1,2] == "type" "$costs_path sheet $comp_name B1 is $(comp_sheet[1,2]) instead of 'type'"
        @assert [comp_sheet[i,1] for i in 1:3] == ["Attribute", "Unit", "Description"] "$costs_path sheet $comp_sheet A1:A3 is $([comp_sheet[i,1] for i in 1:3]) instead of [Attribute, Unit, Description]"
        if n_cols > length(costs_attributes[comp_name]) + 2
            @assert Set([string(comp_sheet[i,j]) for i in 4:n_rows for j in length(costs_attributes[comp_name])+3:n_cols]) == Set(["missing"])
            n_cols = length(costs_attributes[comp_name]) + 2
        end
        @assert [comp_sheet[1,j] for j in 3:n_cols] == costs_attributes[comp_name] "$costs_path sheet $comp_name B3:end3 is $([comp_sheet[1,j] for j in 2:4]) instead of $(costs_attributes[comp_name])"
        for j in 2:n_cols
            costs_units[comp_name][comp_sheet[1,j]] = comp_sheet[2,j]
        end
        for i in 4:n_rows
            if Set([comp_sheet[i,j] for j in 2:n_cols]) != Set(["missing"])  # The row is not empty
                @assert !in("missing", [comp_sheet[i,j] for j in 2:n_cols])  "$costs_path sheet $comp_name row $i. Some data is missing: $([comp_sheet[i,j] for j in 2:n_cols])"
                comp_type = comp_sheet[i,2]
                costs_data[comp_name][comp_type] = Dict(comp_sheet[1,j] => comp_sheet[i,j] for j in 3:n_cols)
            end
        end
    end

    # Build a dictionary containing all the required attributes for the matpower (.m) file, in their required units (be careful about when unit conversions are needed)

    matpower_data = Dict("bus" => Dict(), "busdc" => Dict(), "branch" => Dict(), "branch_currents" => Dict(), "branchdc" => Dict(), 
                         "convdc" => Dict(), "emission_factors" => Dict(), "gen" => Dict(), "gencost" => Dict(), "ndgen" => Dict(), 
                         "load_extra" => Dict(), "storage" => Dict(), "storage_extra" => Dict())
    @assert keys(matpower_data) == keys(default_data)

    # bus
    for bus_id in keys(model_data["bus"])
        model_bus_data = model_data["bus"][bus_id]
        matpower_bus_data = Dict()
        for attribute in keys(default_data["bus"])
            matpower_bus_data[attribute] = default_data["bus"][attribute]["value"]
        end
        matpower_bus_data["bus_i"] = bus_id
        @assert default_data["bus"]["baseKV"]["unit"] == model_units["bus"]["base voltage"]
        matpower_bus_data["baseKV"] = model_bus_data["base voltage"]
        matpower_data["bus"][bus_id] = matpower_bus_data
    end

    # busdc
    for busdc_id in keys(model_data["busdc"])
        model_busdc_data = model_data["busdc"][busdc_id]
        matpower_busdc_data = Dict()
        for attribute in keys(default_data["busdc"])
            matpower_busdc_data[attribute] = default_data["busdc"][attribute]["value"]
        end
        matpower_busdc_data["busdc_i"] = busdc_id
        @assert default_data["busdc"]["basekVdc"]["unit"] == model_units["busdc"]["base voltage"]
        matpower_busdc_data["basekVdc"] = model_busdc_data["base voltage"]
        matpower_data["busdc"][busdc_id] = matpower_busdc_data
    end
    
    # branch
    for branch_id in keys(model_data["branch"])
        model_branch_data = model_data["branch"][branch_id]
        matpower_branch_data = Dict()
        for attribute in keys(default_data["branch"])
            matpower_branch_data[attribute] = default_data["branch"][attribute]["value"]
        end
        from_bus = model_branch_data["from bus id"]
        to_bus = model_branch_data["to bus id"]
        from_voltage = model_data["bus"][from_bus]["base voltage"]
        to_voltage = model_data["bus"][to_bus]["base voltage"]
        matpower_branch_data["fbus"] = from_bus
        matpower_branch_data["tbus"] = to_bus

        model_rating = model_branch_data["power rating"]
        model_unit = model_units["branch"]["power rating"]
        default_unit = default_data["branch"]["rateA"]["unit"]
        @assert default_unit == "MVA"
        if model_unit == "kA"
            matpower_branch_data["rateA"] = model_rating * from_voltage * sqrt(3)
        else
            @assert model_unit == "MW"
            matpower_branch_data["rateA"] = model_rating
        end
        matpower_branch_data["rateB"] = matpower_branch_data["rateA"]
        matpower_branch_data["rateC"] = matpower_branch_data["rateA"]

        model_resistance = model_branch_data["resistance"]
        model_unit = model_units["branch"]["resistance"]
        default_unit = default_data["branch"]["r"]["unit"]
        @assert default_unit in ["pu", "p.u."]
        if model_unit == "Ohm"
            @assert from_voltage == to_voltage  "Branch $branch_id has different voltages at both terminals ($from_voltage kV at bus $from_bus and $to_voltage kV at bus $to_bus). It is a transformer: its resistance must be provided in pu (not in Ohm)."
            resistance_ohm = model_resistance
        else
            @assert model_unit == "pu"
            resistance_ohm = pu_to_base(float(model_resistance), "Ohm", matpower_branch_data["rateA"], from_voltage)  # Base power in the model data is the power rating instead of base_mva.
        end
        matpower_branch_data["r"] = base_to_pu(float(resistance_ohm), "Ohm", base_mva, from_voltage)

        model_reactance = model_branch_data["reactance"]
        model_unit = model_units["branch"]["reactance"]
        default_unit = default_data["branch"]["x"]["unit"]
        @assert default_unit in ["pu", "p.u."]
        if model_unit == "Ohm"
            @assert from_voltage == to_voltage  "Branch $branch_id has different voltages at both terminals ($from_voltage kV at bus $from_bus and $to_voltage kV at bus $to_bus). It is a transformer: its reactance must be provided in pu (not in Ohm)."
            reactance_ohm = model_reactance
        else
            reactance_ohm = pu_to_base(float(model_reactance), "Ohm", matpower_branch_data["rateA"], from_voltage)  # Base power in the model data is the power rating instead of base_mva.
        end
        matpower_branch_data["x"] = base_to_pu(float(reactance_ohm), "Ohm", base_mva, from_voltage)

        matpower_data["branch"][branch_id] = matpower_branch_data
    end

    # branch_currents
    for branch_id in keys(model_data["branch"])
        matpower_data["branch_currents"][branch_id] = Dict("c_rating_a" => default_data["branch_currents"]["c_rating_a"]["value"])
    end

    # branch_dc
    for branchdc_id in keys(model_data["branchdc"])
        model_branchdc_data = model_data["branchdc"][branchdc_id]
        matpower_branchdc_data = Dict()
        for attribute in keys(default_data["branchdc"])
            matpower_branchdc_data[attribute] = default_data["branchdc"][attribute]["value"]
        end
        from_busdc = model_branchdc_data["from bus id"]
        to_busdc = model_branchdc_data["to bus id"]
        from_voltage = model_data["busdc"][from_busdc]["base voltage"]
        to_voltage = model_data["busdc"][to_busdc]["base voltage"]
        @assert from_voltage == to_voltage  "Branchdc $branch_id has different voltages at both terminals ($from_voltage kV at busdc $from_busdc and $to_voltage kV at busdc $to_busdc)"
        matpower_branchdc_data["fbusdc"] = from_busdc
        matpower_branchdc_data["tbusdc"] = to_busdc

        model_rating = model_branchdc_data["power rating"]
        model_unit = model_units["branchdc"]["power rating"]
        default_unit = default_data["branchdc"]["rateA"]["unit"]
        @assert default_unit == "MVA"
        if model_unit == "kA"
            matpower_branchdc_data["rateA"] = model_rating * from_voltage * 2  # TODO check pu conversion methodology in DC
        else
            @assert model_unit == "MW"
            matpower_branchdc_data["rateA"] = model_rating
        end
        matpower_branchdc_data["rateB"] = matpower_branchdc_data["rateA"]
        matpower_branchdc_data["rateC"] = matpower_branchdc_data["rateA"]

        model_resistance = model_branchdc_data["resistance"]
        model_unit = model_units["branchdc"]["resistance"]
        default_unit = default_data["branchdc"]["r"]["unit"]
        @assert default_unit in ["pu", "p.u."]
        if model_unit == "Ohm"
            # P = 2*V*I. V=R*I. P=2*V^2/R. R_base=2*V_base^2/P_base. TODO check pu conversion methodology.
            resistance_ohm = model_resistance
        else
            @assert model_unit == "pu"
            resistance_ohm = model_resistance * 2 * from_voltage^2 / matpower_branchdc_data["rateA"]  # Base power in the model data is the power rating instead of base_mva.
        end
        matpower_branchdc_data["r"] = resistance_ohm / (2 * from_voltage^2 / base_mva)            

        matpower_data["branchdc"][branchdc_id] = matpower_branchdc_data
    end

    # convdc
    for conv_id in keys(model_data["convdc"])
        model_conv_data = model_data["convdc"][conv_id]
        matpower_conv_data = Dict()
        for attribute in keys(default_data["convdc"])
            matpower_conv_data[attribute] = default_data["convdc"][attribute]["value"]
        end
        matpower_conv_data["busdc_i"] = model_conv_data["DC bus id"]
        matpower_conv_data["busac_i"] = model_conv_data["AC bus id"]
        matpower_conv_data["Pacmax"] = base_to_pu(float(model_conv_data["power rating"]), "MW", base_mva, 0)
        matpower_conv_data["Pacmin"] = - matpower_conv_data["Pacmax"]
        matpower_conv_data["conv_confi"] = model_conv_data["configuration"]
        matpower_conv_data["connect_at"] = model_conv_data["connection type"]

        matpower_data["convdc"][conv_id] = matpower_conv_data
    end

    # gen, gencost, ndgen, emission factors
    for gen_id in keys(model_data["gen"])
        model_gen_data = model_data["gen"][gen_id]
        gen_type = model_gen_data["type"]
        costs_gen_data = costs_data["gen"][gen_type]
        cost_production = costs_gen_data["production cost"]
        cost_curtailment = costs_gen_data["curtailment cost"]
        if cost_curtailment == 0  # Dispatchable generator
            matpower_gen_data = Dict()
            for attribute in keys(default_data["gen"])
                matpower_gen_data[attribute] = default_data["gen"][attribute]["value"]
            end
            matpower_gen_data["bus"] = model_gen_data["bus id"]
            matpower_gen_data["Pmax"] = model_gen_data["power rating"]
            matpower_data["gen"][gen_id] = matpower_gen_data

            matpower_gencost_data = Dict()
            for attribute in keys(default_data["gencost"])
                matpower_gencost_data[attribute] = default_data["gencost"][attribute]["value"]
            end
            matpower_gencost_data["c(1)"] = cost_production
            matpower_data["gencost"][gen_id] = matpower_gencost_data
        else  # Non dispatchable generator
            matpower_ndgen_data = Dict()
            for attribute in keys(default_data["ndgen"])
                matpower_ndgen_data[attribute] = default_data["ndgen"][attribute]["value"]
            end
            matpower_ndgen_data["gen_bus"] = model_gen_data["bus id"]
            matpower_ndgen_data["pref"] = model_gen_data["power rating"]
            matpower_ndgen_data["cost_gen"] = cost_production
            matpower_ndgen_data["cost_curt"] = cost_curtailment
            matpower_data["ndgen"][gen_id] = matpower_ndgen_data
        end
        matpower_data["emission_factors"][gen_id] = Dict("CO2" => costs_gen_data["CO2 emission"])  # Other types of emissions (NOx, SO2, COV, particles, ...) could be added
    end

    # load (bus), load_extra
    for load_id in keys(model_data["load"])
        model_load_data = model_data["load"][load_id]
        bus_id = model_load_data["bus id"]
        matpower_data["bus"][bus_id]["Pd"] = model_load_data["power rating"]
        load_type = model_load_data["type"]
        costs_load_data = costs_data["load"][load_type]
        if Set([model_load_data["max power shift up"], model_load_data["max power shift down"], model_load_data["max voluntary reduced power"]]) != Set([0])  # Flexible load
            matpower_load_extra_data = Dict()
            for attribute in keys(default_data["load_extra"])
                matpower_load_extra_data[attribute] = default_data["load_extra"][attribute]["value"]
            end
            matpower_load_extra_data["load_id"] = load_id
            matpower_load_extra_data["pshift_up_rel_max"] = model_load_data["max power shift up"]
            matpower_load_extra_data["pshift_down_rel_max"] = model_load_data["max power shift down"]
            matpower_load_extra_data["pred_rel_max"] = model_load_data["max voluntary reduced power"]
            matpower_load_extra_data["cost_shift"] = costs_load_data["shifting cost"]
            matpower_load_extra_data["cost_red"] = costs_load_data["reduction cost"]
            matpower_load_extra_data["cost_curt"] = costs_load_data["curtailment cost"]
            matpower_data["load_extra"][load_id] = matpower_load_extra_data
        end
    end

    # storage, storage_extra
    for storage_id in keys(model_data["storage"])
        model_storage_data = model_data["storage"][storage_id]

        matpower_storage_data = Dict()
        for attribute in keys(default_data["storage"])
            matpower_storage_data[attribute] = default_data["storage"][attribute]["value"]
        end
        matpower_storage_data["storage_bus"] = model_storage_data["bus id"]
        matpower_storage_data["energy"] = model_storage_data["initial energy"]
        matpower_storage_data["energy_rating"] = model_storage_data["energy rating"]
        matpower_storage_data["charge_rating"] = model_storage_data["consumption rating"]
        matpower_storage_data["discharge_rating"] = model_storage_data["production rating"]
        matpower_storage_data["charge_efficiency"] = model_storage_data["consumption efficiency"]
        matpower_storage_data["discharge_efficiency"] = model_storage_data["production efficiency"]
        matpower_storage_data["thermal_rating"] = max(matpower_storage_data["charge_rating"], matpower_storage_data["discharge_rating"])
        matpower_data["storage"][storage_id] = matpower_storage_data

        matpower_storage_extra_data = Dict()
        for attribute in keys(default_data["storage_extra"])
            matpower_storage_extra_data[attribute] = default_data["storage_extra"][attribute]["value"]
        end
        matpower_storage_extra_data["self_discharge_rate"] = model_storage_data["self discharge rate"]
        matpower_data["storage_extra"][storage_id] = matpower_storage_extra_data
    end

    # Check no missing value
    for matrix_name in keys(matpower_data)
        for id in keys(matpower_data[matrix_name])
            for attribute in keys(matpower_data[matrix_name][id])
                @assert (string(matpower_data[matrix_name][id][attribute]) != "missing") "matpower_data[$matrix_name][$id][$attribute] is missing"
            end
        end
    end
        
    # Builds the matpower (.m) file

    simulation_dir = joinpath(work_dir, "simulation_interface")
    mkpath(simulation_dir)  # Create the folder (and its parents) if it does not exist yet
    matpower_path = joinpath(simulation_dir, "$macro_scenario.m")

    open(matpower_path, "w") do f
        write(f, "function mpc = $macro_scenario()\n\n")
        write(f, "%% MATPOWER Case Format : Version 2\nmpc.version = '2';\n\n")
        write(f, "%% system MVA base\nmpc.baseMVA = $base_mva;\n\n")
        write(f, "mpc.time_elapsed = 1.0;\n")
        for matrix_name in keys(matpower_data)
            if length(matpower_data[matrix_name]) > 0
                sorted_ids = sort([id for id in keys(matpower_data[matrix_name])])
                attributes = sorted_default_attributes[matrix_name]
                write(f, "\n%column_names% $(join(attributes, " "))\n")
                write(f, "mpc.$matrix_name = [\n")
                matrix_data = matpower_data[matrix_name]
                for id in sorted_ids
                    row_values = [matrix_data[id][attribute] for attribute in attributes]
                    write(f, join(row_values, " ") * ";\n")
                end
                write(f, "];\n")
            end
        end
        write(f, "end\n")
    end
    matlab_tool_path = joinpath(pwd(), "src", "matlab_tools")
    if isdir(joinpath(matlab_tool_path, "availability_series"))
        rm(joinpath(matlab_tool_path, "availability_series"), recursive=true)
    end
    cp(matpower_path, joinpath(matlab_tool_path, "grid_model.m"), force=true)
    return model_data
end

function build_csv(micro_scenario_dir::String, sheet::XLSX.Worksheet, readme_data::Dict, model_data, base_mva)
    @assert occursin("|", sheet.name) "Not expected sheet name $(sheet.name)"
    (comp_name, attribute) = split(sheet.name, "|")
    mkpath(joinpath(micro_scenario_dir, comp_name))  # Create the folder (and its parents) if it does not exist yet
    micro_scenario = basename(micro_scenario_dir)
    n_comp = readme_data["Nb of components"]
    n_hours = readme_data["Nb of hours"]
    @assert n_comp > 0  "There should be at least 1 component for $comp_name in micro-scenario $micro_scenario"
    @assert n_hours > 0  "There should be at least 1 hour in micro-scenario $micro_scenario"
    unit = readme_data["Unit"]

    # Check columns & rows names
    model_ids = Set(keys(model_data[comp_name]))
    series_ids = sheet[1, 2:n_comp+1]
    @assert n_comp == length(series_ids)  "Number of $comp_name ids in micro-scenario $micro_scenario is inconsistent in the associated shhet ($(length(series_ids))) and the ReadMe ($n_comp)"
    @assert n_comp == length(model_ids)  "Number of $comp_name ids in micro-scenario $micro_scenario ($n_comp) is not coherent with the associated model ($(length(model_ids)))"
    @assert Set(series_ids) == Set(1:n_comp)  "$comp_name ids in micro-scenario $micro_scenario should be 1:$n_comp instead of $series_ids"
    @assert model_ids == Set(1:n_comp)  "$comp_name ids in micro-scenario $micro_scenario (1:$n_comp) are not coherent with the associated model ($model_ids).\nMissing ids in the model: $(setdiff(Set(1:n_comp), model_ids)).\nExcess ids in the model: $(setdiff(model_ids, Set(1:n_comp)))."
    @assert (sheet[1,1] == "Time\\Id")  "Cell A1 of sheet $comp_name|$attribute of micro-scenario $micro_scenario should be 'Time\\Id', not $(sheet[1,1])"
    @assert [sheet[i,1] for i in 2:n_hours+1] == 1:n_hours  "Hours ids in sheet $comp_name|$attribute of micro-scenario $micro_scenario should be 1:$n_hours, not $(sheet[2:n_hours+1,1])"

    if attribute in ["inflow", "outflow"]
        attribute = "stationary_energy_$attribute"
    end

    # Create CSV
    csv_path = joinpath(micro_scenario_dir, comp_name, "$attribute.csv")
    # series_data = DataFrames.DataFrame(Dict("$id"=>[] for id in 1:n_comp))
    series_data = DataFrames.DataFrame()

    for id in 1:n_comp
        if unit == "MW"
            unit_conversion = 1/base_mva
        else
            @assert unit == "%"
            if comp_name in ["load", "gen"]
                power_rating = model_data[comp_name][id]["power rating"]
            else
                @assert comp_name == "storage"
                if attribute == "inflow"
                    power_rating = model_data["storage"][id]["production rating"] / model_data["storage"][id]["production efficiency"]
                else
                    @assert attribute == "outflow"
                    power_rating = model_data["storage"][id]["consumption rating"] * model_data["storage"][id]["consumption efficiency"]
                end
            end
            unit_conversion = power_rating/base_mva
        end
        series_data[!,"$id"] = [sheet[i,id+1]*unit_conversion for i in 2:n_hours+1]
    end
    CSV.write(csv_path, series_data, writeheader=true)
    return Dict("micro-scenario" => micro_scenario, "comp_name" => comp_name, "attribute" => attribute, "n_comp" => n_comp, "n_hours" => n_hours)
end

function build_power_series(work_dir::String, base_mva::Int, model_data::Dict)
    user_dir = joinpath(work_dir, "user_interface")
    @assert isdir(user_dir) "The work_dir ($work_dir) should contain a folder 'user_interface'"
    user_inputs_dir = joinpath(user_dir, "inputs")
    @assert isdir(user_inputs_dir) "The user_dir ($user_dir) should contain a folder 'inputs'"

    micro_scenarios = Vector{String}()

    for file_name in readdir(user_inputs_dir)
        if occursin("_series.xlsx", file_name)
            push!(micro_scenarios, file_name[1:length(file_name)-12])
        end
    end
    @assert length(micro_scenarios) > 0  "The files containing the power time series should end with '_series.xlsx'. No such file has been found in $user_inputs_dir."
    n_hours = 0
    error_message_info = ""
    for micro_scenario in micro_scenarios
        micro_scenario_dir = joinpath(work_dir, "simulation_interface", "Input_series", "Power", micro_scenario)
        series_path = joinpath(user_inputs_dir, "$micro_scenario"*"_series.xlsx")
        @assert isfile(series_path)  "$series_path is not a file"
        series_file = XLSX.readxlsx(series_path)
        sheet_names = XLSX.sheetnames(series_file)

        @assert "ReadMe" in sheet_names  "file $series_path should have a sheet 'ReadMe'"
        readme_sheet = series_file["ReadMe"]
        @assert [readme_sheet[7,j] for j in 1:4] == ["Attribute", "Nb of components", "Nb of hours", "Unit"]  "ReadMe sheet A7:D7 should be [Attribute, Nb of components, Nb of hours, Unit] while it is $(readme_sheet[7,1:4])"
        @assert [readme_sheet[i,1] for i in 8:11] == ["load|pd", "gen|pmax", "storage|inflow", "storage|outflow"]  "ReadMe sheet A8:A11 should be [load|pd, gen|pmax, storage|inflow, storage|outflow] while it is $(readme_sheet[8:11,1])"
        
        # Build 1 csv per component attribute with time series
        @assert all([readme_sheet[i,2] > 0 for i in 8:9])  "Some time series should be provided both for load|pd and gen|Pmax in micro-scenario $micro_scenario. Check sheet ReadMe in $series_path."
        for i in 8:11
            if readme_sheet[i,2] > 0
                sheet_name = readme_sheet[i,1]
                sheet = series_file[sheet_name]
                csv_info = build_csv(micro_scenario_dir, sheet, Dict(readme_sheet[7,j] => readme_sheet[i,j] for j in 2:4), model_data, base_mva)
                csv_hours = csv_info["n_hours"]
                if n_hours == 0
                    n_hours = csv_hours
                    error_message_info = csv_info
                else
                    @assert n_hours == csv_hours  "Incoherent number of hours in csv files\n$csv_info\n$error_message_info"
                end
            end
        end
    end
    return Dict("n_series" => length(micro_scenarios), "n_hours" => n_hours)
end

function build_availability_series(work_dir::String, n_hours::Int)
    reliability_path = joinpath(work_dir, "user_interface", "inputs", "reliability_data.xlsx")
    @assert isfile(reliability_path)  "$reliability_inputs_path is not a file"
    reliability_file = XLSX.readxlsx(reliability_path)
    sheet_names = XLSX.sheetnames(reliability_file)
    @assert "user_inputs" in sheet_names  "file $reliability_path should have a sheet 'user_inputs'"
    reliability_sheet = reliability_file["user_inputs"]

    reliability_data = DataFrames.DataFrame(reliability_sheet[:],:auto)
    matlab_tool_path = joinpath(pwd(), "src", "matlab_tools")
    reliability_csv_path = joinpath(matlab_tool_path, "reliability_data.csv")
    CSV.write(reliability_csv_path, reliability_data, writeheader=false)

    # series_data[!,"$id"] = [sheet[i,id+1]*unit_conversion for i in 2:n_hours+1]
    

    # TODO build reliability_data.csv in matlab_tool_path from reliability_data.xlsx in work

    println("Run build_availability_series.m in $matlab_tool_path. Then write 'y' and press twice ENTER.")
    a = readline();  # TODO automatically run src/matlab_tools/build_availability_series.m or recode it in Julia
    println("Building availability series")

    availability_series_dir = joinpath(matlab_tool_path, "availability_series")
    for microscenario in  readdir(availability_series_dir)
        micro_dir = joinpath(availability_series_dir, microscenario)
        for comp_name in readdir(micro_dir)
            comp_dir = joinpath(micro_dir, comp_name)
            for csv_name in readdir(comp_dir)
                csv_path = joinpath(comp_dir, csv_name)
                @assert occursin(".csv", csv_name)  "$matlab_tool_path/availability_series/$microscenario/$comp_name/$csv_name is not a CSV file"
                csv_data = CSV.File(csv_path, delim=',') |> DataFrames.DataFrame
                csv_hours = size(csv_data)[1]
                @assert csv_hours >= n_hours  "$csv_hours hours are available in $csv_path. You cannot have $n_hours hours."
                used_data = csv_data[1:n_hours,:]
                new_comp_dir = joinpath(work_dir, "simulation_interface", "Input_series", "Availability", microscenario, comp_name)
                mkpath(new_comp_dir)
                new_csv_path = joinpath(new_comp_dir, csv_name)
                CSV.write(new_csv_path, used_data, writeheader=true)
            end
        end
    end
    # rm(joinpath(matlab_tool_path, "grid_model.m"))  # This file will be deleted after results post-processing
    rm(joinpath(matlab_tool_path, "reliability_data.csv"))
    rm(joinpath(matlab_tool_path, "availability_series"), recursive=true)
end

#=
# Code to test the functions in this script:
println("Builds simulation inputs")
build_simulation_inputs(joinpath(pwd(), "studies\\simple_use_case"))
println("Simulation inputs building finished")
=#