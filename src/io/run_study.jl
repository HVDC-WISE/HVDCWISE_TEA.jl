


function run_study(work_dir::String, hours::Int, base_mva::Int, optimizer, setting)
    simulation_dir = joinpath(work_dir, "simulation_interface")
    if isdir(simulation_dir)
        rm(simulation_dir, recursive=true)
    end

    println("Build simulation inputs")
    build_simulation_inputs(work_dir, base_mva)

    println("Run simulation")
    run_simulation(work_dir, hours, optimizer, setting)

    println("Run results post-processing")
    build_user_results(work_dir, base_mva)

    println("Study finished")
end


function run_simulation(work_dir::String, hours::Int, optimizer, setting)
    #Paths
    simulation_dir = joinpath(work_dir, "simulation_interface")
    path2grid = find_grid_path(simulation_dir)
    path2data = joinpath(simulation_dir, "Input_series")
    
    ## Solve the multiperiod OPF problem
    run_tea(path2grid, path2data, hours, _PM.DCPPowerModel, optimizer; setting = setting)
end


function find_grid_path(simulation_dir::String)
    matpower_name = ""
    for file_name in readdir(simulation_dir)
        file_path = joinpath(simulation_dir, file_name)
        if isfile(file_path)
            macro_scenario = file_name[1:length(file_name)-2]
            if file_name == "$macro_scenario.m"
                @assert (matpower_name == "")  "2 files could contain the grid model data in $simulation_dir: $macro_scenario.m and $matpower_name. Please delete or rename the wrong file."
                matpower_name = file_name
            end
        end
    end
    @assert (matpower_name != "")  "No matpower file has been found in $simulation_dir."
    return joinpath(simulation_dir, matpower_name)
end

#=
# Code to test the functions in this script:

using Distributed
using HVDCWISE_TEA
import PowerModels as _PM
import Ipopt

const _HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA)));  # Root directory of HVDCWISE_TEA package

## Workers setup

nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using Ipopt
    HVDCWISE_TEA.silence()
end

work_dir = joinpath(_HWTEA_dir, "studies\\simple_use_case")
hours_per_subsimulation = 96  # The yearly problem is split into subproblems of this size (1 week is 168h)
@time run_study(work_dir, hours_per_subsimulation)
=#