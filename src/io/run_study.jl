# Manual steps to run a simulation
# - Julia run build_grid_model() & build_power_series()
# - copy matpower file into contingency_generation/
# - Matlab run contingency_generation/FinalScriptHvdcWise.m
# - Julia run build_availability_series()
# - Julia run run_simulation()

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


# import HVDCWISE_TEA as _HWTEA

function run_study_old(work_dir::String, case_name::String)
    base_mva = 100
    # build_raw_inputs(work_dir) -> the .m at the root and 1 folder per micro with the .csv
    build_raw_inputs(work_dir, work_dir, case_name, base_mva, true)
    #
    # run_tea

    path2grid = joinpath(work_dir, "$case_name.m")
    path2data = joinpath(work_dir, case_name)

    optimizer = _HWTEA.optimizer_with_attributes(Ipopt.Optimizer)
    setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

    run_tea(path2grid, path2data, _PM.DCPPowerModel, optimizer; setting = setting)  # FIXME csv in pu or in MW ? # FIXME uncomment

    # build_outputs_from_csv
    results_dir = joinpath(work_dir, case_name, "results")
    build_outputs_from_csv(results_dir, base_mva, true) # FIXME uncomment
    #
end

# @time run_study(joinpath(_HWTEA_dir, "studies\\IEEE39\\3_hours_without_storage"), "IEEE39")

function run_study(work_dir::String, hours::Int, n_availability_series::Int, base_mva::Int=100)
    println("Build simulation inputs")
    build_simulation_inputs(work_dir, n_availability_series, base_mva)

    println("Run simulation")
    run_simulation(work_dir, hours, base_mva)

    println("Simulation finished")
    # TODO add post-processing here
end

function run_simulation(work_dir::String, hours::Int, base_mva::Int=100)
    #Paths
    simulation_dir = joinpath(work_dir, "simulation_interface")
    path2grid = find_grid_path(simulation_dir)
    path2data = joinpath(simulation_dir, "Input_series")

    ## Solver parameters
    optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);
    
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

#
# Code to test the functions in this script:
work_dir = joinpath(_HWTEA_dir, "studies\\simple_use_case")
hours_per_subsimulation = 96  # The yearly problem is split into subproblems of this size (1 week is 168h)
n_availability_series = 3
@time run_study(work_dir, hours_per_subsimulation, n_availability_series)
#
