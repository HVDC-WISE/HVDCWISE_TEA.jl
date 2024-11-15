t0 = time()
using Distributed
using HVDCWISE_TEA
import PowerModels as _PM
import HiGHS
import Ipopt
using XLSX

## Workers setup
nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using HiGHS
    HVDCWISE_TEA.silence()
end

# --- Util functions ---

function solver_data(solver)
    if solver == "highs"
        setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);
        setting_opt = Dict("presolve" => "on", "solver" => "ipm", "run_crossover" => "off", "ipm_iteration_limit" => 3000, "time_limit" => 3600.0, "output_flag" => false)
        optimizer = HVDCWISE_TEA.optimizer_with_attributes(HiGHS.Optimizer, setting_opt...)
    else
        @assert solver == "ipopt"
        optimizer = HVDCWISE_TEA.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);
    end
    return Dict("optimizer" => optimizer, "setting" => setting)
end

function compute_time(work_dir::String, hours_per_subsimulation::Int, base_MVA::Int, optimizer, setting, matlab_octave_path::String)
    t1 = time()
    run_study(work_dir, subperiod, 100, optimizer, setting, matlab_octave_path)
    seconds = trunc(Int, time()-t1)
    minutes = trunc(Int, seconds / 60)
    println("Simulation $(basename(work_dir)) finished in $seconds s ($minutes mn)")

    inputs_dir = joinpath(work_dir, "user_interface", "inputs")
    
    micro_scenarios = Vector{String}()
    for file_name in readdir(inputs_dir)
        if occursin("_series.xlsx", file_name)
            push!(micro_scenarios, file_name)
        end
    end
    n_power_series = length(micro_scenarios)
    power_series_file = XLSX.readxlsx(joinpath(inputs_dir, micro_scenarios[1]))
    series_length = power_series_file["ReadMe"][8,3]

    reliability_file = XLSX.readxlsx(joinpath(inputs_dir, "reliability_data.xlsx"))
    reliability_sheet = reliability_file["user_inputs"]
    n_availability_series = reliability_sheet[2,2]

    XLSX.openxlsx(joinpath(main_dir, "Computation_time.xlsx"), mode="rw") do xf
        sheet = xf["Results"]
        row = XLSX.get_dimension(sheet).stop.row_number + 1
        sheet[row, 1] = work_dir  # Folder
        sheet[row, 2] = series_length  # Time series length (n hours)
        sheet[row, 3] = n_power_series  # Number of micro-scenarios for power series
        sheet[row, 4] = n_availability_series  # Number of micro-scenarios for availability series
        sheet[row, 5] = subperiod  # Time series length of each multi-period OPF
        sheet[row, 6] = seconds
    end
end

# --- Main script ---

# Directory containing user inputs in the subfolder "user_interface/inputs"
HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA)))  # Root path of the project
work_dir = joinpath(HWTEA_dir, "studies\\simple_use_case")  # Folder 'studies' is ignored by git

# The yearly problem is split into subproblems of this size (1 week is 168h, 1 month is 720h, 1 year is 8760h)
hours_per_subsimulation = 168

# Base power (in MVA) for the per unit conversion
base_MVA = 100

## Solver parameters
solver = solver_data("highs")
optimizer = solver["optimizer"]
setting = solver["setting"]

# Verify that Matlab or Octave is installed
matlab_octave_path = detect_matlab_or_octave()  # You can replace this line by the path of your Matlab/Octave launcher. for example: "C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe"

# Run the study
@time run_study(work_dir, hours_per_subsimulation, base_MVA, optimizer, setting, matlab_octave_path)

#=
println("Julia initialization time: $(trunc(Int, time()-t0)) s")
for fifo in readdir(main_dir)
    work_dir = joinpath(main_dir, fifo)
    if isdir(work_dir)
        compute_time(work_dir, hours_per_subsimulation, base_MVA, optimizer, setting, matlab_octave_path)
    end
end
=#