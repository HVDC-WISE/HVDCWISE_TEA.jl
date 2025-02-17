using Distributed

## Workers setup
nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using HiGHS
    HVDCWISE_TEA.silence()
end

# Directory containing user inputs in the subfolder "user_interface/inputs"
main_dir = @__DIR__  # Directory of the current file
work_dir = joinpath(main_dir, "studies\\study_2")  # Folder 'studies' is ignored by git

# Directory of another macro-scenario. Used to ensure consistency of availability time series between macro-scenarios. If no pre-existing macro-scenario: previous_work_dir = ""
previous_work_dir = joinpath(main_dir, "studies\\study_1")

# Number of contingency time series to generated. N_microscenarios = n_availability_series * N_power_series
n_availability_series = 1

# The yearly problem is split into subproblems of this size (1 week is 168h, 1 month is 720h, 1 year is 8760h)
hours_per_subsimulation = 168

# Base power (in MVA) for the per unit conversion
base_MVA = 100

## Solver parameters
setting_opt = Dict("presolve" => "on", "solver" => "ipm", "run_crossover" => "off", "ipm_iteration_limit" => 3000, "time_limit" => 3600.0, "output_flag" => false)
optimizer = HVDCWISE_TEA.optimizer_with_attributes(HiGHS.Optimizer, setting_opt...)
setting = Dict("output" => Dict("branch_flows" => true, "duals" => false), "conv_losses_mp" => false);

# Verify that Matlab or Octave is installed. Matlab and the "Statistics and Machine Learning Toolbox" are required if you want to consider correlations between pole failures in DC components.
octave_path = HVDCWISE_TEA.detect_octave()  # You can replace this line by the path of your Octave or Matlab launcher. For example: "C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe"

# Run the study
@time HVDCWISE_TEA.run_study(work_dir, previous_work_dir, n_availability_series, hours_per_subsimulation, base_MVA, optimizer, setting, octave_path)
