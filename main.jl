using Distributed
using HVDCWISE_TEA
import HiGHS

## Workers setup
nprocs() > 1 && rmprocs(workers())
addprocs(Sys.CPU_THREADS ÷ 2; exeflags = "--project=$(Base.active_project())")
@everywhere begin
    using HVDCWISE_TEA
    using HiGHS
    HVDCWISE_TEA.silence()
end

# Directory containing user inputs in the subfolder "user_interface/inputs"
HWTEA_dir = dirname(dirname(pathof(HVDCWISE_TEA)))  # Root path of the project
work_dir = joinpath(HWTEA_dir, "studies\\simple_use_case")  # Folder 'studies' is ignored by git

# The yearly problem is split into subproblems of this size (1 week is 168h, 1 month is 720h, 1 year is 8760h)
hours_per_subsimulation = 168

# Base power (in MVA) for the per unit conversion
base_MVA = 100

## Solver parameters
setting_opt = Dict("presolve" => "on", "solver" => "ipm", "run_crossover" => "off", "ipm_iteration_limit" => 3000, "time_limit" => 3600.0, "output_flag" => false)
optimizer = HVDCWISE_TEA.optimizer_with_attributes(HiGHS.Optimizer, setting_opt...)
setting = Dict("output" => Dict("branch_flows" => true, "duals" =>false), "conv_losses_mp" => false);

# Verify that Matlab or Octave is installed
matlab_octave_path = detect_matlab_or_octave()  # You can replace this line by the path of your Matlab/Octave launcher. for example: "C:/Users/n.barla/AppData/Local/Programs/GNU Octave/Octave-9.2.0/octave-launch.exe"

# Run the study
@time run_study(work_dir, hours_per_subsimulation, base_MVA, optimizer, setting, matlab_octave_path)
