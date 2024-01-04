include("../src/io/load_case.jl")

# Load and run test case 1
results_1 = load_case("Test01")

if results_1 !== nothing
    println("---------------- RESULTS TEST CASE 1 ----------------")
    for network_id in 1:length(results_1["solution"]["nw"])
        print("timestep : $network_id"  )
        _PM.print_summary(results_1["solution"]["nw"]["$network_id"])
    end
else
    println("Error loading Test01. Check CSV and .m files location.")
end

# Load and run test case 4
results_4 = load_case("Test04")

if results_4 !== nothing
    println("---------------- RESULTS TEST CASE 4 ----------------")
    for network_id in 1:length(results_4["solution"]["nw"])
        print("timestep : $network_id"  )
        _PM.print_summary(results_4["solution"]["nw"]["$network_id"])
    end
else
    println("Error loading Test04. Check CSV and .m files location.")
end