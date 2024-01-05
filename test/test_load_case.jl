using Test
include("../src/io/load_case.jl")

@testset "Test case 01 : Single AC Line" begin
    # Load and run test case 1
    results_1 = load_case("test01")

    if results_1 !== nothing
        @test results_1["termination_status"] == _HWTEA.OPTIMAL
        @test results_1["objective"] ≈ 30900.0 rtol = 1e-3
        # 8 MW generated and consumed on first timestep, no curtailment (8 MW demanded by load)
        @test results_1["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 8.0 atol=1e-3
        @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pt"] ≈ -8.0 atol=1e-3
        @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pf"] ≈ 8.0 atol=1e-3
        @test results_1["solution"]["nw"]["1"]["load"]["1"]["pflex"] ≈ 8.0 atol=1e-3
        @test results_1["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0.0 atol=1e-3
        # 10 MW generated and consumed on second timestep, 3MW curtailed (13 MW demanded by load)
        @test results_1["solution"]["nw"]["2"]["gen"]["1"]["pg"] ≈ 10.0 atol=1e-3
        @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pt"] ≈ -10.0 atol=1e-3
        @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pf"] ≈ 10.0 atol=1e-3
        @test results_1["solution"]["nw"]["2"]["load"]["1"]["pflex"] ≈ 10.0 atol=1e-3
        @test results_1["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 3.0 atol=1e-3
    else
        println("Error loading test01. Check CSV and .m files name and location.")
    end
end

@testset "Test case 04 : AC/DC Point to Point with Single AC Line" begin
    # Load and run test case 4
    results_4 = load_case("test04")

    if results_4 !== nothing
        println("---------------- RESULTS TEST CASE 4 ----------------")
        for network_id in 1:length(results_4["solution"]["nw"])
            print("timestep : $network_id"  )
            _PM.print_summary(results_4["solution"]["nw"]["$network_id"])
        end
    else
        println("Error loading test04. Check CSV and .m files name and location.")
    end
end