using Test
include("../src/io/load_case.jl")
@testset "Tests Cases" begin
    @testset "Test case 01 : Single AC Line" begin
        # Load and run test case 1
        global results_1 = load_case("test01")

        if results_1 !== nothing
            @test results_1["termination_status"] == _HWTEA.OPTIMAL
            @test results_1["objective"] ≈ 30900.0 rtol = 1e-3
            # 8 MW generated and consumed on first timestep, no curtailment (8 MW demanded by load)
            @test results_1["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 8.0 atol=1e-3
            @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pt"] ≈ -8.0 atol=1e-3
            @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pf"] ≈ 8.0 atol=1e-3
            @test results_1["solution"]["nw"]["1"]["load"]["1"]["pflex"] ≈ 8.0 atol=1e-3
            @test results_1["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0.0 atol=1e-3
            # 10 MW generated and consumed on second timestep, 3 MW curtailed (13 MW demanded by load)
            @test results_1["solution"]["nw"]["2"]["gen"]["1"]["pg"] ≈ 10.0 atol=1e-3
            @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pt"] ≈ -10.0 atol=1e-3
            @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pf"] ≈ 10.0 atol=1e-3
            @test results_1["solution"]["nw"]["2"]["load"]["1"]["pflex"] ≈ 10.0 atol=1e-3
            @test results_1["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 3.0 atol=1e-3
        else
            println("Error loading test01. Check CSV and .m files name and location.")
        end
    end
    @testset "Test case 02 : Single AC Line with flexible load and non-dispatchable generators" begin
        # Load and run test case 2
        global results_2 = load_case("test02")

        if results_2 !== nothing
            @test results_2["termination_status"] == _HWTEA.OPTIMAL
            # @test results_2["objective"] ≈ 30900.0 rtol = 1e-3
            # 8 MW of load and 11 MW of available production. 1 MW should be shifted upward.
            # 9 MW of expected consumption & production. Neither reduction nor curtailment.
            @test results_2["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 9.0 atol=1e-3
            @test results_2["solution"]["nw"]["1"]["branch"]["1"]["pt"] ≈ -9.0 atol=1e-3
            @test results_2["solution"]["nw"]["1"]["branch"]["1"]["pf"] ≈ 9.0 atol=1e-3
            @test results_2["solution"]["nw"]["1"]["load"]["1"]["pflex"] ≈ 9.0 atol=1e-3
            # TODO add tests for pred, pshift_up, pshift_down
            @test results_2["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0.0 atol=1e-3
        else
            println("Error loading test02. Check CSV and .m files name and location.")
        end
    end
    @testset "Test case 04 : AC/DC Point to Point with Single AC Line" begin
        # Load and run test case 4
        global results_4 = load_case("test04")

        if results_4 !== nothing
            @test results_4["termination_status"] == _HWTEA.OPTIMAL
        else
            println("Error loading test04. Check CSV and .m files name and location.")
        end
    end
end