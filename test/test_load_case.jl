import HVDCWISE_TEA as _HWTEA

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package

using Test
include("../src/io/load_case.jl")
@testset "Tests Cases" begin
    #
    @testset "Test case 01 : Single AC Line" begin
        # Load and run test case 1
        global results_1 = load_case("test01")

        if results_1 !== nothing
            baseMVA = 100
            @test results_1["termination_status"] == _HWTEA.OPTIMAL
            @test results_1["objective"] ≈ 9900.0 rtol = 1e-3
            # All results are in per unit, with baseMVA = 100 MVA (specified in the .m file)
            # 8 MW generated and consumed on first timestep, no curtailment (8 MW demanded by load)
            @test results_1["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 8/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pt"] ≈ -8/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["1"]["branch"]["1"]["pf"] ≈ 8/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["1"]["load"]["1"]["pflex"] ≈ 8/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0.0 atol=1e-3
            # 10 MW generated and consumed on second timestep, 3 MW curtailed (13 MW demanded by load)
            @test results_1["solution"]["nw"]["2"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pt"] ≈ -10/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["2"]["branch"]["1"]["pf"] ≈ 10/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["2"]["load"]["1"]["pflex"] ≈ 10/baseMVA atol=1e-3
            @test results_1["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 3/baseMVA atol=1e-3
        else
            println("Error loading test01. Check CSV and .m files name and location.")
        end
    end
    #=
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
    =#
    @testset "Test case 05a : 1 storage + 1 generator + 1 load" begin
        # Load and run test case 5a
        global results_5 = load_case("test05a")

        ## baseMVA reading
        file = joinpath(_HWTEA_dir, "test\\data\\test05a\\.m_files\\test05a.m")
        file_string_data = read(open(file), String)
        file_dict_data = _IM.parse_matlab_string(file_string_data, extended=true)[1];
        baseMVA = file_dict_data["mpc.baseMVA"]  # 100 MVA

        if results_5 !== nothing
            @test results_5["termination_status"] == _HWTEA.OPTIMAL
            @test results_5["objective"] ≈ 10450 rtol = 1e-3
            # t=1: Demand=8, Available=10 and storage is available -> 1 MW consumed by the storage & 8 MW by the load (pgcurt=1 MW)
            @test results_5["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 9/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["1"]["gen"]["1"]["pgcurt"] ≈ 1/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["1"]["storage"]["1"]["ps"] ≈ -1/baseMVA atol=1e-3  # FIXME what is the storage attribute for storage net production ?

            # t=2: Demand=12, Available=10 but storage not available -> 2 MW of load shedding
            @test results_5["solution"]["nw"]["2"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["2"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 2/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["2"]["storage"]["1"]["ps"] ≈ 0/baseMVA atol=1e-3

            # t=3: Demand=12, Available=10 and storage is available -> 1 MW produced by the storage & 1 MW of load shedding
            @test results_5["solution"]["nw"]["3"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["3"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["3"]["load"]["1"]["pcurt"] ≈ 1/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["3"]["storage"]["1"]["ps"] ≈ 1/baseMVA atol=1e-3
        else
            println("Error loading test05a. Check CSV and .m files name and location.")
        end
    end
end