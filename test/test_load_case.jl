import HVDCWISE_TEA as _HWTEA

const _HWTEA_dir = dirname(dirname(pathof(_HWTEA))) # Root directory of HVDCWISE_TEA package
const save_results = true  # To save test results in Excel files

using Test
include("../src/io/load_case.jl")

@testset "Tests Cases" begin
    #
    @testset "Test case 01 : Single AC Line" begin
        # Load and run test case 1
        baseMVA = 100
        global results_1 = load_case("test01", baseMVA, true)

        if results_1 !== nothing
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
    #
    @testset "Test case 02 : Single AC Line with flexible load and non-dispatchable generators" begin
        # Load and run test case 2
        baseMVA = 100
        global results_2 = load_case("test02", baseMVA, true)

        if results_2 !== nothing
            @test results_2["termination_status"] == _HWTEA.OPTIMAL
            @test results_2["objective"] ≈ 411.0 rtol = 1e-3
            # t=1
            # 8 MW of load and 11 MW of available production. 1 MW should be shifted upward (load shift is limited to +/- 1 MW).
            # 9 MW of expected consumption & production. Neither load reduction nor load curtailment.
            @test results_2["solution"]["nw"]["1"]["gen"]["1"]["pg"] + results_2["solution"]["nw"]["1"]["gen"]["2"]["pg"] ≈ 9.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["1"]["branch"]["1"]["pt"] ≈ -9.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["1"]["branch"]["1"]["pf"] ≈ 9.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["1"]["load"]["1"]["pflex"] ≈ 9.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["1"]["load"]["1"]["pshift_up"] ≈ 1.0/baseMVA atol=1e-3
            # t=3
            # 8 MW of load and 5 MW of available production. 1 MW should be shifted downward (load shift is limited to +/- 1 MW).
            # 2 MW should be reduced (beyond 2 MW the load would be curtailed).
            # Expected consumption & production: 5 MW.
            @test results_2["solution"]["nw"]["3"]["gen"]["1"]["pg"] + results_2["solution"]["nw"]["3"]["gen"]["2"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["3"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["3"]["load"]["1"]["pshift_down"] ≈ 1.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["3"]["load"]["1"]["pred"] ≈ 2.0/baseMVA atol=1e-3
            # t=2
            # 8 MW of load and 7 MW of available production. Shifting downward would be compensated by load curtialment in hour 2, so no load shift and 1 MW of load reduction.
            # Expected consumption & production: 7 MW.
            @test results_2["solution"]["nw"]["2"]["gen"]["1"]["pg"] + results_2["solution"]["nw"]["2"]["gen"]["2"]["pg"] ≈ 7.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["2"]["load"]["1"]["pflex"] ≈ 7.0/baseMVA atol=1e-3
            @test results_2["solution"]["nw"]["2"]["load"]["1"]["pred"] ≈ 1.0/baseMVA atol=1e-3
        else
            println("Error loading test02. Check CSV and .m files name and location.")
        end
    end
    #
    @testset "Test case 04 : AC/DC Point to Point with Single AC Line" begin
        # Load and run test case 4
        baseMVA = 100
        global results_4 = load_case("test04", baseMVA, true)

        if results_4 !== nothing
            output_series_4 = results_4["solution"]["nw"]
            @test results_4["termination_status"] == _HWTEA.OPTIMAL
            # t=1. 8 MW of load and 10 MW of available production.
            # Both lines are fully available: there is no load shedding.
            @test output_series_4["1"]["gen"]["1"]["pg"] ≈ 8.0/baseMVA atol=1e-3
            @test output_series_4["1"]["load"]["1"]["pflex"] ≈ 8.0/baseMVA atol=1e-3
            @test output_series_4["1"]["load"]["1"]["pcurt"] ≈ 0.0/baseMVA atol=1e-3
            # t=2. 13 MW of load and 10 MW of available production.
            # Both lines are fully available: there is 3 MW of load shedding.
            @test output_series_4["2"]["gen"]["1"]["pg"] ≈ 10.0/baseMVA atol=1e-3
            @test output_series_4["2"]["load"]["1"]["pflex"] ≈ 10.0/baseMVA atol=1e-3
            @test output_series_4["2"]["load"]["1"]["pcurt"] ≈ 3.0/baseMVA atol=1e-3

            # t=3. 8 MW of load and 10 MW of available production.
            # The AC line is unavailable so only 5 MW can be consumed by the load.
            @test output_series_4["3"]["gen"]["1"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["3"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["3"]["load"]["1"]["pcurt"] ≈ 3.0/baseMVA atol=1e-3

            # t=4. 8 MW of load and 10 MW of available production.
            # 1 pole of the DC line is unavailable so its power rating is halved to 2.5 MW: only 7.5 MW can be consumed by the load.
            @test output_series_4["4"]["gen"]["1"]["pg"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["4"]["load"]["1"]["pflex"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["4"]["load"]["1"]["pcurt"] ≈ 0.5/baseMVA atol=1e-3

            # t=5. 8 MW of load and 10 MW of available production.
            # 1 pole of a converter is unavailable so its power rating is halved to 2.5 MW: only 7.5 MW can be consumed by the load.
            @test output_series_4["5"]["gen"]["1"]["pg"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["5"]["load"]["1"]["pflex"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["5"]["load"]["1"]["pcurt"] ≈ 0.5/baseMVA atol=1e-3

            # t=6. 8 MW of load and 10 MW of available production.
            # The negative pole of the DC line & the positive pole of a converter are unavailable so the whole DC transmission is unavailable: only 5 MW can be consumed by the load.
            @test output_series_4["6"]["gen"]["1"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["6"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["6"]["load"]["1"]["pcurt"] ≈ 3.0/baseMVA atol=1e-3
        else
            println("Error loading test04. Check CSV and .m files name and location.")
        end
    end
    #
    @testset "Test case 05 : 1 storage + 1 generator + 1 load" begin
        # Load and run test case 5
        baseMVA = 100
        global results_5 = load_case("test05", baseMVA, true)

        if results_5 !== nothing
            @test results_5["termination_status"] == _HWTEA.OPTIMAL
            @test results_5["objective"] ≈ 10450 rtol = 1e-3
            # t=1: Demand=8, Available=10 and storage is available -> 1 MW consumed by the storage & 8 MW by the load (pgcurt=1 MW)
            @test results_5["solution"]["nw"]["1"]["gen"]["1"]["pg"] ≈ 9/baseMVA atol=1e-3
            # @test results_5["solution"]["nw"]["1"]["gen"]["1"]["pgcurt"] ≈ 1/baseMVA atol=1e-3  # no pgcurt because dispatchable generator
            @test results_5["solution"]["nw"]["1"]["load"]["1"]["pcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["1"]["storage"]["1"]["ps"] ≈ -1/baseMVA atol=1e-3  # ps=sc-sd (storage power = chargin power - discharging power)
            @test results_5["solution"]["nw"]["1"]["storage"]["1"]["sd"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd

            # t=2: Demand=12, Available=10 but storage not available -> 2 MW of load shedding
            @test results_5["solution"]["nw"]["2"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            # @test results_5["solution"]["nw"]["2"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 2/baseMVA atol=1e-3
            @test !haskey(results_5["solution"]["nw"]["2"], "storage")  # no storage in the results because status=0

            # t=3: Demand=12, Available=10 and storage is available -> 1 MW produced by the storage & 1 MW of load shedding
            @test results_5["solution"]["nw"]["3"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            # @test results_5["solution"]["nw"]["3"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["3"]["load"]["1"]["pcurt"] ≈ 1/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["3"]["storage"]["1"]["ps"] ≈ 1/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["1"]["storage"]["1"]["sc"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd
        else
            println("Error loading test05. Check CSV and .m files name and location.")
        end
    end
    #
end