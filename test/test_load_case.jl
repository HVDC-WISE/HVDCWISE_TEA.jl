using Distributed
using HVDCWISE_TEA
import HVDCWISE_TEA as _HWTEA
import PowerModels as _PM

using Test

const save_results = true  # To save test results in Excel files


@testset "Tests Cases" begin
    #
    @testset "Test case 01 : Single AC Line" begin
        # Load and run test case 1
        baseMVA = 100
        global results_1 = load_case("test01", baseMVA, true)

        if results_1 !== nothing
            @test results_1["termination_status"] == _HWTEA.OPTIMAL
            @test results_1["objective"] ≈ 9900.0 rtol = 1e-3
            output_series_1 = results_1["solution"]["nw"]
            # All results are in per unit, with baseMVA = 100 MVA (specified in the .m file)
            # 8 MW generated and consumed on first timestep, no curtailment (8 MW demanded by load)
            @test output_series_1["1"]["gen"]["1"]["pg"] ≈ 8/baseMVA atol=1e-3
            @test output_series_1["1"]["branch"]["1"]["pt"] ≈ -8/baseMVA atol=1e-3
            @test output_series_1["1"]["branch"]["1"]["pf"] ≈ 8/baseMVA atol=1e-3
            @test output_series_1["1"]["load"]["1"]["pflex"] ≈ 8/baseMVA atol=1e-3
            @test output_series_1["1"]["load"]["1"]["pred"] ≈ 0.0 atol=1e-3
            # 10 MW generated and consumed on second timestep, 3 MW curtailed (13 MW demanded by load)
            @test output_series_1["2"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            @test output_series_1["2"]["branch"]["1"]["pt"] ≈ -10/baseMVA atol=1e-3
            @test output_series_1["2"]["branch"]["1"]["pf"] ≈ 10/baseMVA atol=1e-3
            @test output_series_1["2"]["load"]["1"]["pflex"] ≈ 10/baseMVA atol=1e-3
            @test output_series_1["2"]["load"]["1"]["pred"] ≈ 3/baseMVA atol=1e-3
        else
            println("Error loading test01. Check CSV and .m files name and location.")
        end
    end
    #
    @testset "Test case 02 : Single AC Line with flexible load and non-dispatchable generator" begin
        # Load and run test case 2
        baseMVA = 100
        global results_2 = load_case("test02", baseMVA, true)

        if results_2 !== nothing
            @test results_2["termination_status"] == _HWTEA.OPTIMAL
            @test results_2["objective"] ≈ 411.0 rtol = 1e-3
            output_series_2 = results_2["solution"]["nw"]
            # t=1
            # 8 MW of load and 11 MW of available production. 1 MW should be shifted upward (load shift is limited to +/- 1 MW).
            # 9 MW of expected consumption & production. Neither load reduction nor load curtailment.
            @test output_series_2["1"]["gen"]["1"]["pg"] + output_series_2["1"]["gen"]["2"]["pg"] ≈ 9.0/baseMVA atol=1e-3
            @test output_series_2["1"]["branch"]["1"]["pt"] ≈ -9.0/baseMVA atol=1e-3
            @test output_series_2["1"]["branch"]["1"]["pf"] ≈ 9.0/baseMVA atol=1e-3
            @test output_series_2["1"]["load"]["1"]["pflex"] ≈ 9.0/baseMVA atol=1e-3
            @test output_series_2["1"]["load"]["1"]["pshift_up"] ≈ 1.0/baseMVA atol=1e-3
            # t=3
            # 8 MW of load and 5 MW of available production. 1 MW should be shifted downward (load shift is limited to +/- 1 MW).
            # 2 MW should be reduced (beyond 2 MW the load would be curtailed).
            # Expected consumption & production: 5 MW.
            @test output_series_2["3"]["gen"]["1"]["pg"] + output_series_2["3"]["gen"]["2"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_2["3"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_2["3"]["load"]["1"]["pshift_down"] ≈ 1.0/baseMVA atol=1e-3
            @test output_series_2["3"]["load"]["1"]["pred"] ≈ 2.0/baseMVA atol=1e-3
            # t=2
            # 8 MW of load and 7 MW of available production. Shifting downward would be compensated by load curtialment in hour 2, so no load shift and 1 MW of load reduction.
            # Expected consumption & production: 7 MW.
            @test output_series_2["2"]["gen"]["1"]["pg"] + output_series_2["2"]["gen"]["2"]["pg"] ≈ 7.0/baseMVA atol=1e-3
            @test output_series_2["2"]["load"]["1"]["pflex"] ≈ 7.0/baseMVA atol=1e-3
            @test output_series_2["2"]["load"]["1"]["pred"] ≈ 1.0/baseMVA atol=1e-3
        else
            println("Error loading test02. Check CSV and .m files name and location.")
        end
    end
    #
    @testset "Test case 04 : AC/DC Point to Point + Single AC Line" begin
        # Load and run test case 4
        baseMVA = 100
        global results_4 = load_case("test04", baseMVA, true)

        if results_4 !== nothing
            @test results_4["termination_status"] == _HWTEA.OPTIMAL
            @test results_4["objective"] ≈ 32150.0 rtol = 1e-3
            output_series_4 = results_4["solution"]["nw"]
            # t=1. 8 MW of load and 10 MW of available production.
            # Both lines are fully available: there is no load shedding.
            @test output_series_4["1"]["gen"]["1"]["pg"] ≈ 8.0/baseMVA atol=1e-3
            @test output_series_4["1"]["load"]["1"]["pflex"] ≈ 8.0/baseMVA atol=1e-3
            @test output_series_4["1"]["load"]["1"]["pred"] ≈ 0.0/baseMVA atol=1e-3
            # t=2. 13 MW of load and 10 MW of available production.
            # Both lines are fully available: there is 3 MW of load shedding.
            @test output_series_4["2"]["gen"]["1"]["pg"] ≈ 10.0/baseMVA atol=1e-3
            @test output_series_4["2"]["load"]["1"]["pflex"] ≈ 10.0/baseMVA atol=1e-3
            @test output_series_4["2"]["load"]["1"]["pred"] ≈ 3.0/baseMVA atol=1e-3

            # t=3. 8 MW of load and 10 MW of available production.
            # The AC line is unavailable so only 5 MW can be consumed by the load.
            @test output_series_4["3"]["gen"]["1"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["3"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["3"]["load"]["1"]["pred"] ≈ 3.0/baseMVA atol=1e-3

            # t=4. 8 MW of load and 10 MW of available production.
            # 1 pole of the DC line is unavailable so its power rating is halved to 2.5 MW: only 7.5 MW can be consumed by the load.
            @test output_series_4["4"]["gen"]["1"]["pg"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["4"]["load"]["1"]["pflex"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["4"]["load"]["1"]["pred"] ≈ 0.5/baseMVA atol=1e-3

            # t=5. 8 MW of load and 10 MW of available production.
            # 1 pole of a converter is unavailable so its power rating is halved to 2.5 MW: only 7.5 MW can be consumed by the load.
            @test output_series_4["5"]["gen"]["1"]["pg"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["5"]["load"]["1"]["pflex"] ≈ 7.5/baseMVA atol=1e-3
            @test output_series_4["5"]["load"]["1"]["pred"] ≈ 0.5/baseMVA atol=1e-3

            # t=6. 8 MW of load and 10 MW of available production.
            # The negative pole of the DC line & the positive pole of a converter are unavailable so the whole DC transmission is unavailable: only 5 MW can be consumed by the load.
            @test output_series_4["6"]["gen"]["1"]["pg"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["6"]["load"]["1"]["pflex"] ≈ 5.0/baseMVA atol=1e-3
            @test output_series_4["6"]["load"]["1"]["pred"] ≈ 3.0/baseMVA atol=1e-3
        else
            println("Error loading test04. Check CSV and .m files name and location.")
        end
    end
    #=
    @testset "Test case 05 : 1 storage + 1 generator + 1 load" begin
        # Load and run test case 5
        baseMVA = 100
        global results_5 = load_case("test05", baseMVA, true)

        if results_5 !== nothing
            @test results_5["termination_status"] == _HWTEA.OPTIMAL
            @test results_5["objective"] ≈ 10450 rtol = 1e-3
            output_series_5 = results_5["solution"]["nw"]
            # t=1: Demand=8, Available=10 and storage is available -> 1 MW consumed by the storage & 8 MW by the load (pgcurt=1 MW)
            @test output_series_5["1"]["gen"]["1"]["pg"] ≈ 9/baseMVA atol=1e-3
            # @test output_series_5["1"]["gen"]["1"]["pgcurt"] ≈ 1/baseMVA atol=1e-3  # no pgcurt because dispatchable generator
            @test output_series_5["1"]["load"]["1"]["pcurt"] ≈ 0/baseMVA atol=1e-3
            @test output_series_5["1"]["storage"]["1"]["ps"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd (storage power = chargin power - discharging power)
            @test output_series_5["1"]["storage"]["1"]["sc"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd

            # t=2: Demand=12, Available=10 but storage not available -> 2 MW of load shedding
            @test output_series_5["2"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            # @test output_series_5["2"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test results_5["solution"]["nw"]["2"]["load"]["1"]["pcurt"] ≈ 2/baseMVA atol=1e-3
            @test !haskey(output_series_5["2"], "storage")  # no storage in the results because status=0

            # t=3: Demand=12, Available=10 and storage is available -> 1 MW produced by the storage & 1 MW of load shedding
            @test output_series_5["3"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            # @test output_series_5["3"]["gen"]["1"]["pgcurt"] ≈ 0/baseMVA atol=1e-3
            @test output_series_5["3"]["load"]["1"]["pcurt"] ≈ 1/baseMVA atol=1e-3
            @test output_series_5["3"]["storage"]["1"]["ps"] ≈ -1/baseMVA atol=1e-3
            @test output_series_5["1"]["storage"]["1"]["sd"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd
        else
            println("Error loading test05. Check CSV and .m files name and location.")
        end
    end
    =#
    @testset "Test case 05b : 1 storage + 1 generator + 1 load" begin
        # Load and run test case 5
        baseMVA = 100
        global results_5 = load_case("test05b", baseMVA, true)

        if results_5 !== nothing
            @test results_5["termination_status"] == _HWTEA.OPTIMAL
            @test results_5["objective"] ≈ 4500 rtol = 1e-2  # 4508 instead of 4500 due to losses
            output_series_5 = results_5["solution"]["nw"]
            for t in ["1", "2"]
            # t = 1 & 2: Demand=9 and Available=10 -> 1 MW consumed by the storage & 9 MW by the load (no load shedding)
                @test output_series_5[t]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
                @test output_series_5[t]["load"]["1"]["pred"] ≈ 0/baseMVA atol=1e-3
                @test output_series_5[t]["storage"]["1"]["ps"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd (storage power = chargin power - discharging power)
                # @test output_series_5[t]["storage"]["1"]["sc"] ≈ 1/baseMVA atol=1e-3  # ps=sc-sd
            end

            # t=3: Demand=13 and Available=10 -> 2 MW produced by the storage (because its production power is limited to 2 MW) & 1 MW of load shedding
            @test output_series_5["3"]["gen"]["1"]["pg"] ≈ 10/baseMVA atol=1e-3
            @test output_series_5["3"]["load"]["1"]["pred"] ≈ 1/baseMVA atol=1e-3
            @test output_series_5["3"]["storage"]["1"]["ps"] ≈ -2/baseMVA atol=1e-3
            # @test output_series_5["1"]["storage"]["1"]["sd"] ≈ 2/baseMVA atol=1e-3  # ps=sc-sd
        else
            println("Error loading test05b. Check CSV and .m files name and location.")
        end
    end
    #=
    @testset "Test size" begin
        # Load and run test case
        baseMVA = 100
        N = 100
        S = 100

        results_size = load_case("test_size", baseMVA, true)

        if results_size !== nothing
            @test results_size["termination_status"] == _HWTEA.OPTIMAL
            @test results_size["objective"] ≈ 4500*N*S rtol = 1e-3
            output_series_size = results_size["solution"]["nw"]
            for s=1:S
                # t=3*S-2: Demand=9*N & Available=10*N -> N MW consumed by the storage & 9*N MW by the load (pgcurt=0 MW)
                @test sum([output_series_size["$(3*s-2)"]["load"]["$n"]["pflex"] for n=1:N])/N ≈ 9/baseMVA atol=1e-3  # FIXME
                @test sum([output_series_size["$(3*s-2)"]["storage"]["$n"]["ps"] for n=1:N])/N ≈ 1/baseMVA atol=1e-3  # FIXME
                # t=3*S: Demand=13*N & Available=10*N -> 2*N MW produced by the storage, N MW of load shedding
                @test sum([output_series_size["$(3*s)"]["load"]["$n"]["pflex"] for n=1:N])/N ≈ 12/baseMVA atol=1e-3  # FIXME
                @test sum([output_series_size["$(3*s)"]["storage"]["$n"]["ps"] for n=1:N])/N ≈ -2/baseMVA atol=1e-3  # FIXME
            end
        else
            println("Error loading test_size. Check CSV and .m files name and location.")
        end
    end
    =#
end