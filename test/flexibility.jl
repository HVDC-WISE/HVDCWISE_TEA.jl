# Test slack resources of the optimisation model

# The purpose of these tests is to verify if slack resources, such as load
# and non-dispatchable generation, relax their working point in the correct way 
# if system constraints becomes active. 

## Load model

@testset "Load model" begin

    @testset "Load curtailment" begin
        
        """
        Case where no load offers flexibility and, thus, only forced load curtailment is possible.
        The AC load at the AC/DC converter of area 3 is set such that it exceeds sum of the generation
        capacity of the area (i.e. 300 MW) and of the HVDC interconnection (i.e. 50 MW). Both reach
        saturation. A single optimisation period is considered.
        """
        # Read and process data
        file = joinpath(_HWTEA_dir, "test/data/case5/case5_3grids_MC - curtailment.m")
        data = _HWTEA.parse_data(file)
        # Increase demand at load of area 3
        data["load"]["9"]["pd"] = 400 / data["baseMVA"]
        # Solve optimisation problem
        result = _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, lp_optimizer, setting = s)

        @test result["solution"]["load"]["9"]["pflex"]      ≈ 3.488     rtol=1e-3
        @test result["solution"]["load"]["9"]["pcurt"]      ≈ 0.512     atol=1e-3
        @test result["solution"]["gen"]["5"]["pg"]          ≈ 3.0       rtol=1e-3
        @test result["solution"]["convdc"]["3"]["pdc"][1]   ≈ 0.5       atol=1e-3
        @test result["objective"]                           ≈ 156410.4  rtol=1e-3    
    end
end

## Non-dispatchable generation model

@testset "Generator model" begin

    @testset "Renewable generation curtailment" begin

        """
        Case in which a large renewable plant is connected at area 2, provinding enough power
        to satisfy the entire demand of area 2 (i.e. 165 MW) and to saturate the 100 MVA AC/DC
        converter connected to area 2.
        """
        # Read and process data
        file = joinpath(_HWTEA_dir, "test/data/case5/case5_3grids_MC - curtailment.m")
        data = _HWTEA.parse_data(file)
        # Increase renewable generation at the AC terminal of converter 2
        data["gen"]["6"]["pmax"] = 300 / data["baseMVA"]
        # Solve optimisation problem
        result = _HWTEA.solve_mc_acdcopf(data, _PM.DCPPowerModel, lp_optimizer, setting = s)
        # Get index of load connected area 2
        mask = _PM.component_table(data, "load", ["load_bus"])[:, 2]
        mask = _PM.component_table(data, "bus", ["area"])[mask, 2] .== 2

        @test sum(_PM.component_table(result["solution"], "load", ["pflex"])[mask, 2])  ≈ 1.65      rtol=1e-3
        @test result["solution"]["gen"]["3"]["pg"]                                      ≈ 0.1       atol=1e-3
        @test result["solution"]["gen"]["4"]["pg"]                                      ≈ 0.1       atol=1e-3
        @test result["solution"]["gen"]["6"]["pg"]                                      ≈ 2.450     rtol=1e-3
        @test result["solution"]["gen"]["6"]["pgcurt"]                                  ≈ 0.550     atol=1e-3
        @test sum(result["solution"]["convdc"]["2"]["pgrid"])                           ≈ 1.0       atol=1e-3
        @test result["objective"]                                                       ≈ 55273.6   rtol=1e-3 
    end
end