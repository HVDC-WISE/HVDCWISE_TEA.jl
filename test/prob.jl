# Test problem functions

# The purpose of the tests contained in this file is to detect if anything has accidentally
# changed in the problem functions. Accordingly, only termination status and objective value
# are tested. For testing specific features, it is better to write ad-hoc tests in separate files.


@testset "Problem" begin

    @testset "mcdcopf DCP" begin
        @testset "case5_2grids_MC" begin

            file = joinpath(_HWTEA_dir, "test/data/grids/acdc/case5_3grids_MC.m")
            result_dcp = _HWTEA.solve_mc_acdcopf(file, _PM.DCPPowerModel, lp_optimizer)
            @test result_dcp["termination_status"] == _HWTEA.OPTIMAL
            @test result_dcp["objective"] ≈ 823.0 rtol = 1e-3
        end
    end
end
