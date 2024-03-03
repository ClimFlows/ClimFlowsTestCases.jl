import ClimFlowsTestCases as TC
using Test

function check_case(Case::Type, F=Float32)
    case = TC.default_testcase(Case, F)
    @info TC.describe(case)
    check_values(case, F)
    true
end

function check_values(case::TC.TestCaseSW, F)
    lon, lat = map(F, (0.0, 0.0))
    @info TC.initial_flow(lon, lat, case)
end

function check_values(case::TC.TestCaseHPE, F)
    lon, lat, p = map(F, (0.0, 0.0, 1e5))
    @info TC.initial_surface(lon, lat, case)
    @info TC.initial_flow(lon, lat, p, case)
end

@testset "ClimFlowsTestCases.jl" begin
    @test check_case(TC.Williamson91{6})
    @test check_case(TC.Jablonowski06)
    @test check_case(TC.DCMIP{42})
end
