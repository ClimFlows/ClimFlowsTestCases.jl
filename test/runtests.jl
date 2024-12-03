import ClimFlowsTestCases as TC
using Test

function check_case(Case::Type, F=Float32)
    case = Case(F; toto=1, u0=20)
    @info TC.describe(case)
    check_values(case, F)
    true
end

function check_values(case::TC.TestCaseSW, F)
    lon, lat = map(F, (0.0, 0.0))
    @info case(lon, lat)
end

function check_values(case::TC.TestCaseHPE, F)
    lon, lat, p = map(F, (0.0, 0.0, 1e5))
    @info TC.initial(case, lon, lat)
    @info case(lon, lat, p)
end

@testset "ClimFlowsTestCases.jl" begin
    @test check_case(TC.Williamson91{2})
    @test check_case(TC.Williamson91{6})
    @test check_case(TC.Jablonowski06)
    @test check_case(TC.DCMIP{42})
end
