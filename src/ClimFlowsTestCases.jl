module ClimFlowsTestCases

"""
    abstract type TestCase end

Parent type for test case types.
"""
abstract type TestCase end
# @inline Base.getproperty(case::TestCase, name) = getproperty(getfield(case,:params), name)

"""
Parent type for shallow-water test cases.
"""
abstract type TestCaseSW <: TestCase end

"""
Parent type for hydrostatic test cases.
"""
abstract type TestCaseHPE <: TestCase end

"""
    params = default_params(TestCaseType)
Returns a named tuple of default parameters for a certain test case.
More about named tuples [here](https://stackoverflow.com/questions/60883704/how-to-manipulate-named-tuples).
"""
function default_params() end

"""
    case = default_testcase(Case::Type, F=Float64)
    w91_6 = default_testcase(Williamson91{6}, Float32)

Returns `case`, the test case with default parameters for case `Case`.
"""
function default_testcase(Case::Type{TC}, F::Type{FF}=Float64) where { TC<:TestCase, FF<:Real }
    Case(map(F, default_params(Case)))
end

"""
    ulon, ulat, gH = initial_flow(lon, lat, case::TestCaseSW)

For a shallow-water test case, returns the initial geopotential thickness and velocity
at given latitude and longitude.

    gH, ulon, ulat, q = initial_flow(lon, lat, p, case::TestCaseHPE)

For a hydrostatic test case, returns the initial geopotential, velocity
and composition at given latitude, longitude and pressure.
"""
function initial_flow end

"""
    ps, Phis = initial_surface(lon, lat, case::TestCaseHPE)

For a hydrostatic test case, returns the initial surface pressure and geopotential
at given latitude and longitude.
"""
function initial_surface end

"""
    string = describe(case::TestCase)

Plain text description of `case`, possibly multi-line.
"""
function describe end

Base.show(io::IO, case::TestCase) = print(io, describe(case))

# Implementation details that may change between non-breaking versions.
module priv

using ..ClimFlowsTestCases: TestCase, TestCaseSW, TestCaseHPE
import ..ClimFlowsTestCases: default_params, describe, initial_flow, initial_surface

# We make heavy use of the destructuring syntax
#       (; a, b, c) = x
# equivalent to
#       a, b, c = x.a, x.b, x.c

include("julia/williamson91.jl")
include("julia/DCMIP2012.jl")
include("julia/Jablonowski06.jl")

end # priv

using .priv

# it would be simpler to export these types from priv, but it would not make them visible via TAB completion
"""
    struct Williamson91{N} <: TestCaseSW

Williamson (1991) test cases. Currently only N=6 is implemented.
"""
const Williamson91 = priv.Williamson91

"""
    struct DCMIP{N} <: TestCaseHPE

DCMIP 2012 test cases. Currently only N=42 (test case 4.2) is implemented.
"""
const DCMIP = priv.DCMIP


"""
    struct Jablonowski06{P} <: TestCaseHPE

Jablonowski (2006) and DCMIP 4.2 test cases (dry vs. moist). The only difference is specific humidity is included for the moist case.
"""
const Jablonowski06 = priv.Jablonowski06

end # module
