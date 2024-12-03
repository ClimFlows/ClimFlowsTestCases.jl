module ClimFlowsTestCases

"""
    abstract type TestCase end

Parent type for test case types. Derived types must implement a constructor
accepting an optional precision and keyword arguments to override default ones. For example:

    case = Williams91{6}(Float32; R0=3e-3)

Derived types must also possess a field `params` which is a named tuple containing the parameters, and implement
[`initial`](@ref) and [`describe`](@ref). 
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
    gh, ulon, ulat = initial(case::TestCaseSW, lon, lat)

For a shallow-water test case, returns the initial geopotential thickness and velocity
at given latitude and longitude.

    ps, Phis = initial(case::TestCaseHPE, lon, lat)

For a hydrostatic test case, returns the initial surface pressure and geopotential
at given latitude and longitude.

    Phi, ulon, ulat, q = initial(case::TestCaseHPE, lon, lat, p)

For a hydrostatic test case, returns the initial geopotential, velocity
and composition at given latitude, longitude and pressure.
"""
function initial end

"""
    string = describe(case::TestCase)

Plain text description of `case`, possibly multi-line.
"""
function describe end

Base.show(io::IO, case::TestCase) = print(io, describe(case))
(case::TestCase)(lon, lat) = initial(case, lon, lat)
(case::TestCaseHPE)(lon, lat, p) = initial(case, lon, lat, p)

# Implementation details that may change between non-breaking versions.
module priv

using ..ClimFlowsTestCases: TestCase, TestCaseSW, TestCaseHPE
import ..ClimFlowsTestCases: describe, initial

# We make heavy use of the destructuring syntax
#       (; a, b, c) = x
# equivalent to
#       a, b, c = x.a, x.b, x.c

include("julia/williamson91.jl")
include("julia/DCMIP2012.jl")
include("julia/Jablonowski06.jl")
include("julia/isothermal.jl")

override(F, defaults, params) = map(F, merge(defaults, params)[propertynames(defaults)])

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

"""
    struct Isothermal <: TestCaseHPE
    
Return `testcase` corresponding to an initial state of rest with exponential decrease of pressure with geopotential.
For a perfect gas, this coincides with the classic isothermal profile, and psurf * vsurf = R*T with T the constant temperature.
The surface geopotential defines a wiggly mountain. Numerics based on terrain-following coordinates typically struggle with such case
and generate spurious motion.
"""
const Isothermal = priv.Isothermal

end # module
