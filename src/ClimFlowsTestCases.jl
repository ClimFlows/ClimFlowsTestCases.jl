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

function default_testcase(Case::Type{TC}, F::Type{FF}=Float64) where { TC<:TestCase, FF<:Real }
    Case(map(F, default_params(Case)))
end

function initial_flow end
function initial_surface end
function describe end

Base.show(io::IO, case::TestCase) = print(io, describe(case))

"""
    module priv

This module contains implementation detail that may change between non-breaking versions.
"""
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
const Williamson91 = priv.Williamson91
const DCMIP = priv.DCMIP
const Jablonowski06 = priv.Jablonowski06

end # module
