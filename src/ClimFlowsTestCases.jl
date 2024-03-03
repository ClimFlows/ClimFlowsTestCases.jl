module ClimFlowsTestCases

"""
This module exports into ClimFlowTestCases those functions and types that are part of the API.
Anything else is implementation detail and may change between non-breaking version.
"""
module priv

export TestCase, TestCaseSW, TestCaseHPE
export Williamson91, DCMIP, Jablonowski06
export default_params, default_testcase, describe, initial_flow, initial_surface

"""
    abstract type TestCase end

Parent type for test case types.
"""
abstract type TestCase end
# @inline Base.getproperty(case::TestCase, name) = getproperty(getfield(case,:params), name)

abstract type TestCaseSW <: TestCase end
abstract type TestCaseHPE <: TestCase end

"""
    params = default_params(TestCaseType)
Returns a named tuple of default parameters for a certain test case.
More about named tuples [here](https://stackoverflow.com/questions/60883704/how-to-manipulate-named-tuples).
"""
function default_params() end

# We make heavy use of the destructuring syntax
#       (; a, b, c) = x
# equivalent to
#       a, b, c = x.a, x.b, x.c

function default_testcase(Case::Type{TC}, F::Type{FF}) where { TC<:TestCase, FF<:Real }
    Case(map(F, default_params(Case)))
end

#============== Williamson (1991) test cases ===============#

"""
    Williamson91{N}
Type for Williamson (1991) test cases. Currently only N=6 is implemented.
"""
struct Williamson91{N,P} <: TestCase
    params::P
    # The syntax params[(...)] extracts only the parameters relevant for the test case
    Williamson91{6}(params) = let p=params[(:R0, :Omega, :gH0, :K, :n)]
        new{6, typeof(p)}(p)
    end
    Williamson91{2}(params) = let p=params[(:R0, :Omega, :gH0, :u0, :alpha)]
        new{2, typeof(p)}(p)
    end
end

W91_6 = Williamson91{6}
default_params(::Type{W91_6}) = (R0=6.4e6, Omega=2*pi/86400., gH0=1.0e4, K=7.848e-6, n=4)

function describe(case::W91_6)
    (; n, K, Omega) = case.params
    nu = (n*(n+3)*K-2Omega)/((1+n)*(2+n))
    println("Williamson 1991 test case number 6 (Rossby-Haurwitz wave).")
    println("Angular velocity of Rossby-Haurwitz wave : $nu")
    println("Time to advance by 1/2 period : $(pi/nu/n/3600)h")
end

function initial_flow(lon, lat, case::W91_6)
    (; n, K, Omega, R0, gH0) = case.params
    u = R0*K*( cos(lat) + (cos(lat)^(n-1))*(n*sin(lat)^2 - cos(lat)^2) *cos(n*lon) )
    v = -R0*K*n*(cos(lat)^(n-1))*sin(lat)*sin(n*lon)
    A = 0.5*K*(2Omega + K)*(cos(lat))^2 + 0.25*(K^2)*((cos(lat))^(2*n))*((n+1)*cos(lat)^2
        + (2*n^2 - n - 2) - 2*n^2*cos(lat)^(-2) )
    B = ( (2*(Omega + K)*K )/( (n+1)*(n+2)) ) *cos(lat)^n * ( (n^2 + 2*n + 2) - (n+1)^2*cos(lat)^2 )
    C = 0.25*(K^2)*(cos(lat)^(2*n))*((n+1)*cos(lat)^2 - (n+2) )
    gH = gH0 + (R0^2)*(A  + B*cos(n*lon) + C*cos(2*n*lon))

    return u, v, gH
end

const W91_2 = Williamson91{2}
default_params(::Type{W91_2}) = (R0=6.4e6, Omega=2*pi/86400., gH0=2.94e4, u0=2*pi*6.4e6/(12*86400), alpha=0.0)

function describe(case::W91_2)
    println("Williamson 1991 test case number 2 (Global steady state nonlinear geostrophic flow).")
end

function initial_flow(lon, lat, case::W91_2)
    (; R0, Omega, gH0, u0, alpha) = case.params
    u = u0*( cos(lat)*cos(alpha) + cos(lon)*sin(lat)*sin(alpha) )
    v = -u0*sin(lon)*sin(alpha)
    gH = gH0 - ( R0*Omega*u0 + u0^2/2 )*( sin(lat)*cos(alpha) - cos(lon)*cos(lat)*sin(alpha) )^2
    return u, v, gH
end

#================= DCMIP4.2 : (optionnaly sheared) solid-body rotation ===============#

struct DCMIP{N,P} <: TestCase
    params :: P
    DCMIP{42}(params) = let p=params[(:p0, :u0, :pv0, :c)]
        new{42, typeof(p)}(p)
    end
end

const DCMIP42 = DCMIP{42}
function default_params(::Type{DCMIP42})
    T0, Rd = 300.0, 287.0
    return (p0=1e5, u0=20.0, pv0=Rd*T0, c=0.0)
end

function describe(case::DCMIP42)
    (; p0, u0, pv0, c) = case.params
    println("DCMIP 4.2 test case with reference values p=$p0, u=$u0, rho=$(p0/pv0), du/dPhi=$(c*u0)")
end

function initial_surface(lon, lat, case::DCMIP42)
    p0 = case.params.p0
    Phi, _, _, _ = initial_flow(lon, lat, p0, case)
    return p0, Phi
end

function initial_flow(lon, lat, p, case::DCMIP42)
    (; p0, u0, pv0, c) = case.params
    # Eq numbers refer to DCMIP-TestCaseDocument_v1.6_23Jul2012
    # (81) => Pv(lat) = Pv_eq - c u_eq^2 sin^2 lat  (c u_eq ~ du/dPhi at the surface)
    # (82) => u(lat) = u_eq cos(lat) sqrt(Pv/Pv_eq + 2c*Phi * (Pv_eq/Pv) )
    # (84) => Phi(p,lat) = Pv*log(p_eq/p) - (Pv/2Pv_eq)*u_eq^2*sin^2 lat
    # sanity check : dPhi/dp = -Pv / p = -v
    slat, clat = sin(lat), cos(lat)
    pv   = pv0 - c*u0^2*slat^2
    Phi  = pv*log(p0/p) - (pv/2pv0)*u0^2*slat^2
    ulon, z = u0*clat*sqrt(pv/pv0 + 2c*Phi*pv0/pv ), zero(pv)
    return Phi, ulon, z, z
end

#============== Jablonowski (2006) baroclinic instability test cases ===============#

"""
    Jablonowski06
Type for Jablonowski (2006) and DCMIP 4.2 test cases (dry vs. moist). The only difference is specific humidity is included for the moist case.
"""
struct Jablonowski06{P} <: TestCase
    params::P
    function Jablonowski06(params)
        p=params[(:Uplanet, :gH0, :pV0, :delta_pV, :ps, :u0, :eta0, :etat, :etas,
                :up, :lonc, :latc, :q0, :pw, :latw)]
        new{typeof(p)}(p)
    end
end

const J06 = Jablonowski06
function default_params(::Type{J06})
    params = ( ps=1e5, u0=35.0, eta0=0.252, etat=0.2, etas=1.0, up=1.0, lonc=pi/9, latc=2pi/9, q0=0.021, pw=340*1e2, latw=2*pi/9 )
    R0, Omega, g, T0, DeltaT, Gamma, Rd = 6.4e6, 2pi/86400., 9.80616, 288.0, 4.8e5, 0.005, 287.0
    return ( Uplanet=R0*Omega, gH0=T0*g/Gamma, pV0=Rd*T0, delta_pV=Rd * DeltaT, params... )
end

function describe(case::J06)
    println("Jablonowski06 test case for baroclinic instability with parameters $(case.params).")
end

function initial_surface(lon, lat, case::J06)
    ps = case.params.ps
    Phis, _, _, _ = initial_flow(lon, lat, ps, case)
    return ps, Phis
end

function initial_flow(lon, lat, p, case::J06)
    (; Uplanet, gH0, pV0, delta_pV, ps, u0, eta0, etat, etas, latc, lonc, up, q0, pw, latw ) = case.params
    F=typeof(p)

    # pressure coordinate
    eta = p/ps
    etav = (eta - eta0) * F(pi)/2

    gH_prime = u0 * cos(etav)^(3//2) * ( (-2 * sin(lat)^6 * (cos(lat)^2 + 1//3) + 10//63) * u0 * cos(etav)^(3//2) + (8//5 * cos(lat)^3 * (sin(lat)^2 + 2//3) - F(pi)/4) * Uplanet )

    if etat <= eta <= etas
        # original expression :
        # gH_mean = T0 * g / Gamma * (1 - eta^(Rd * Gamma / g))

        # Rd * Gamma / g = Rd.T0 / T0*g/Gamma
        # Thus the basic profile has only two parameters,
        # Rd.T0=p0/rh0 and T0*g/Gamma = Phi0 (both homogeneous to a gepotential).
        gH_mean = gH0 * (1 - eta^(pV0/gH0))
    elseif eta < etat
        # gH_mean = T0 * g / Gamma * (1 - eta^(Rd * Gamma / g)) - Rd * DeltaT * ( (log(eta/etat) + 137/60) * etat^5 - 5*etat^4*eta + 5*etat^3*eta^2 - 10/3*etat^2*eta^3 + 5/4*etat*eta^4 - 1/5*eta^5)
        # Here we have another parameter Rd * DeltaT = Delta(p/rho)
        gH_mean = gH0 * (1 - eta^(pV0/gH0)) - delta_pV * ( (log(eta/etat) + 137//60) * etat^5 - 5*etat^4*eta + 5*etat^3*eta^2 - 10//3*etat^2*eta^3 + 5//4*etat*eta^4 - 1//5*eta^5)
    end

    gH = gH_mean + gH_prime

    # humidity
    q = q0 * exp(-(lat/latw)^4) * exp(-((eta - 1) * ps / pw)^2)

    # wind (purely zonal, balanced + perturbation)
    r = acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
    u = u0 * cos(etav)^(3//2) * sin(2 * lat)^2 + up*exp(-(10r)^2)
    v = 0

    return gH, u, v, q
end

end # priv

using .priv

end # module

