#============== Isothermal test case ===============#

struct Isothermal{P} <: TestCaseHPE
    params::P
    function Isothermal(F=Float64; user...)
        defaults = (psurf=1e5, vsurf=1.0, lon0=3pi/2, lat0=0, Phi0=1e4, zeta0=pi/16, R0=3pi/4)
        params = override(F, defaults, user)
        new{typeof(params)}(params)
    end
end

function describe(case::Isothermal)
    params = join(("    $first=$second" for (;first, second) in pairs(case.params)), "\n")
    "Isothermal state of rest over mountain, with parameters: \n$params"
end

function initial(case::Isothermal, lon, lat, p)
    # p = ps * exp(-Phi/(psurf*vsurf) => dp = -p/vsurf
    (; psurf, vsurf) = case.params
    z = zero(psurf)
    return -psurf*vsurf*log(p/psurf), z, z, z
end

function initial(case::Isothermal, lon, lat)
    (; psurf, vsurf, lon0, lat0, Phi0, zeta0, R0) = case.params
    # borrowed from DCMIP 2012 test case document, v1.6, Eq (63) page 27
    θ = min(R0, acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0)))
    Phi = Phi0 * (1+cos(pi*θ/R0))*cos(pi*θ/zeta0)^2  
    return psurf*exp(-Phi/(psurf*vsurf)), Phi 
end
