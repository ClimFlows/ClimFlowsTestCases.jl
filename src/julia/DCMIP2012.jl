#================= DCMIP4.2 : (optionnaly sheared) solid-body rotation ===============#

"""
    struct DCMIP{N} <: TestCaseHPE

Type for DCMIP 2012 test cases. Currently only N=42 (test case 4.2) is implemented.
"""
struct DCMIP{N,P} <: TestCaseHPE
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
    "DCMIP 4.2 test case with reference values p=$p0, u=$u0, rho=$(p0/pv0), du/dPhi=$(c*u0)"
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

