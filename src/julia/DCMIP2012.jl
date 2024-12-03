#================= DCMIP4.2 : (optionnaly sheared) solid-body rotation ===============#

struct DCMIP{N,P} <: TestCaseHPE
    params :: P
    function DCMIP{42}(F=Float64; user...) 
        T0, Rd = 300.0, 287.0
        p = override(F, (p0=1e5, u0=20.0, pv0=Rd*T0, c=0.0), user)
        new{42, typeof(p)}(p)
    end
end

const DCMIP42 = DCMIP{42}

function describe(case::DCMIP42)
    (; p0, u0, pv0, c) = case.params
    "DCMIP 4.2 test case with reference values p=$p0, u=$u0, rho=$(p0/pv0), du/dPhi=$(c*u0)"
end

function initial(case::DCMIP42, lon, lat)
    p0 = case.params.p0
    Phi, _, _, _ = initial(case, lon, lat, p0)
    return p0, Phi
end

function initial(case::DCMIP42, lon, lat, p)
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
