
struct DCMIP{N,P} <: TestCaseHPE
    params :: P

    function DCMIP{42}(F=Float64; user...) 
        T0, Rd = 300.0, 287.0
        p = override(F, (p0=1e5, u0=20.0, pv0=Rd*T0, c=0.0), user)
        new{42, typeof(p)}(p)
    end

    function DCMIP{20}(F=Float64; user...) 
        # parameters
        # λₘ        lambda_m
        # ϕₘ        phi_m
        # Rₘ        radius_m
        # ζₘ        zeta_m
        # h₀        Phi_m = g h₀
        # p₀        p0
        # T₀        pv0 = RT₀
        # Γ         dpv_dPhi = R Γ / g
        T0, Rd , Gamma, g = 300.0, 287.0, 0.0065, 9.81 # non-generic parameters, from DCMIP 2012 document v1.6_23 p.26
        mountain = (lambda_m=3pi/2, phi_m=0.0, radius_m=3pi/4, zeta_m=pi/16, Phi_m=2000*g)
        p = override(F, (p0=1e5, pv0=Rd*T0, dpv_dPhi=Rd*Gamma/g, mountain...), user)
        new{20, typeof(p)}(p)
    end

    function DCMIP{21}(F=Float64; user...) 
        # parameters
        # λ_c        lambda_m
        # ϕ_c        phi_m
        # d          delta_m
        # ξ          xi_m
        # h_0        Phi_m = g h_0
        # p_eq       p0
        # T_eq       pv0 = RT₀
        # u_eq       u0

        # non-generic parameters
        X = 500
        Rd, radius, g = 287.0, 6.37122e6/X, 9.81 # from DCMIP 2012 document v1.6_23 p.4
        T0, d, xi, h0 = 300.0, 5000.0, 400.0, 250.0 # from DCMIP 2012 document v1.6_23 p.30
        mountain = (lambda_m=pi/4, phi_m=0.0, delta_m=d/radius, xi_m=xi/radius, Phi_m=h0*g)
        p = override(F, (p0=1e5, pv0=Rd*T0, u0=20.0, mountain...), user)
        new{21, typeof(p)}(p)
    end

end

#================= DCMIP2.0 : State of rest over orography ===============#

const DCMIP20 = DCMIP{20}

function describe(case::DCMIP20)
    "DCMIP 2.0 test case: state of rest over orography"
end

#================= DCMIP2.1 : Mountain waves ===============#

const DCMIP21 = DCMIP{21}

function describe(case::DCMIP21)
    "DCMIP 2.1 test case: orographic wave"
end

function initial(case::DCMIP21, lon, lat)
    (; lambda_m, phi_m, delta_m, xi_m, Phi_m, p0, pv0, u0) = case.params
    r = dist((lon, lat), (lambda_m, phi_m)) 
    Phis = Phi_m * exp(-r^2/delta_m^2)*cos(pi*r/xi_m)^2
    ps = p0*exp(-(Phis + sin(lat)^2*u0^2/2)/pv0)
    return ps, Phis
end

function initial(case::DCMIP21, lon, lat, p)
    (; p0, pv0, u0) = case.params
    Phi = pv0*log(p0/p) - sin(lat)*u0^2/2    
    ulon = u0*cos(lat)
    z = zero(Phi)
    return Phi, ulon, z, z
end

#================= DCMIP4.2 : Baroclinic instability ===============#

const DCMIP42 = DCMIP{42}

function describe(case::DCMIP42)
    (; p0, u0, pv0, c) = case.params
    "DCMIP 4.2 test case: baroclinic instability with reference values p=$p0, u=$u0, rho=$(p0/pv0), du/dPhi=$(c*u0)"
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
