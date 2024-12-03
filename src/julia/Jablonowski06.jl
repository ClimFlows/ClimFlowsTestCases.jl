#============== Jablonowski (2006) baroclinic instability test cases ===============#

struct Jablonowski06{P} <: TestCaseHPE
    params::P
    function Jablonowski06(F=Float64; user...)
        other = ( ps=1e5, u0=35.0, eta0=0.252, etat=0.2, etas=1.0, up=1.0, lonc=pi/9, latc=2pi/9, q0=0.021, pw=340*1e2, latw=2*pi/9 )
        R0, Omega, g, T0, DeltaT, Gamma, Rd = 6.4e6, 2pi/86400., 9.80616, 288.0, 4.8e5, 0.005, 287.0
        defaults = merge( (Uplanet=R0*Omega, gH0=T0*g/Gamma, pV0=Rd*T0, delta_pV=Rd * DeltaT), other)
        p = override(F, defaults, user)
        new{typeof(p)}(p)
    end
end

const J06 = Jablonowski06

function describe(case::J06)
    params = join(("    $first=$second" for (;first, second) in pairs(case.params)), "\n")
    "Jablonowski06 test case for baroclinic instability with parameters: \n$params"
end

function initial(case::J06, lon, lat)
    ps = case.params.ps
    Phis, _, _, _ = initial(case, lon, lat, ps)
    return ps, Phis
end

function initial(case::J06, lon, lat, p)
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
        # Rd.T0=p0/rh0 and T0*g/Gamma = Phi0 (both homogeneous to a geopotential).
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
