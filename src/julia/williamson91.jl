"""
    Williamson91{N}
Type for Williamson (1991) test cases. Currently only N=6 is implemented.
"""
struct Williamson91{N,P} <: TestCaseSW
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
    """
    Williamson 1991 test case number 6 (Rossby-Haurwitz wave).
    Angular velocity of Rossby-Haurwitz wave : $nu
    Time to advance by 1/2 period : $(pi/nu/n/3600)h
    """
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

describe(::W91_2) = "Williamson 1991 test case number 2 (Global steady state nonlinear geostrophic flow)."

function initial_flow(lon, lat, case::W91_2)
    (; R0, Omega, gH0, u0, alpha) = case.params
    u = u0*( cos(lat)*cos(alpha) + cos(lon)*sin(lat)*sin(alpha) )
    v = -u0*sin(lon)*sin(alpha)
    gH = gH0 - ( R0*Omega*u0 + u0^2/2 )*( sin(lat)*cos(alpha) - cos(lon)*cos(lat)*sin(alpha) )^2
    return u, v, gH
end
