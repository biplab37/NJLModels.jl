# This file contains code for the analytic structure of mesons

@doc raw"""
    Π0_meson_analytic(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)

This function stiches two Riemann scheet together. The physical sheet on the upeer half plane and the second sheet on the lower half.
"""
function Π0_meson_analytic(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)
    pauli(p) = 1 - numberF(T, mu, En(p, m)) - numberF(T, -mu, En(p, m))
    ep(p) = En(p, m)
    repart = 1 / (2 * param.Gs) - 2*Nc * Nf * (integrate(p -> (1 - m^2 / ep(p)^2)^(NM - 0.5) * p^2 * pauli(p) * (1 / (ep(p) + ω / 2) + 1 / (ep(p) - ω / 2)), 0, param.Λ)) / (4 * π^2)
    if imag(ω) > 0
        return repart
    elseif imag(ω) < 0
        return repart - 2im * imagpart_meson_q0_C(T, mu, ω, m, NM, Nc=Nc, Nf=Nf)
    else
        return realpart_meson_q0(T, mu, ω, m, NM, param, Nc=Nc, Nf=Nf)
    end
end

@doc raw"""
    find_meson_mass(T, mu, NM, param::Parameters)

This function calculates the mass of mesons. Note that NM = 0.5 for pion and NM=1.5 for sigma mesons.
Returns both the mass and the widths in a vector. When there is a bound state it just returns 0.0 for
the width.
"""
function find_meson_mass(T, mu, NM, param::Parameters)
    m, ome = massgap(T, mu, param).zero
    mu += ome
    repart(ω) = Π0_meson_analytic(T, mu, ω, m, NM, param)

    if repart(0.0) * repart(2 * m) < 0 # there is a bound state with 0<M<2m
        return [bisection(repart, 0.0, 2 * m), 0.0]
    end

    function ff!(F, x)
        term = repart(x[1] - 0.5im * x[2]) # Look for pole in the second Reimann sheet.
        F[1] = real(term)
        F[2] = imag(term)
    end
    return nlsolve(ff!, [2 * m, 0.1]).zero
end

find_pion_mass(T, mu, param) = find_meson_mass(T, mu, 0.5, param)
find_sigma_mass(T, mu, param) = find_meson_mass(T, mu, 1.5, param)

export find_pion_mass, find_sigma_mass