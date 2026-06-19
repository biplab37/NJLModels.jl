# rdf+pnjl

function Π0_meson_rdf(T, mu, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    pauli(p) = -f_polyakov(T, mu, En(p, m), Phi, Phibar) - fbar_polyakov(T, mu, En(p, m), Phi, Phibar)
    vacuum_term = 1 / (2 * param.Gs) - Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 / En(p, m), 0, param.Λ)) / (1 * π^2)
    medium_term = -Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 * pauli(p) / En(p, m), 0, 2.0)) / (1 * π^2)
    return vacuum_term + medium_term
end

function imagpart_meson_q0_C_rdf(T, mu, ω, m, NM, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    # NM=1/2 for pion and 3/2 for sigma
    factor = Nc * Nf / 8π
    pauli_term = 1 -f_polyakov(T, mu, ω / 2, Phi, Phibar) - fbar_polyakov(T, mu, ω / 2, Phi, Phibar)
    return factor * ω^2 * ((1 - (4 * m^2 / ω^2))^NM) * pauli_term
end

function imagpart_meson_q0_rdf(T, mu, ω, m, NM, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    if ω^2 < 4 * m^2
        return 0.0
    end
    return imagpart_meson_q0_C_rdf(T, mu, ω, m, NM, Phi, Phibar, Nc=Nc, Nf=Nf)
end

function imagpart_pi_q0_rdf(T, mu, ω, m, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    return imagpart_meson_q0_rdf(T, mu, ω, m, 0.5, Phi, Phibar, Nc=Nc, Nf=Nf)
end

function imagpart_sigma_q0_rdf(T, mu, ω, m, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    return imagpart_meson_q0_rdf(T, mu, ω, m, 1.5, Phi, Phibar, Nc=Nc, Nf=Nf)
end

function realpart_meson_q0_rdf(T, mu, ω, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    impart(x) = imagpart_meson_q0_rdf(T, mu, x, m, NM, Phi, Phibar, Nc=Nc, Nf=Nf)
    impart_vac(x) = imagpart_meson_q0_rdf(0.0, 0.0, x, m, NM, 0.0, 0.0, Nc=Nc, Nf=Nf)
    cutoff = 2 * sqrt(param.Λ^2 + m^2)
    vac_part = realpart_kramers_kronig(impart_vac, ω, cutoff)
    med_part = realpart_kramers_kronig(impart, ω, 2.0)
    return Π0_meson_rdf(T, mu, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf) - vac_part - med_part
end

function phase_shift_meson_q0_rdf(T, mu, ω, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    m = massgap_m(T, mu, param)
    impi = imagpart_meson_q0_rdf(T, mu, ω, m, NM, Phi, Phibar, Nc=Nc, Nf=Nf)
    repi = realpart_meson_q0_rdf(T, mu, ω, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    return atan(impi, repi)
end

phase_shift_pi_q0_rdf(T, mu, ω, param, Phi=0.0, Phibar=0.0; Nc=3, Nf=2) = phase_shift_meson_q0_rdf(T, mu, ω, 0.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)
phase_shift_sigma_q0_rdf(T, mu, ω, param, Phi=0.0, Phibar=0.0; Nc=3, Nf=2) = phase_shift_meson_q0_rdf(T, mu, ω, 1.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)

@doc raw"""
    Π0_meson_analytic_rdf(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)

This function stiches two Riemann scheet together. The physical sheet on the upeer half plane and the second sheet on the lower half.
"""
function Π0_meson_analytic_rdf(T, mu, ω, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    if imag(ω) == 0
        return realpart_meson_q0_rdf(T, mu, ω, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    end
    pauli(p) = -f_polyakov(T, mu, En(p, m), Phi, Phibar) - fbar_polyakov(T, mu, En(p, m), Phi, Phibar)
    ep(p) = En(p, m)
    repart_vacuum = 1 / (2 * param.Gs) - 2 * Nc * Nf * (integrate(p -> (1 - m^2 / ep(p)^2)^(NM - 0.5) * p^2 * (1 / (ep(p) + ω / 2) + 1 / (ep(p) - ω / 2)), 0, param.Λ)) / (4 * π^2)
    repart_medium = -2 * Nc * Nf * (integrate(p -> (1 - m^2 / ep(p)^2)^(NM - 0.5) * p^2 * pauli(p) * (1 / (ep(p) + ω / 2) + 1 / (ep(p) - ω / 2)), 0, 2.0)) / (4 * π^2)
    if imag(ω) > 0
        return repart_vacuum + repart_medium
    elseif imag(ω) < 0
        return repart_vacuum + repart_medium - 2im * sum(imagpart_meson_q0_rdf(T, mu, ω, m, NM, Phi, Phibar, Nc=Nc, Nf=Nf))
    end
end

@doc raw"""
    find_meson_mass(T, mu, NM, param::Parameters)

This function calculates the mass of mesons. Note that NM = 0.5 for pion and NM=1.5 for sigma mesons.
Returns both the mass and the widths in a vector. When there is a bound state it just returns 0.0 for
the width.
"""
function find_meson_mass_rdf(T, mu, m, NM, param::Parameters, Phi=0.0, Phibar=0.0)
    repart(ω) = Π0_meson_analytic_rdf(T, mu, ω, m, NM, param)

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

find_pion_mass_rdf(T, mu, m, param, Phi=0.0, Phibar=0.0) = find_meson_mass_rdf(T, mu, m, 0.5, param, Phi, Phibar)
find_sigma_mass_rdf(T, mu, m, param, Phi=0.0, Phibar=0.0) = find_meson_mass_rdf(T, mu, m, 1.5, param, Phi, Phibar)

export find_pion_mass_rdf, find_sigma_mass_rdf


## Finite temperature and chemical potential 
#TODO: