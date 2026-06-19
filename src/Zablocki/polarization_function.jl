# about the factor 1/4 from the dubinin thesis.

function Π0_meson(T, mu, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    pauli(p) = -FD_dist(T, mu, En(p, m), Phi, Phibar) - FD_dist(T, -mu, En(p, m), Phi, Phibar)
    vacuum_term = 1 / (2 * param.Gs) - Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 / En(p, m), 0, param.Λ)) / (1 * π^2)
    medium_term = -Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 * pauli(p) / En(p, m), 0, 2.0)) / (1 * π^2)
    return vacuum_term + medium_term
end

function imagpart_meson_q0(T, mu, ω::Real, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    # NM=1/2 for pion and 3/2 ofr sigma
    if ω^2 < 4 * m^2
        return 0.0
    end
    vacuum_term, medium_term = imagpart_meson_q0_C(T, mu, ω, m, NM, Phi, Phibar, Nc=Nc, Nf=Nf)
    if ω^2 > 4 * (param.Λ^2 + m^2)
        vacuum_term = 0.0
    end
    return vacuum_term + medium_term
end

function imagpart_meson_q0_C(T, mu, ω, m, NM, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    # NM=1/2 for pion and 3/2 for sigma
    factor = Nc * Nf / 8π
    pauli_term = -FD_dist(T, mu, ω / 2, Phi, Phibar) - FD_dist(T, -mu, ω / 2, Phi, Phibar)
    return factor * ω^2 * ((1 - (4 * m^2 / ω^2))^NM) * [1, pauli_term]
end

function imagpart_pi_q0(T, mu, ω, m, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    return imagpart_meson_q0(T, mu, ω, m, 0.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)
end

function imagpart_sigma_q0(T, mu, ω, m, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    return imagpart_meson_q0(T, mu, ω, m, 1.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)
end

function realpart_meson_q0(T, mu, ω, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    impart(x) = imagpart_meson_q0(T, mu, x, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    cutoff = 2 * sqrt(param.Λ^2 + m^2)
    integrand(ν) = 2 * ν * impart(ν) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    return Π0_meson(T, mu, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf) - integrate(integrand, 0.0, 2.0)
end

function fullrealpart_meson_q0(T, mu, ω, m, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    impart(x) = imagpart_meson_q0(T, mu, x, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    cutoff = 2 * sqrt(param.Λ^2 + m^2)
    integrand(ν) = 2 * ν * impart(ν) * (PrincipalValue(ν^2 - ω^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function phase_shift_meson_q0(T, mu, ω, NM, param::Parameters, Phi=0.0, Phibar=0.0; Nc=3, Nf=2)
    m = massgap_m(T, mu, param)
    impi = imagpart_meson_q0(T, mu, ω, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    repi = realpart_meson_q0(T, mu, ω, m, NM, param, Phi, Phibar, Nc=Nc, Nf=Nf)
    return atan(impi, repi)
end

phase_shift_pi_q0(T, mu, ω, param, Phi=0.0, Phibar=0.0; Nc=3, Nf=2) = phase_shift_meson_q0(T, mu, ω, 0.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)
phase_shift_sigma_q0(T, mu, ω, param, Phi=0.0, Phibar=0.0; Nc=3, Nf=2) = phase_shift_meson_q0(T, mu, ω, 1.5, param, Phi, Phibar, Nc=Nc, Nf=Nf)
