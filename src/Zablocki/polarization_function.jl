function Π0_meson(T, mu, m, NM, param::Parameters; Nc=3, Nf=2)
    pauli(p) = 1 - numberF(T, mu, En(p, m)) - numberF(T, -mu, En(p, m))
    return 1 / (2 * param.Gs) - Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 * pauli(p) / En(p, m), 0, param.Λ)) / (2 * π^2)
end

function imagpart_meson_q0(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)
    # NM=1/2 for pion and 3/2 for sigma
    if ω^2 < 4 * m^2 || ω^2 > 4 * (param.Λ^2 + m^2)
        return 0.0
    end
    factor = Nc / 4π
    pauli_term = 1 - numberF(T, mu, ω / 2) - numberF(T, -mu, ω / 2)
    return factor * ω^2 * ((1 - (4 * m^2 / ω^2))^NM) * pauli_term
end

function imagpart_pi_q0(T, mu, ω, m, param::Parameters; Nc=3, Nf=2)
    return imagpart_meson_q0(T, mu, ω, m, 0.5, param, Nc=Nc, Nf=Nf)
end

function imagpart_sigma_q0(T, mu, ω, m, param::Parameters; Nc=3, Nf=2)
    return imagpart_meson_q0(T, mu, ω, m, 1.5, param, Nc=Nc, Nf=Nf)
end

function realpart_meson_q0(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)
    impart(x) = imagpart_meson_q0(T, mu, x, m, NM, param, Nc=Nc, Nf=Nf)
    cutoff = 2 * sqrt(param.Λ^2 + m^2)
    integrand(ν) = 2 * ν * impart(ν) * (PrincipalValue(ν^2 - ω^2) - PrincipalValue(ν^2)) / π
    return Π0_meson(T, mu, m, NM, param) - integrate(integrand, 0.0, cutoff)
end

function fullrealpart_meson_q0(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)
    impart(x) = imagpart_meson_q0(T, mu, x, m, NM, param, Nc=Nc, Nf=Nf)
    cutoff = 2 * sqrt(param.Λ^2 + m^2)
    integrand(ν) = 2 * ν * impart(ν) * (PrincipalValue(ν^2 - ω^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function phase_shift_meson_q0(T, mu, ω, NM, param::Parameters; Nc=3, Nf=2)
    m = massgap_m(T, mu, param)
    impi = imagpart_meson_q0(T, mu, ω, m, NM, param, Nc=Nc, Nf=Nf)
    repi = realpart_meson_q0(T, mu, ω, m, NM, param, Nc=Nc, Nf=Nf)
    return atan(impi, repi)
end

phase_shift_pi_q0(T, mu, ω, param; Nc=3, Nf=2) = phase_shift_meson_q0(T, mu, ω, 0.5, param, Nc=Nc, Nf=Nf)
phase_shift_sigma_q0(T, mu, ω, param; Nc=3, Nf=2) = phase_shift_meson_q0(T, mu, ω, 1.5, param, Nc=Nc, Nf=Nf)