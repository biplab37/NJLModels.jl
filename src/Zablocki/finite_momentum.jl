## This file contains code related to the calculation of polarization function
## at finite external momenta

function imagpart_meson_q(ω, T, mu, q, gap, Nm, param::Parameters)
    if q == 0
        return imagpart_meson_q0(T, mu, ω, gap[1], Nm + 1 / 2, param)
    end
    m = gap[1]
    mus = gap[2] + mu # μ⁰ = μ + ω₀
    s = ω^2 - q^2

    if 0 <= s <= 4 * m^2 || s >= 4 * (param.Λ^2 + m^2)
        return 0.0
    end

    disc = q * sqrt(1 - (4 * m^2 / s))
    Ep_plus = 0.5 * (ω + disc)
    Ep_minus = 0.5 * (ω - disc)

    factor = 3 * 2 * (s - (4 * m^2 * Nm)) / (8 * π * q)

    if s > 4 * m^2
        if ω > 0
            return factor * JM_plus(T, mus, Ep_plus, Ep_minus)
        else
            return factor * JM_minus(T, mus, Ep_plus, Ep_minus)
        end
    else # s<0
        return factor * JM_Landau(T, mus, Ep_plus, Ep_minus)
    end
end

function JM(T, mu, Ep_plus, Ep_minus)
    NN(x) = numberF(T, mu, x)
    return T * log((NN(-Ep_minus) * NN(Ep_minus)) / (NN(-Ep_plus) * NN(Ep_plus)))
end

JM_plus(T, mu, Ep_plus, Ep_minus) = JM(T, mu, Ep_plus, Ep_minus)
JM_minus(T, mu, Ep_plus, Ep_minus) = JM(T, -mu, Ep_plus, Ep_minus)

function JM_Landau(T, mu, Ep_plus, Ep_minus)
    NM(x) = numberF(T, mu, x)
    NP(x) = numberF(T, -mu, x)
    return T * log(NP(Ep_minus) * NM(Ep_minus) / (NP(-Ep_plus) * NM(-Ep_plus)))
end

imagpart_pi_q(ω, T, mu, q, gap, param) = imagpart_meson_q(ω, T, mu, q, gap, 0, param)
imagpart_sigma_q(ω, T, mu, q, gap, param) = imagpart_meson_q(ω, T, mu, q, gap, 1, param)

function realpart_meson_q(ω, T, mu, q, gap, Nm, param)
    m = gap[1]
    impart(x, y) = imagpart_meson_q(x, T, mu, y, gap, Nm, param)
    cutoff = 2 * sqrt(param.Λ^2 + m^2 + q^2 / 4)
    function integrand(ν)
        return 2 * ν * (impart(ν, q) * PrincipalValue(ν^2 - ω^2) - impart(ν, 0.0) * PrincipalValue(ν^2)) / π
    end
    return Π0_meson(T, mu + gap[2], m, Nm + 1 / 2, param) - integrate(integrand, 0.0, cutoff)
end

function fullrealpart_meson_q(ω, T, mu, q, gap, Nm, param)
    m = gap[1]
    impart(x, y) = imagpart_meson_q(x, T, mu, y, gap, Nm, param)
    cutoff = 2 * sqrt(param.Λ^2 + m^2 + q^2 / 4)
    function integrand(ν)
        return 2 * ν * impart(ν, q) * PrincipalValue(ν^2 - ω^2) / π
    end
    return 1 / (2 * param.Gs) - integrate(integrand, 0.0, cutoff)
end

function phase_shift_meson_q(T, mu, ω, q, Nm, param::Parameters)
    gap = massgap(T, mu, param).zero
    impi = imagpart_meson_q(ω, T, mu, q, gap, Nm, param)
    repi = realpart_meson_q(ω, T, mu, q, gap, Nm, param)
    return atan(impi, repi)
end

phase_shift_pi_q(T, mu, ω, q, param) = phase_shift_meson_q(T, mu, ω, q, 0, param)
phase_shift_sigma_q(T, mu, ω, q, param) = phase_shift_meson_q(T, mu, ω, q, 1, param)