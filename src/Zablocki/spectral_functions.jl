## This file contains code related to computing spectral functions

function spectral_function_meson_q(T, mu, ω, q, Nm, param::Parameters)
    gap = massgap(T, mu, param).zero

    return spectral_function_meson_q(T, mu, ω, q, gap, Nm, param)
end

function spectral_function_meson_q(T, mu, ω, q, gap, Nm, param::Parameters)
    impart = imagpart_meson_q(ω, T, mu, q, gap, Nm, param)
    repart = realpart_meson_q(ω, T, mu, q, gap, Nm, param)

    return impart / (π * (impart^2 + repart^2))
end

spectral_function_pi_q(T, mu, ω, q, gap, param::Parameters) = spectral_function_meson_q(T, mu, ω, q, gap, 0, param)
spectral_function_sigma_q(T, mu, ω, q, gap, param::Parameters) = spectral_function_meson_q(T, mu, ω, q, gap, 1, param)
