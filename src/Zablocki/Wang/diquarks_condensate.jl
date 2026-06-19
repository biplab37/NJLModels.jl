## This file contains code for diquarks in the superconducting phase, where Δ ≠ 0.
# the definitions here follows the Wang, Wang, Rischke paper.

function Ekk(k, m)
    return sqrt(k^2 + m^2)
end

function xi_k(k, m, mu, sign)
    return sign * Ekk(k, m) - mu
end

function eps_k(k, m, mu, Δ, sign)
    return sqrt(xi_k(k, m, mu, sign)^2 + Δ^2)
end

function ckpk(k, p, m, e, ep, u)
    # u = cos(theta)
    return 1 + e * ep * ((k * p * u + k^2) + m^2) / (Ekk(k, m) * Ekk(sqrt(k^2 + p^2 + 2k * p * u), m))
end

function integrand_Pi_0_1(T, mu, ω, k, p, u, m, Δ, part_function::Function)
    p_plus_k = sqrt(p^2 + k^2 + 2 * p * k * u)
    eps1_sign(sign) = eps_k(k, m, mu, Δ, sign)
    eps2_sign(sign) = eps_k(p_plus_k, m, mu, Δ, sign)
    xi1_sign(sign) = xi_k(k, m, mu, sign)
    xi2_sign(sign) = xi_k(p_plus_k, m, mu, sign)
    ff(x) = FD_dist(T, 0.0, x)
    prefactor(e, ep) = k^2 * ckpk(k, p, m, e, ep, u) / (2 * π^2)
    term1(e, ep, e1, ep1) = (ep1 * eps1_sign(ep) + xi1_sign(ep)) * (1 - ff(ep1 * eps1_sign(ep)) - ff(xi2_sign(e))) * part_function(ω - ep1 * eps1_sign(ep) - xi2_sign(e)) / (2 * ep1 * eps1_sign(ep))
    term2(e, ep, e1, ep1) = (e1 * eps2_sign(e) + xi2_sign(e)) * (1 - ff(e1 * eps2_sign(e)) - ff(xi2_sign(ep))) * part_function(ω - e1 * eps2_sign(e) - xi2_sign(ep)) / (2 * e1 * eps2_sign(e))

    signs = [+1, -1]
    return 0.5 * sum(prefactor(e, ep) * (term1(e, ep, e1, ep1) + term2(e, ep, e1, ep1)) for e = signs, ep = signs, e1 = signs, ep1 = signs)
end

function integrand_Pi_0_3(T, mu, ω, k, p, u, m, Δ, part_function::Function)
    p_plus_k = sqrt(p^2 + k^2 + 2 * p * k * u)
    eps1_sign(sign) = eps_k(k, m, mu, Δ, sign)
    eps2_sign(sign) = eps_k(p_plus_k, m, mu, Δ, sign)
    xi1_sign(sign) = xi_k(k, m, mu, sign)
    xi2_sign(sign) = xi_k(p_plus_k, m, mu, sign)
    ff(x) = FD_dist(T, 0.0, x)
    prefactor(e, ep) = k^2 * ckpk(k, p, m, e, ep, u) / (π^2)
    term(e, ep, e1, ep1) = (ep1 * eps1_sign(ep) + xi1_sign(ep)) * (e1 * eps2_sign(e) + xi2_sign(e)) * (1 - ff(ep1 * eps1_sign(ep)) - ff(e1 * eps2_sign(e))) * part_function(ω - ep1 * eps1_sign(ep) - e1 * eps2_sign(e)) / (4 * ep1 * e1 * eps1_sign(ep) * eps2_sign(e))

    signs = [+1, -1]
    return 0.5 * sum(prefactor(e, ep) * (term(e, ep, e1, ep1)) for e = signs, ep = signs, e1 = signs, ep1 = signs)
end

function integrand_Pi_1_3(T, mu, ω, k, p, u, m, Δ, part_function::Function)
    p_plus_k = sqrt(p^2 + k^2 + 2 * p * k * u)
    eps1_sign(sign) = eps_k(k, m, mu, Δ, sign)
    eps2_sign(sign) = eps_k(p_plus_k, m, mu, Δ, sign)
    xi1_sign(sign) = xi_k(k, m, mu, sign)
    xi2_sign(sign) = xi_k(p_plus_k, m, mu, sign)
    ff(x) = FD_dist(T, 0.0, x)
    prefactor(e, ep) = k^2 * ckpk(k, p, m, e, ep, u) / (4 * π^2)
    term(e, ep, e1, ep1) = Δ^2 * (1 - ff(e1 * eps1_sign(e)) - ff(ep1 * eps2_sign(ep))) * part_function(ω - ep1 * eps2_sign(ep) - e1 * eps1_sign(e)) / (ep1 * e1 * eps1_sign(e) * eps2_sign(ep))

    signs = [+1, -1]
    return 0.5 * sum(prefactor(e, ep) * (term(e, ep, e1, ep1)) for e = signs, ep = signs, e1 = signs, ep1 = signs)
end

function real_part_function(x; η=1e-6)
    return x / (x^2 + η^2)
end

function imag_part_function(x; η=1e-4)
    return -η / (x^2 + η^2)
end

function Pi_1(T, mu, ω, p, m, Δ, part_function, param)
    return integrate(x -> integrand_Pi_0_1(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0], rtol=1e-4, maxevals=100000)
end
function Pi_R_3(T, mu, ω, p, m, Δ, part_function, param)
    return 0.5 * (integrate(x -> integrand_Pi_0_3(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0]) +
                  integrate(x -> integrand_Pi_1_3(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0]))
end
function Pi_I_3(T, mu, ω, p, m, Δ, part_function, param)
    return 0.5 * (integrate(x -> integrand_Pi_0_3(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0]) -
                  integrate(x -> integrand_Pi_1_3(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0]))
end

function Pi_0_3(T, mu, ω, p, m, Δ, part_function, param)
    return integrate(x -> integrand_Pi_0_3(T, mu, ω, x[1], p, x[2], m, Δ, part_function), [0.0, -1.0], [param.Λ, 1.0])
end

function spectral_function_1(T, mu, ω, p, param)
    m, _, Δ = massgap_full(T, mu, param).zero
    if ω < 2m
        return 0.0
    end
    # if abs(Δ) < 5e-2
    #     return spectral_function_without_condensate(T, mu, ω, p, param)
    # end

    repart = -(1 / (4 * param.GD)) - Pi_1(T, mu, ω, p, m, Δ, real_part_function, param)
    impart = -Pi_1(T, mu, ω, p, m, Δ, imag_part_function, param)

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_R_3(T, mu, ω, p, param)
    m, _, Δ = massgap_full(T, mu, param).zero

    repart = -(1 / (4 * param.GD)) - Pi_R_3(T, mu, ω, p, m, Δ, real_part_function, param)
    impart = -Pi_R_3(T, mu, ω, p, m, Δ, imag_part_function, param)

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_I_3(T, mu, ω, p, param)
    m, _, Δ = massgap_full(T, mu, param).zero

    repart = -(1 / (4 * param.GD)) - Pi_I_3(T, mu, ω, p, m, Δ, real_part_function, param)
    impart = -Pi_I_3(T, mu, ω, p, m, Δ, imag_part_function, param)

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_without_condensate(T, mu, ω, p, param)
    m, ome = massgap(T, mu, param).zero
    mu += ome
    repart = realpart_D_normal(T, mu, ω, p, m, param)
    impart = imagpart_D_normal(T, mu, ω, p, m, param)

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_without_condensate(T, mu, ω, p, m, param)
    repart = realpart_D_normal(T, mu, ω, p, m, param)
    impart = imagpart_D_normal(T, mu, ω, p, m, param)

    return impart / (π * (repart^2 + impart^2))
end

function diquark_propagator_without_condensate(T, mu, p0, p, m, param)
    spectral(ω) = spectral_function_without_condensate(T, mu, ω, p, m, param)
    cutoff = sqrt(param.Λ^2 + m^2 + p^2 / 4)
    return UsefulFunctions.PVintegral(spectral, -2 * (cutoff + mu), 2 * (cutoff - mu), p0, integrate)
end

function diquark_propagator_without_condensate(T, mu, p0, p, param)
    m, ome = massgap(T, mu, param).zero
    mu += ome
    return diquark_propagator_without_condensate(T, mu, p0, p, m, param)
end

function diquark_propagator_without_condensate1(T, mu, p0, p, param)
    m, _, Δ = massgap_full(T, mu, param).zero

    repart = -(1 / (4 * param.GD)) - Pi_I_3(T, mu, p0, p, m, Δ, real_part_function, param)
    impart = -Pi_I_3(T, mu, p0, p, m, Δ, imag_part_function, param)

    return -repart / ((repart^2 + impart^2))
end
