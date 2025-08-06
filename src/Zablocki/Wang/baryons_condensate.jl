
function quark_propagator_without_condensate11(p0, k0, k, mu, m)
    return ((p0 - k0 - mu) + m) * PrincipalValue((p0 - k0 + mu)^2 - (k^2 + m^2), 1e-3)
end

function quark_propagator_without_condensate22(p0, k0, k, mu, m)
    return (-(p0 - k0 - mu) + m) / ((p0 - k0 + mu)^2 - (k^2 + m^2))
end

function baryon_self_energy_without_condensate(T, mu, p0, param)
    m, ome = massgap(T, mu, param).zero
    mu += ome
    integrand(k0, k) = k^2 * quark_propagator_without_condensate11(p0, k0, k, mu, m) * diquark_propagator_without_condensate(T, mu, k0, k, m, param)
    return integrate(x -> integrand(x[1], x[2]), [-2.0, 0], [0.0, param.Λ], rtol=1e-1)
end

function J00(T, mu, param) # This integral is very costly numerically, it can be proven to be almost vainshing.
    m = massgap(T, mu, param).zero[1]
    func(k0) = quark_propagator_without_condensate11(0.0, k0, 0.0, mu, m) * diquark_propagator_without_condensate(T, mu, k0, 0.0, m, param)
    return integrate(func, -2 * param.Λ, 2 * param.Λ)
end

function imag_part_baryon_k(T, mu, p0, k, m, param)
    term1 = (1 + (m - 2 * mu) / (sqrt(k^2 + m^2))) * diquark_propagator_without_condensate(T, mu, p0 + mu - sqrt(k^2 + m^2), k, m, param)
    term2 = (1 - (m - 2 * mu) / (sqrt(k^2 + m^2))) * diquark_propagator_without_condensate(T, mu, p0 + mu + sqrt(k^2 + m^2), k, m, param)

    return 0.25 * (term1 + term2)
end

function real_part_baryon_k(T, mu, p0, k, m, param)
    func(x, q) = imag_part_baryon_k(T, mu, x, q, m, param)

    return realpart_kramers_kronig_q_1(func, p0, k, -2, 2, rtol=1e-2, maxevals=200, err=true)
end

function matsubara_summed_propagator_D(T, mu, p0, ω, k, m, param)
    term(sign) = 0.5 * (1 + sign * (m - 2 * mu) / sqrt(k^2 + m^2)) * (numberB(T, 0.0, ω) - numberB(T, -mu, p0 - sign * sqrt(k^2 + m^2))) / (p0 + mu - sign * sqrt(k^2 + m^2) - ω)

    return -0.25 * (term(+1) + term(-1))
end

function baryon_self_energy_integrand(T, mu, p0, ω, k, m, param)
    return k^2 * (matsubara_summed_propagator_D(T, mu, p0, ω, k, m, param) * spectral_function_without_condensate(T, mu, ω, k, m, param)) / (2 * π^2)
end

function baryon_self_energy_without_condensate(T, mu, p0, m, param)
    func(x) = baryon_self_energy_integrand(T, mu, p0, x[1], x[2], m, param)

    return integrate(func, [-5.0, 5.0], [0.0, 2.5 * param.Λ])
end