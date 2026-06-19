## handles the integral with diquark condensate better.


# with nested integral
function _im_Pi_0_1(T, mu, ω, p, m, Δ, param)
    integrand(k, u) = integrand_Pi_0_1(T, mu, ω, k, p, u, m, Δ, imag_part_function)
    u_integrand(k) = integrate(u -> integrand(k, u), -1.0, 1.0)
    return integrate(u_integrand, 0.0, param.Λ)
end