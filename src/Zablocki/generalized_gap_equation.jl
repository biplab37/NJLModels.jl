### This file contains code to calculate the solution of the gap equation including the fluctuations.

function potential_mean_field(T, mu, m, param)

    scalar_part = (m - param.m0)^2 / (4 * param.Gs)
    factor = 2 * 2 * 3 / (2 * π^2)
    divergentpart(p, m) = p^2 * En(p, m)
    mediumpart(p, m) = T * p^2 * (log(1 + exp(-(En(p, m) - mu) / T)) + log(1 + exp(-(En(p, m) + mu) / T)))

    div_part = integrate(q -> divergentpart(q, m), 0.0, param.Λ)
    # div_part_zero = integrate(q -> divergentpart(q, m2), 0.0, param.Λ)
    medium_term = integrate(q -> mediumpart(q, m), 0.0, 10.0)

    return scalar_part - factor * (div_part + medium_term)
end

function _integrand_potential_m(phase_shift::Function, T, mu, ω, q, m, degeneracy, param)
    return -degeneracy * phase_shift(T, mu, ω, q, m, param) * q^2 * numberB(T, 0.0, ω) / (2 * π^3)
end

function _potential_fluctuation_m(integrand_func::Function, T, mu, m, param; rtol=1e-2, LD=1.0)
    return integrate(x->integrand_func(T, mu, x[1], x[2], m, param), [0.01, 0.01], [5 * param.Λ, LD * param.Λ], rtol=rtol)
end

_integrand_potential_pi_m(T, mu, om, q, m, param) = _integrand_potential_m(phase_shift_pi_q_m, T, mu, om, q, m, 3, param)
_integrand_potential_sigma_m(T, mu, om, q, m, param) = _integrand_potential_m(phase_shift_sigma_q_m, T, mu, om, q, m, 1, param)

_potential_fluctuation_pi_m(T, mu, m, param; rtol=1e-2) = _potential_fluctuation_m(_integrand_potential_pi_m, T, mu, m, param, rtol=rtol)
_potential_fluctuation_sigma_m(T, mu, m, param; rtol=1e-2) = _potential_fluctuation_m(_integrand_potential_sigma_m, T, mu, m, param, rtol=rtol)

function potential_fluctuations(T, mu, m, param; rtol=1e-2)
    return _potential_fluctuation_pi_m(T, mu, m, param; rtol=rtol) + _potential_fluctuation_sigma_m(T, mu, m, param; rtol=rtol)
end

function potential_total(T, mu, m, param; rtol=1e-2)
    return potential_mean_field(T, mu, m, param) + potential_fluctuations(T, mu, m, param; rtol=rtol)
end

function gap_mean_field(T, mu, param)
    der(m) = UsefulFunctions._derivative(x->potential_mean_field(T, mu, x, param), m)

    return UsefulFunctions.bisection(der, 0.0, param.Λ)
end

function gap_fluctuations(T, mu, param; rtol=1e-2)
    der(m) = UsefulFunctions._derivative(x->potential_fluctuations(T, mu, x, param; rtol=rtol), m)

    return UsefulFunctions.NewtonRaphson(der, 0.3)
end

export potential_mean_field, potential_fluctuations, potential_total
