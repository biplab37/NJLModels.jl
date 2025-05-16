# This file contains code to calculate the effect on fluctuation pressure due to momentum dependence of phase shift.
function generalized_beth(theta)
    return theta - (sin(2 * theta) / 2)
end

function integrand_pressure(phase_shift::Function, T, mu, ω, q, degeneracy, param)
    return degeneracy * phase_shift(T, mu, ω, q, param) * q^2 * numberB(T, 0.0, ω) / (2 * π^3)
end

function integrand_pressure_generalized(phase_shift::Function, T, mu, ω, q, degeneracy, param)
    generalized_func(T, mu, ω, q, param) = generalized_beth(phase_shift(T, mu, ω, q, param))
    return integrand_pressure(generalized_func, T, mu, ω, q, degeneracy, param)
end

function phase_shift_LD(phase_shift::Function, T, mu, ω, q, param)
    if ω > q
        return 0.0
    end
    return phase_shift(T, mu, ω, q, param)
end

function phase_zero(phase_shift::Function, T, mu, ω, q, param)
    if ω <= q
        return 0.0
    end
    return phase_shift(T, mu, sqrt(ω^2 - q^2), 0.001, param)
end

function pressure_fluctuation(integrand_func::Function, T, mu, param; rtol=1e-2, LD=1.0)
    integrand(x) = integrand_func(T, mu, x[1], x[2], param)

    return integrate(integrand, [0.01, 0.01], [5 * param.Λ, LD * param.Λ], rtol=rtol)
end

# Pseudo-scalar
integrand_pressure_pi(T, mu, ω, q, param) = integrand_pressure(phase_shift_pi_q, T, mu, ω, q, 3, param)
integrand_pressure_pi_LD(T, mu, ω, q, param) = (ω >= q) ? 0.0 : integrand_pressure(phase_shift_pi_q, T, mu, ω, q, 3, param)

integrand_pressure_pi_generalized(T, mu, ω, q, param) = integrand_pressure_generalized(phase_shift_pi_q, T, mu, ω, q, 3, param)
integrand_pressure_pi_generalized_LD(T, mu, ω, q, param) = (ω >= q) ? 0.0 : integrand_pressure_generalized(phase_shift_pi_q, T, mu, ω, q, 3, param)

pressure_fluctuation_pi(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_pi, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_pi_LD(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_pi_LD, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_pi_generalized(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_pi_generalized, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_pi_generalized_LD(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_pi_generalized_LD, T, mu, param, rtol=rtol, LD=LD)

# Scalar
phase_shift_sigma_LD(T, mu, ω, q, param) = phase_shift_LD(phase_shift_sigma_q, T, mu, ω, q, param)
integrand_pressure_sigma(T, mu, ω, q, param) = integrand_pressure(phase_shift_sigma_q, T, mu, ω, q, 1, param)
integrand_pressure_sigma_LD(T, mu, ω, q, param) = integrand_pressure(phase_shift_sigma_LD, T, mu, ω, q, 1, param)

integrand_pressure_sigma_generalized(T, mu, ω, q, param) = integrand_pressure_generalized(phase_shift_sigma_q, T, mu, ω, q, 1, param)
integrand_pressure_sigma_generalized_LD(T, mu, ω, q, param) = integrand_pressure_generalized(phase_shift_sigma_LD, T, mu, ω, q, 1, param)

pressure_fluctuation_sigma(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_sigma, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_sigma_LD(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_sigma_LD, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_sigma_generalized(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_sigma_generalized, T, mu, param, rtol=rtol, LD=LD)
pressure_fluctuation_sigma_generalized_LD(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_pressure_sigma_generalized_LD, T, mu, param, rtol=rtol, LD=LD)

## Boosted approximation
phase_shift_zero_pi(T, mu, ω, q, param) = phase_zero(phase_shift_pi_q, T, mu, ω, q, param)
integrand_boosted_pi(T, mu, ω, q, param) = integrand_pressure(phase_shift_zero_pi, T, mu, ω, q, 3, param)
pressure_boosted_pi(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_boosted_pi, T, mu, param, rtol=rtol, LD=LD)

phase_shift_zero_sigma(T, mu, ω, q, param) = phase_zero(phase_shift_sigma_q, T, mu, ω, q, param)
integrand_boosted_sigma(T, mu, ω, q, param) = integrand_pressure(phase_shift_zero_sigma, T, mu, ω, q, 1, param)
pressure_boosted_sigma(T, mu, param; rtol=1e-2, LD=1.0) = pressure_fluctuation(integrand_boosted_sigma, T, mu, param, rtol=rtol, LD=LD)
