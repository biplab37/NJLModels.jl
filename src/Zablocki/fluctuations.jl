# This file contains code to calculate the effect on fluctuation pressure due to momentum dependence of phase shift.

function integrand_pressure(T, mu, ω, q, param)
    return 3 * phase_shift_pi_q(T, mu, ω, q, param) * q^2 * (1 / (exp((ω - mu) / T) - 1.0) + 1 / (exp((ω + mu) / T) - 1)) / (4 * π^3)
end

function pressure_fluctuation(T, mu, param; rtol=1e-2)
    integrand(x) = integrand_pressure(T, mu, x[1], x[2], param)
    return integrate(integrand, [0.0, 0.0], [sqrt(5) * param.Λ, param.Λ], rtol=rtol)
end

function pressure_fluctuation(trange::AbstractArray, mu, param)
    pres = zeros(length(trange))

    Threads.@threads for i in eachindex(trange)
        pres[i] = pressure_fluctuation(trange[i], mu, param)[1]
    end

    return pres
end

function phase_zero(T, mu, ω, q, param)
    if ω <= q
        return 0.0
    end
    return phase_shift_sigma_q(T, mu, sqrt(ω^2 - q^2), 0.001, param)
end

function pressure_fluctuation_boosted(T, mu, param)
    integrand(x) = 3 * phase_zero(T, mu, x[1], x[2], param) * x[2]^2 * (1 / ((exp((x[1] - mu) / T) - 1.0)) + 1 / ((exp((x[1] + mu) / T) - 1.0))) / (4 * π^3)
    return integrate(integrand, [0.0, 0.0], [sqrt(5) * param.Λ, param.Λ], rtol=1e-2)
end

function pressure_fluctuation_boosted(trange::AbstractArray, mu, param)
    pressure_boost = zeros(length(trange))
    Threads.@threads for i in eachindex(trange)
        pressure_boost[i] = pressure_fluctuation_boosted(trange[i], mu, param)[1]
    end
    return pressure_boost
end

export pressure_fluctuation, pressure_fluctuation_boosted
