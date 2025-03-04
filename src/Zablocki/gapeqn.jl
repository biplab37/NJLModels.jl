## Normal Phase Δ_MF = 0

function En(p, m)
    return sqrt(p^2 + m^2)
end

function integrand_m(T, μ, m, ω, param::Parameters)
    return p -> 4 * 6 * param.Gs * (p^2 / 2π^2) * m *
                (1 - numberF(T, -μ, En(p, m) + ω) - numberF(T, μ, En(p, m) - ω)) / En(p, m)
end

function integrand_μ(T, μ, m, ω, param::Parameters)
    return p -> 4 * 6 * param.Gv * (p^2 / 2π^2) * m *
                (numberF(T, μ, En(p, m) - ω) - numberF(T, -μ, En(p, m) + ω))
end

function gapeqns(T, μ, param::Parameters)
    function gap(m, ω)
        return [
            m - param.m0 - integrate(integrand_m(T, μ, m, ω, param), 0, param.Λ),
            ω + integrate(integrand_μ(T, μ, m, ω, param), 0, param.Λ),
        ]
    end
    return x -> gap(x[1], x[2])
end

function gapeqn_m(T, mu, param::Parameters)
    return m -> m - param.m0 - integrate(integrand_m(T, mu, m, 0, param), 0, param.Λ)
end

"""
    massgap(T, μ, param::Parameters)
    massgap(trange::AbstractRange, μ, param::Parameters; initial_guess = [0.4, 0.0])
    massgap(T, μrange::AbstractRange, param::Parameters; initial_guess = [0.4, 0.0])
    massgap(trange::AbstractRange, μrange::AbstractRange, param::Parameters; initial_guess = [0.4, 0.0])

Solves the gap equation in normal phase for a given temperature and chemical potential. 
Returns the quark mass and the ω condensate.
"""
function massgap(T, μ, param::Parameters)
    return nlsolve(gapeqns(T, μ, param), [0.3, 0.01])
end

function massgap(trange::AbstractRange, μ, param::Parameters; initial_guess=[0.01, 0.01])
    result_m = zeros(length(trange))
    result_ω = zeros(length(trange))
    sol = initial_guess
    for (i, t) in enumerate(trange)
        sol = nlsolve(gapeqns(t, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

function massgap(T, μrange::AbstractRange, param::Parameters; initial_guess=[0.2, 0.0])
    result_m = zeros(length(μrange))
    result_ω = zeros(length(μrange))
    sol = initial_guess
    for (i, μ) in enumerate(μrange)
        sol = nlsolve(gapeqns(T, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

function massgap(trange::AbstractRange, μrange::AbstractRange, param::Parameters; initial_guess=[0.4, 0.0])
    if length(trange) != length(μrange)
        throw(DimensionMismatch("temperature list and the chemical potential list must have the same length"))
    end
    result_m = zeros(length(trange))
    result_ω = zeros(length(trange))
    sol = initial_guess
    for (i, (t, μ)) in enumerate(zip(trange, μrange))
        sol = nlsolve(gapeqns(t, μ, param), sol).zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
    end
    return result_m, result_ω
end

function massgap_m(T, mu, param::Parameters)
    return fzero(gapeqn_m(T, mu, param), 0.3)
end

export massgap
