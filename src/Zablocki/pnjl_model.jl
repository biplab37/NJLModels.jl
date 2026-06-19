# Polyakov loop included NJL model

Base.@kwdef struct PolyakovParameters
    T0::Float64 = 0.27
    a0::Float64 = 6.75
    a1::Float64 = -1.95
    a2::Float64 = 2.625
    a3::Float64 = -7.44
    b3::Float64 = 0.75
    b4::Float64 = 7.5
end

function U_Polyakov_derivs(T, Φ, Φbar, pyp::PolyakovParameters)
    if T < 1e-5
        return 0.0, 0.0
    end
    # Scale b2 by T^4 to be numerically stable at small T
    # b2(T) = a0 + a1 * (T0/T) + a2 * (T0/T)^2 + a3 * (T0/T)^3
    b2_scaled = pyp.a0 * T^4 + pyp.a1 * pyp.T0 * T^3 + pyp.a2 * pyp.T0^2 * T^2 + pyp.a3 * pyp.T0^3 * T
    T4 = T^4

    dU_dΦ = -0.5 * b2_scaled * Φbar - 0.5 * pyp.b3 * T4 * Φ^2 + 0.5 * pyp.b4 * T4 * Φbar^2 * Φ
    dU_dΦbar = -0.5 * b2_scaled * Φ - 0.5 * pyp.b3 * T4 * Φbar^2 + 0.5 * pyp.b4 * T4 * Φ^2 * Φbar

    return dU_dΦ, dU_dΦbar
end

function f_polyakov(T, μ_star, E, Φ, Φbar)
    if T < 1e-6
        return E - μ_star > 0 ? 0.0 : 1.0
    end
    d = (E - μ_star) / T
    if real(d) >= 0
        x = exp(-d)
        return (Φ * x + 2 * Φbar * x^2 + x^3) / (1 + 3 * Φ * x + 3 * Φbar * x^2 + x^3)
    else
        y = exp(d)
        return (Φ * y^2 + 2 * Φbar * y + 1) / (y^3 + 3 * Φ * y^2 + 3 * Φbar * y + 1)
    end
end

function fbar_polyakov(T, μ_star, E, Φ, Φbar)
    if T < 1e-6
        return E + μ_star > 0 ? 0.0 : 1.0
    end
    d = (E + μ_star) / T
    if real(d) >= 0
        x = exp(-d)
        return (Φbar * x + 2 * Φ * x^2 + x^3) / (1 + 3 * Φbar * x + 3 * Φ * x^2 + x^3)
    else
        y = exp(d)
        return (Φbar * y^2 + 2 * Φ * y + 1) / (y^3 + 3 * Φbar * y^2 + 3 * Φ * y + 1)
    end
end

function f_polyakov_C(T, μ_star, E, Φ, Φbar)
    if T < 1e-6
        return E - μ_star > 0 ? 0.0 : 1.0
    end
    d = (E - μ_star) / T
    if real(d) >= 0
        x = exp(-d)
        return (Φ * x + 2 * Φbar * x^2 + x^3) / (1 + 3 * Φ * x + 3 * Φbar * x^2 + x^3)
    else
        y = exp(d)
        return (Φ * y^2 + 2 * Φbar * y + 1) / (y^3 + 3 * Φ * y^2 + 3 * Φbar * y + 1)
    end
end


function I1_polyakov(T, μ_star, E, Φ, Φbar)
    if T < 1e-6
        return 0.0
    end

    # Term 1: exp(-beta*(E - μ_star)) / g+(E)
    d1 = (E - μ_star) / T
    if d1 >= 0
        x1 = exp(-d1)
        term1 = x1 / (1 + 3 * Φ * x1 + 3 * Φbar * x1^2 + x1^3)
    else
        y1 = exp(d1)
        term1 = y1^2 / (y1^3 + 3 * Φ * y1^2 + 3 * Φbar * y1 + 1)
    end

    # Term 2: exp(-2*beta*(E + μ_star)) / g-(E)
    d2 = (E + μ_star) / T
    if d2 >= 0
        x2 = exp(-d2)
        term2 = x2^2 / (1 + 3 * Φbar * x2 + 3 * Φ * x2^2 + x2^3)
    else
        y2 = exp(d2)
        term2 = y2 / (y2^3 + 3 * Φbar * y2^2 + 3 * Φ * y2 + 1)
    end

    return term1 + term2
end

function I2_polyakov(T, μ_star, E, Φ, Φbar)
    if T < 1e-6
        return 0.0
    end

    # Term 1: exp(-2*beta*(E - μ_star)) / g+(E)
    d1 = (E - μ_star) / T
    if d1 >= 0
        x1 = exp(-d1)
        term1 = x1^2 / (1 + 3 * Φ * x1 + 3 * Φbar * x1^2 + x1^3)
    else
        y1 = exp(d1)
        term1 = y1 / (y1^3 + 3 * Φ * y1^2 + 3 * Φbar * y1 + 1)
    end

    # Term 2: exp(-beta*(E + μ_star)) / g-(E)
    d2 = (E + μ_star) / T
    if d2 >= 0
        x2 = exp(-d2)
        term2 = x2 / (1 + 3 * Φbar * x2 + 3 * Φ * x2^2 + x2^3)
    else
        y2 = exp(d2)
        term2 = y2^2 / (y2^3 + 3 * Φbar * y2^2 + 3 * Φ * y2 + 1)
    end

    return term1 + term2
end

function _log_Z_phi_plus(T, mu, E, Φ, Φbar)
    if T < 1e-6
        return E - mu > 0 ? 0.0 : 3*(E + mu)
    end
    x = (E - mu)/T
    if x>=0
        return T*log(1 + 3 * Φbar * exp(-x) + 3 * Φ * exp(-2 * x) + exp(-3 * x))
    else
        return _log_Z_phi_minus(T, mu, -E, Φ, Φbar) + 3 * (E + mu) 
    end
end

function _log_Z_phi_minus(T, mu, E, Φ, Φbar)
    if T < 1e-6
        return E + mu > 0 ? 0.0 : 3*(E - mu)
    end
    x = (E + mu)/T
    if x>=0
        return T*log(1 + 3 * Φ * exp(-x) + 3 * Φbar * exp(-2 * x) + exp(-3 * x))
    else
        return _log_Z_phi_plus(T, mu, -E, Φ, Φbar) + 3 * (E - mu)
    end
end

function safe_En(p, m)
    abs_m = abs(m)
    if abs_m > 1e10
        return abs_m
    end
    return sqrt(p^2 + m^2)
end

function safe_m_over_E(p, m)
    abs_m = abs(m)
    if abs_m > 1e10
        return sign(m)
    end
    E = sqrt(p^2 + m^2)
    if E == 0.0
        return 0.0
    end
    return m / E
end

function integrand_m_pnjl(T, μ, m, ω, Φ, Φbar, param::Parameters)
    μ_star = μ + ω
    return p -> 12 * param.Gs * (p^2 / π^2) * safe_m_over_E(p, m) *
                (1 - f_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar) - fbar_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar))
end

function integrand_m_pnjl_pv(T, μ, m, ω, Φ, Φbar, param::Parameters)
    μ_star = μ + ω
    m1 = sqrt(m^2 + param.Λ^2)
    m2 = sqrt(m^2 + 2 * param.Λ^2)
    return p -> 12 * param.Gs * (p^2 / π^2) * (
                    (safe_m_over_E(p, m) - 2 * m / safe_En(p, m1) + m / safe_En(p, m2)) -
                    safe_m_over_E(p, m) * (f_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar) + fbar_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar))
                )
end

function integrand_μ_pnjl(T, μ, m, ω, Φ, Φbar, param::Parameters)
    μ_star = μ + ω
    return p -> 12 * param.Gv * (p^2 / π^2) *
                (f_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar) - fbar_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar))
end

function integrand_Φ_pnjl(T, μ, m, ω, Φ, Φbar)
    μ_star = μ + ω
    return p -> 6 * T * (p^2 / π^2) * I1_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar)
end

function integrand_Φbar_pnjl(T, μ, m, ω, Φ, Φbar)
    μ_star = μ + ω
    return p -> 6 * T * (p^2 / π^2) * I2_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar)
end

function gapeqns_pnjl(T, μ, param::Parameters, pyp::PolyakovParameters; regularization=:sharp, pv_scale=10.0)
    cutoff = (regularization == :pv) ? pv_scale * param.Λ : param.Λ

    function gap(m, ω, Φ, Φbar)
        if isnan(m) || isnan(ω) || isnan(Φ) || isnan(Φbar) || isinf(m) || isinf(ω) || isinf(Φ) || isinf(Φbar)
            return [1e4, 1e4, 1e4, 1e4]
        end
        m_c = clamp(m, -5.0, 5.0)
        ω_c = clamp(ω, -5.0, 5.0)
        Φ_c = clamp(Φ, 0.0, 1.0)
        Φbar_c = clamp(Φbar, 0.0, 1.0)
        dU_dΦ, dU_dΦbar = U_Polyakov_derivs(T, Φ_c, Φbar_c, pyp)

        int_m = (regularization == :pv) ?
                integrate(integrand_m_pnjl_pv(T, μ, m_c, ω_c, Φ_c, Φbar_c, param), 0, cutoff) :
                integrate(integrand_m_pnjl(T, μ, m_c, ω_c, Φ_c, Φbar_c, param), 0, cutoff)

        return [
            m - param.m0 - int_m,
            ω + integrate(integrand_μ_pnjl(T, μ, m_c, ω_c, Φ_c, Φbar_c, param), 0, cutoff),
            dU_dΦ - integrate(integrand_Φ_pnjl(T, μ, m_c, ω_c, Φ_c, Φbar_c), 0, cutoff),
            dU_dΦbar - integrate(integrand_Φbar_pnjl(T, μ, m_c, ω_c, Φ_c, Φbar_c), 0, cutoff)
        ]
    end
    return x -> gap(x[1], x[2], x[3], x[4])
end

"""
    massgap_pnjl(T, μ, param::Parameters, pyp::PolyakovParameters = PolyakovParameters(); initial_guess=[0.3, 0.01, 0.1, 0.1], regularization=:sharp, pv_scale=10.0)

Solves the PNJL gap equations for a given temperature and chemical potential.
Supports either `:sharp` or `:pv` (Pauli-Villars) regularization scheme.
Returns the nlsolve solution object. The variables are:
- `m`: Constituent quark mass
- `ω`: Vector condensate
- `Φ`: Polyakov loop
- `Φbar`: Conjugate Polyakov loop
"""
function massgap_pnjl(T, μ, param::Parameters, pyp::PolyakovParameters=PolyakovParameters(); initial_guess=[0.3, 0.01, 0.1, 0.1], regularization=:sharp, pv_scale=10.0)
    return mcpsolve(gapeqns_pnjl(T, μ, param, pyp; regularization=regularization, pv_scale=pv_scale), [0.0, 0.0, -1, -1], [0.5, 0.2, 0.5, 0.5], initial_guess, inplace=false)
end

"""
    massgap_pnjl(trange::AbstractVector, μ, param::Parameters, pyp::PolyakovParameters = PolyakovParameters(); initial_guess=[0.3, 0.01, 0.1, 0.1], regularization=:sharp, pv_scale=10.0)

Solves the PNJL gap equations sequentially for a vector of temperatures, feeding the solution of the previous temperature as the initial guess for the next.
Supports either `:sharp` or `:pv` (Pauli-Villars) regularization scheme.
Returns four vectors: `(result_m, result_ω, result_Φ, result_Φbar)`.
"""
function massgap_pnjl(trange::AbstractVector, μ, param::Parameters, pyp::PolyakovParameters=PolyakovParameters(); initial_guess=[0.3, 0.01, 0.1, 0.1], regularization=:sharp, pv_scale=10.0)
    result_m = zeros(length(trange))
    result_ω = zeros(length(trange))
    result_Φ = zeros(length(trange))
    result_Φbar = zeros(length(trange))
    sol = initial_guess
    for (i, t) in enumerate(trange)
        res = nlsolve(gapeqns_pnjl(t, μ, param, pyp; regularization=regularization, pv_scale=pv_scale), sol)
        sol = res.zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
        result_Φ[i] = sol[3]
        result_Φbar[i] = sol[4]
    end
    return result_m, result_ω, result_Φ, result_Φbar
end

function massgap_pnjl(t, μrange::AbstractVector, param::Parameters, pyp::PolyakovParameters=PolyakovParameters(); initial_guess=[0.3, 0.01, 0.1, 0.1], regularization=:sharp, pv_scale=10.0)
    result_m = zeros(length(μrange))
    result_ω = zeros(length(μrange))
    result_Φ = zeros(length(μrange))
    result_Φbar = zeros(length(μrange))
    sol = initial_guess
    for (i, μ) in enumerate(μrange)
        res = nlsolve(gapeqns_pnjl(t, μ, param, pyp; regularization=regularization, pv_scale=pv_scale), sol)
        sol = res.zero
        result_m[i] = sol[1]
        result_ω[i] = sol[2]
        result_Φ[i] = sol[3]
        result_Φbar[i] = sol[4]
    end
    return result_m, result_ω, result_Φ, result_Φbar
end

export PolyakovParameters, massgap_pnjl
