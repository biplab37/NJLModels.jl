# relativistic density functional

function _get_rdf_couplings(condensate, condensate_0, a, D0)
    A = (1 + a) * condensate_0^2 - condensate^2
    cubrtA = cbrt(A)
    L0 = D0 * cubrtA
    Sigma_MF = 2 * D0 * condensate / (3 * cubrtA^2)
    G_S = D0 * (3 * A + 4 * condensate^2) / (9 * cubrtA^5)
    G_PS = D0 / (3 * cubrtA^2)

    return L0, Sigma_MF, G_S, G_PS
end

function _get_G_PS(condensate, condensate_0, a, D0)
    A = (1 + a) * condensate_0^2 - condensate^2

    return D0 / (3 * cbrt(A)^2)
end

function _Ener(p, m, mu)
    return sqrt(p^2 + m^2) - mu
end

## vacuum gap equations 

function _condensate_integrand(condensate, p, m0, G_PS)
    m = m0 + 2 * G_PS * condensate
    #  return (12 / pi^2) * p^2 * sign(m)
    return (12 / (2 * pi^2)) * p^2 * m / sqrt(p^2 + m^2)
end

function _vacuum_gap(condensate, params)
    m0, a, D0, cutoff = params
    G_PS = _get_G_PS(condensate, condensate, a, D0)
    integrand(p) = _condensate_integrand(condensate, p, m0, G_PS)
    return integrate(p -> integrand(p), 0, cutoff)
end

function _get_vacuum_condensate(params)
    return fzero(condensate -> condensate - _vacuum_gap(condensate, params), 0.02)
end

function _get_pion_pol(condensate, params, G_PS)
    m0, a, D0, cutoff = params
    m = m0 + 2 * G_PS * condensate
    I(k, M) = 12 * k^2 / ((2 * π^2) * En(k, m)) * PrincipalValue(4 * En(k, m)^2 - M^2)

    pion_pole(Mpi) = Mpi^2 - (1 / (2 * G_PS) - condensate / m) / integrate(x -> I(x, Mpi), 0, cutoff)

    return fzero(pion_pole, 0.1)
end

function _get_pion_pol(params)
    condensate = _get_vacuum_condensate(params)
    G_PS = _get_G_PS(condensate, condensate, params[2], params[3])
    return _get_pion_pol(condensate, params, G_PS)
end

function _get_pion_pole(params)
    @show condensate = _get_vacuum_condensate(params)
    @show G_PS = _get_G_PS(condensate, condensate, params[2], params[3])
    @show m = params[1] + 2 * G_PS * condensate
    param = Parameters1(m0=params[1], Λ=params[4], Gs=G_PS)
    rep(x) = realpart_meson_q0(0.01, 0.0, x, m, 1 / 2, param)

    return fzero(rep, 0.1)
end

function _get_f_pi(params)
    Mpi = _get_pion_pol(params)
    @show condensate = _get_vacuum_condensate(params)
    # m = params[1] + _get_G_PS(condensate, condensate, params[2], params[3])
    return sqrt(params[1] * condensate / Mpi^2)
end

function _find_parameters(; condensate=0.039, mpi=0.135, fpi=92, mass=0.4)
    function func!(F, params)
        F[1] = condensate - _get_vacuum_condensate(params)
        F[2] = mpi - _get_pion_pol(params)
        F[3] = fpi - _get_f_pi(params)
        F[4] = mass - params[1] - _get_G_PS(condensate, condensate, params[2], params[3]) * condensate
    end

    return mcpsolve(func!, [0.004, 0.5, 0.1, 0.56], [0.006, 2.0, 2.0, 0.7], [0.005, 2., 0.5, 0.6])
end

function _find_mass(params)
    condensate = _get_vacuum_condensate(params)
    G_PS = _get_G_PS(condensate, condensate, params[2], params[3])
    return params[1] + 2 * G_PS * condensate
end

## medium values
abstract type ParametersRDF end

Base.@kwdef struct ParametersRDF1 <: ParametersRDF
    m0::Float64 = 0.0055
    a::Float64 = 1.5
    D0::Float64 = 0.2
    Λ::Float64 = 0.62
end
Base.@kwdef struct ParametersRDF2 <: ParametersRDF
    m0::Float64 = 0.0042
    a::Float64 = 1.43
    D0::Float64 = 1.39 * (0.573)^2
    Λ::Float64 = 0.573
end

Base.@kwdef struct ParametersRDF3 <: ParametersRDF
    m0::Float64 = 0.0055
    a::Float64 = 3.43
    D0::Float64 = 0.897 * (0.62)^2
    Λ::Float64 = 0.62
end

Base.@kwdef struct ParametersRDF4 <: ParametersRDF
    m0::Float64 = 0.0055
    a::Float64 = 2.1
    D0::Float64 = 0.2955
    Λ::Float64 = 0.59
end

function _get_vacuum_condensate(params::ParametersRDF)
    return fzero(condensate -> condensate - _vacuum_gap(condensate, [params.m0, params.a, params.D0, params.Λ]), 0.05)
end


function _int_condensate_rdf_pnjl(T, mu, condensate, Phi, Phibar, param::ParametersRDF)
    condensate_0 = _get_vacuum_condensate(param)
    G_PS = _get_G_PS(condensate, condensate_0, param.a, param.D0)
    m = param.m0 + 2 * G_PS * condensate
    vacuum_part_integrand(p) = (12 / (2 * pi^2)) * p^2 * m / sqrt(p^2 + m^2)
    vacuum_part = integrate(vacuum_part_integrand, 0, param.Λ)
    medium_part_integrand(p) = -(12 / (2 * pi^2)) * p^2 * m * (f_polyakov(T, mu, En(p, m), Phi, Phibar) + fbar_polyakov(T, mu, En(p, m), Phi, Phibar)) / sqrt(p^2 + m^2)
    # medium_part_integrand(p) = -(12 / (2 * pi^2)) * p^2 * m * (FD_dist(T, mu, En(p, m)) + FD_dist(T, -mu, En(p, m))) / sqrt(p^2 + m^2)
    medium_part = integrate(medium_part_integrand, 0, 2.0)
    return vacuum_part + medium_part
end

function integrand_μ_pnjl(T, μ, m, ω, Φ, Φbar, param::Parameters)
    μ_star = μ + ω
    return p -> 12 * param.Gv * (p^2 / π^2) *
                (f_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar) - fbar_polyakov(T, μ_star, safe_En(p, m), Φ, Φbar))
end


function _gap_equations_rdf(T, mu, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters())
    function gap(condensate, Phi, Phibar)
        if isnan(condensate) || isnan(Phi) || isnan(Phibar) || isinf(condensate) || isinf(Phi) || isinf(Phibar)
            return [1e4, 1e4, 1e4]
        end

        c_c = clamp(condensate, 0.0001, 0.06)
        Phi_c = clamp(Phi, 0.0, 1.0)
        Phibar_c = clamp(Phibar, 0.0, 1.0)
        condensate_0 = _get_vacuum_condensate(param)
        G_PS = _get_G_PS(c_c, condensate_0, param.a, param.D0)
        m = param.m0 + 2 * G_PS * c_c

        m_c = clamp(m, 0.0, 5.0)
        dU_dΦ, dU_dΦbar = U_Polyakov_derivs(T, Phi_c, Phibar_c, pyp)
        return [
            m_c - param.m0 - 2 * G_PS * _int_condensate_rdf_pnjl(T, mu, c_c, Phi_c, Phibar_c, param),
            dU_dΦ - integrate(integrand_Φ_pnjl(T, mu, m_c, 0.0, Phi_c, Phibar_c), 0, 2.0),
            dU_dΦbar - integrate(integrand_Φbar_pnjl(T, mu, m_c, 0.0, Phi_c, Phibar_c), 0, 2.0)
        ]
    end
    return x -> gap(x[1], x[2], x[3])
end

function massgap_rdf(T, μ, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(); initial_guess=[0.02, 0.01, 0.01])
    return mcpsolve(_gap_equations_rdf(T, μ, param, pyp), [-0.01, 0.0, 0.0], [0.02818, 1.0, 1.0], initial_guess, inplace=false)
end

function massgap_rdf(trange::AbstractVector, μ, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(); initial_guess=[0.03, 0.01, 0.1])
    result_m = zeros(length(trange))
    result_Φ = zeros(length(trange))
    result_Φbar = zeros(length(trange))
    sol = initial_guess
    for (i, t) in enumerate(trange)
        res = nlsolve(_gap_equations_rdf(t, μ, param, pyp), sol)
        sol = res.zero
        result_m[i] = sol[1]
        result_Φ[i] = sol[2]
        result_Φbar[i] = sol[3]
    end
    return result_m, result_Φ, result_Φbar
end

function massgap_rdf(t, μrange::AbstractVector, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(); lower=[0.0001, 0.0, 0.0], upper=[0.029, 1.0, 1.0], initial_guess=[0.02, 0.01, 0.01])
    result_m = zeros(length(μrange))
    result_Φ = zeros(length(μrange))
    result_Φbar = zeros(length(μrange))
    sol = initial_guess
    for (i, μ) in enumerate(μrange)
        res = mcpsolve(_gap_equations_rdf(t, μ, param, pyp), lower, upper, sol, inplace=false)
        sol = res.zero
        result_m[i] = sol[1]
        result_Φ[i] = sol[2]
        result_Φbar[i] = sol[3]
    end
    return result_m, result_Φ, result_Φbar
end

function _get_mass(condensate, param::ParametersRDF)
    condensate_0 = _get_vacuum_condensate(param)
    G_PS = _get_G_PS(condensate, condensate_0, param.a, param.D0)
    return param.m0 + 2 * G_PS * condensate
end

function _get_condensate(T, mu, param::ParametersRDF)
    gap(con) = con - _int_condensate_rdf_pnjl(T, mu, con, 0.0, 0.0, param)
    return return gap
end

function _get_mass(T, mu, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(T0=0.208))
    sol = massgap_rdf(T, mu, param, pyp).zero
    return _get_mass(sol[1], param)
end

function gap_rdf(T, mu, param, pyp=PolyakovParameters(T0=0.208))
    sol = massgap_rdf(T, mu, param, pyp).zero
    m = _get_mass(sol[1], param)

    return [m, sol[2], sol[3]]
end

function _get_density(T, mu, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(T0=0.208))
    sol = massgap_rdf(T, mu, param, pyp).zero
    m = _get_mass(sol[1], param)
    integrand(p) = (12 / (2 * pi^2)) * p^2 * (f_polyakov(T, mu, En(p, m), sol[2], sol[3]) - fbar_polyakov(T, mu, En(p, m), sol[2], sol[3]))
    return integrate(integrand, 0.0, 2.0)
end
function _get_density(T, mu, m, Phi, Phibar)
    integrand(p) = (12 / (2 * pi^2)) * p^2 * (f_polyakov(T, mu, En(p, m), Phi, Phibar) - fbar_polyakov(T, mu, En(p, m), Phi, Phibar))
    return integrate(integrand, 0.0, 2.0)
end

function U_Polyakov(T, Phi, Phibar, pyp::PolyakovParameters=PolyakovParameters(T0=0.208))
    b2 = pyp.a0 + pyp.a1 * (pyp.T0 / T) + pyp.a2 * (pyp.T0 / T)^2 + pyp.a3 * (pyp.T0 / T)^3

    return (-0.5*b2*Phi*Phibar - (pyp.b3/6)*(Phi^3 + Phibar^3) + (pyp.b4/4)*(Phi*Phibar)^2)*T^4
end

function _Omega(T, mu, con, Phi, Phibar, param::ParametersRDF, pyp::PolyakovParameters=PolyakovParameters(T0=0.208))
    condensate_0 = _get_vacuum_condensate(param)
    G_PS = _get_G_PS(con, condensate_0, param.a, param.D0)
    m = param.m0 + 2 * G_PS * con
    Omega_vacuum = -12*integrate(p->p^2*En(p, m), 0.0, param.Λ)/(2*pi^2)
    Omega_quark = 12*integrate(p->p^2*(_log_Z_phi_minus(T, mu, En(p, m), Phi, Phibar) + _log_Z_phi_plus(T, mu, En(p, m), Phi, Phibar)), 0.0, 2.0)/(2*pi^2)

    U_pol = U_Polyakov(T, Phi, Phibar, pyp)
    L_RDF0 = param.D0*cbrt((1 + param.a) * condensate_0^2 - con^2)

    return Omega_vacuum + Omega_quark + U_pol + L_RDF0 - 2*G_PS*con^2
end