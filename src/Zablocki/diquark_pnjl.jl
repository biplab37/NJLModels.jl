function Π0_D(T, mu, m, param, Phi, Phibar)
    factor = (3 - 1) * 2 / (2 * π^2)
    function integrand_pi0_vacuum(ep)
        return sqrt(ep^2 - m^2) * ep * (PrincipalValue(ep + mu) + PrincipalValue(ep - mu))
    end
    function integrand_pi0_medium(ep)
        return sqrt(ep^2 - m^2) * ep * (((-2 * FD_dist(T, -mu, ep, Phi, Phibar)) * PrincipalValue(ep + mu)) + (-2 * FD_dist(T, mu, ep, Phi, Phibar)) * PrincipalValue(ep - mu))
    end

    return 1 / (2 * param.GD) - factor * integrate(integrand_pi0_vacuum, m, sqrt(param.Λ^2 + m^2)) - factor * integrate(integrand_pi0_medium, m, 3.0)
end

_NN(T, mu, ω, q, m, sign_up, sign_down, sign_fermi, sign_energy, Phi, Phibar) = FD_dist(T, -sign_down * sign_fermi * mu, sign_energy * E_D(sign_up, sign_down, ω, q, mu, m), Phi, Phibar)

function J_D_pair(sign::Int, T, mu, ω, q, m, Phi, Phibar)
    # mu = 0.0
    return T * log((_NN(T, mu, ω, q, m, -1, sign, +1, -1, Phi, Phibar) * _NN(T, mu, ω, q, m, -1, sign, -1, +1, Phi, Phibar)) / (_NN(T, mu, ω, q, m, +1, sign, +1, -1, Phi, Phibar) * _NN(T, mu, ω, q, m, +1, sign, -1, +1, Phi, Phibar)))
end

function J_D_Landau(sign::Int, T, mu, ω, q, m, Phi, Phibar)
    # mu = 0.0
    return 2 * T * log(_NN(T, mu, ω, q, m, -1, sign, -1, +1, Phi, Phibar) / _NN(T, mu, ω, q, m, +1, sign, +1, -1, Phi, Phibar))
end

function imagpart_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    factor = (3 - 1) * 2 / (8π)
    # ω1 = ω + 2 * mu
    positive_term = 0.0
    if ω * ((ω / 4) + mu) > m^2 - mu^2
        positive_term = (ω + 2mu)^2 * sqrt(1 - 4 * m^2 / (ω + 2mu)^2)
    end

    if (ω + 2 * mu)^2 > 4 * (param.Λ^2 + m^2)
        # return factor * (-2 * FD_dist(T, 0.0, ω / 2, Phi, Phibar)) * (positive_term)
        return 0.0
    end

    # negative_term = 0.0
    # if ω * ((ω / 4) - mu) > m^2 - mu^2
    #     negative_term = (ω - 2mu)^2 * sqrt(1 - 4 * m^2 / (ω - 2mu)^2)
    # end

    return factor * (1 - 2 * FD_dist(T, 0.0, ω / 2, Phi, Phibar)) * (positive_term) #+ negative_term)
end

function imagpart_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
    # ω = ω1 + 2 * mu
    if (ω + 2mu)^2 > 4 * (param.Λ^2 + m^2 + (q^2 / 4))
        return 0.0
    end

    if q == 0.0
        return imagpart_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    end

    factor = 2 / (4 * π * q)

    # postive branch
    sp = s_D(ω, q, mu, +1)

    if sp > 4m^2
        positive = sp * J_D_pair(+1, T, mu, ω, q, m, Phi, Phibar)
    elseif sp < 0
        positive = sp * J_D_Landau(+1, T, mu, ω, q, m, Phi, Phibar)
    else
        positive = 0.0
    end

    #negative branch
    # sm = s_D(ω, q, mu, -1)

    # if sm > 4m^2
    #     negative = sm * J_D_pair(-1, T, mu, ω, q, m)
    # elseif sm < 0
    #     negative = sm * J_D_Landau(-1, T, mu, ω, q, m)
    # else
    #     negative = 0.0
    # end

    return factor * (positive)# + negative)
end

function realpart_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
    impart(x, y) = imagpart_D_normal(T, mu, x, y, m, param, Phi, Phibar)

    cutoff = sqrt(param.Λ^2 + m^2 + (q^2 / 4))
    repart_dependent = realpart_kramers_kronig_q_1(impart, ω, q, -2 * (cutoff + mu), 2 * (cutoff - mu))

    return Π0_D(T, mu, m, param, Phi, Phibar) - repart_dependent
end

function realpart_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    impart(x) = imagpart_D_normal_q0(T, mu, x, m, param, Phi, Phibar)

    cutoff = 2.0
    repart_dependent = realpart_kramers_kronig(impart, ω, cutoff)

    return Π0_D(T, mu, m, param, Phi, Phibar) - repart_dependent
end

function Πq0_D_analytic(T, mu, ω, m, param, Phi, Phibar)
    if imag(ω) == 0
        return realpart_D_normal_q0(T, mu, real(ω), m, param, Phi, Phibar) #- 1im * imagpart_D_normal_q0(T, mu, real(ω), m, param, Phi, Phibar)
    end
    integrand_vac(ep) = 4 * sqrt(ep^2 - m^2) * ep * ((-1) / (ω + 2 * mu + 2 * ep) + (1) / (ω + 2 * mu - 2 * ep)) / π^2
    integrand_med(ep) = 4 * sqrt(ep^2 - m^2) * ep * ((2 * FD_dist(T, -mu, ep, Phi, Phibar)) / (ω + 2 * mu + 2 * ep) + (-2 * FD_dist(T, mu, ep, Phi, Phibar)) / (ω + 2 * mu - 2 * ep)) / π^2
    imagpart = (imag(ω) >= 0.0) ? 0.0 : 1im * (1 - 2 * FD_dist(T, 0.0, ω / 2)) * (ω + 2 * mu)^2 * sqrt(1 - 4 * m^2 / (ω + 2 * mu)^2) / (2 * π)
    # int_vac = integrand_vac(m, 3.0) - 2 * integrand_vac(sqrt(m^2 + param.Λ^2), 3.0) + integrand_vac(sqrt(m^2 + 2 * param.Λ^2), 3.0)
    return 1 / (2 * param.GD) + integrate(integrand_vac, m, sqrt(param.Λ^2 + m^2)) + integrate(integrand_med, m, 2.0) - 2imagpart
end

function find_mass_D_q0(T, mu, m, param, Phi, Phibar, guess=[0.6, 0.01])
    rep(ω) = realpart_D_normal_q0(T, mu, ω - 2mu, m, param, Phi, Phibar)
    if m > mu && rep(0.0) * rep(2(m)) < 0.0
        return [bisection(rep, 0.0, 2(m)) - 2mu, 0.0]
    end
    function ff!(F, x)
        term = Πq0_D_analytic(T, mu, x[1] - 1im * x[2] / 2 - 2mu, m, param, Phi, Phibar)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return mcpsolve(ff!, [0.0, 0.0], [1.0, 1.0], guess).zero - [2mu, 0.0]
end

function spectral_function_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    impart = imagpart_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    if abs(impart) <= 1e-4
        return 0.0
    end
    repart = real(Πq0_D_analytic(T, mu, ω, m, param, Phi, Phibar))

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
    impart = imagpart_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
    if abs(impart) <= 1e-4
        return 0.0
    end
    repart = realpart_D_normal(T, mu, ω, q, m, param, Phi, Phibar)

    return impart / (π * (repart^2 + impart^2))
end

function _has_diquark_bound_state_q0(T, mu, m, param, Phi, Phibar)
    rep(ω) = realpart_D_normal_q0(T, mu, ω, m, param, Phi, Phibar)
    return _has_bound_state(rep, 2 * m)
end
function _has_diquark_bound_state(T, mu, q, m, param, Phi, Phibar)
    rep(ω) = realpart_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
    return _has_bound_state(rep, 2 * sqrt(m^2 + q^2/4))
end

function wave_function_renormalization_diquark(T, mu, m, ed, param, Phi, Phibar, tol=1e-6)
    return tol / abs(realpart_D_normal_q0(T, mu, ed + tol, m, param, Phi, Phibar) - realpart_D_normal_q0(T, mu, ed, m, param, Phi, Phibar))
end

function wave_function_renormalization_diquark(T, mu, q, m, ed, param, Phi, Phibar, tol=1e-6)
    return tol / abs(realpart_D_normal(T, mu, ed + tol, q, m, param, Phi, Phibar) - realpart_D_normal(T, mu, ed, q, m, param, Phi, Phibar))
end

function _f_sum_diquark(T, mu, m, param, Phi, Phibar)
    M = _has_diquark_bound_state_q0(T, mu, m, param, Phi, Phibar)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2 * M * wave_function_renormalization_diquark(T, mu, m, M, param, Phi, Phibar)
    end

    _fsum += integrate(o -> 2 * o * spectral_function_D_normal_q0(T, mu, o, m, param, Phi, Phibar), 0.0, 2.0)
    return _fsum
end

function _f_sum_diquark(T, mu, q, m, param, Phi, Phibar)
    M = _has_diquark_bound_state(T, mu, q, m, param, Phi, Phibar)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2 * M * wave_function_renormalization_diquark(T, mu, q, m, M, param, Phi, Phibar)
    end

    _fsum += integrate(o -> 2 * o * spectral_function_D_normal(T, mu, o, q, m, param, Phi, Phibar), 0.0, 2*(sqrt(param.Λ^2 + m^2 + q^2/4)))
    return _fsum
end
