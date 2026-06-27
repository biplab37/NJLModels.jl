# diquarks in PNJL + RDF model

function Π0_D_rdf(T, mu, m, param, Phi, Phibar)
    factor = (3 - 1) * 2 / (2 * π^2)
    function integrand_pi0_vacuum(ep)
        return sqrt(ep^2 - m^2) * ep * (PrincipalValue(ep + mu) + PrincipalValue(ep - mu))
    end
    function integrand_pi0_medium(ep)
        return sqrt(ep^2 - m^2) * ep * ((-2*fbar_polyakov(T, mu, ep, Phi, Phibar)) * PrincipalValue(ep + mu) + (-2*f_polyakov(T, mu, ep, Phi, Phibar)) * PrincipalValue(ep - mu))
    end
    return 1 / (2 * param.GD) - factor * integrate(integrand_pi0_vacuum, m, sqrt(param.Λ^2 + m^2)) - factor * integrate(integrand_pi0_medium, m, 3.0)
end

function J_D_pair_rdf(xi::Int, T, mu, eplus, eminus, Phi, Phibar)
    Z_plus(e) = _log_Z_phi_plus(T, mu, e, Phi, Phibar)
    Z_minus(e) = _log_Z_phi_minus(T, mu, e, Phi, Phibar)
    # mu = 0.0
    return (Z_plus(xi*eplus) + Z_minus(-xi*eplus) - Z_plus(xi*eminus) - Z_minus(-xi*eminus))/3
end

function J_D_Landau_rdf(T, mu, eplus, eminus, Phi, Phibar)
    # mu = 0.0
    Z_plus(e) = _log_Z_phi_plus(T, mu, e, Phi, Phibar)
    Z_minus(e) = _log_Z_phi_minus(T, mu, e, Phi, Phibar)
    return -(Z_minus(eplus) + Z_plus(-eplus) - Z_minus(eminus) - Z_plus(-eminus))/3
end

function imagpart_D_normal_q0_C_rdf(T, mu, ω, m, param, Phi, Phibar)
    factor = (3 - 1) * 2 / (8π)
    z = ω + 2 * mu
    positive_term = z^2 * sqrt(1 - 4 * m^2 / z^2)

    return factor * (1 - f_polyakov(T, mu, z / 2, Phi, Phibar) - fbar_polyakov(T, mu, z / 2, Phi, Phibar)) * (positive_term)
end

function imagpart_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)
    z = ω + 2 * mu
    if z^2 < 4m^2 || z^2 > 4 * (param.Λ^2 + m^2)
        return 0.0
    end
    return imagpart_D_normal_q0_C_rdf(T, mu, ω, m, param, Phi, Phibar)
end

function imagpart_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    # ω = ω1 + 2 * mu
    z = ω + 2 * mu
    if (z)^2 > 4 * (param.Λ^2 + m^2 + (q^2 / 4))
        return 0.0
    end

    if q == 0.0
        return imagpart_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)
    end

    factor = 2 / (4 * π * q)

    # postive branch
    sp = z^2 - q^2

    if sp > 4m^2
        eplus = 0.5*(z + q*sqrt(1 - 4 * m^2 / sp))
        eminus = 0.5*(z - q*sqrt(1 - 4 * m^2 / sp))
        if z>0
            positive = sp * J_D_pair_rdf(+1, T, mu, eplus, eminus, Phi, Phibar)
        else
            positive = sp * J_D_pair_rdf(-1, T, mu, eplus, eminus, Phi, Phibar)
        end
    elseif sp < 0
        eplus = 0.5*(z + q*sqrt(1 - 4 * m^2 / sp))
        eminus = 0.5*(z - q*sqrt(1 - 4 * m^2 / sp))
        positive = sp * J_D_Landau_rdf(T, mu, eplus, eminus, Phi, Phibar)
    else
        positive = 0.0
    end

    return factor * (positive)
end

function realpart_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    impart(x, y) = imagpart_D_normal_rdf(T, mu, x, y, m, param, Phi, Phibar) 

    cutoff = sqrt(param.Λ^2 + m^2 + (q^2 / 4))
    repart = realpart_kramers_kronig_q_1(impart, ω, q, -2*cutoff-2mu, 2*cutoff-2mu)

    return Π0_D_rdf(T, mu, m, param, Phi, Phibar) - repart
end

function realpart_D_normal_q0_rdf(T, mu, z, m, param, Phi, Phibar)
    impart_vac(x) = imagpart_D_normal_q0_rdf(0.0, 0.0, x, m, param, Phi, Phibar)
    impart(x) = imagpart_D_normal_q0_rdf(T, mu, x, m, param, Phi, Phibar)  - impart_vac(x)

    cutoff = sqrt(param.Λ^2 + m^2)
    repart = realpart_kramers_kronig_1(impart_vac, z, -2*cutoff - 2mu, 2*cutoff- 2mu)
    repart_med = realpart_kramers_kronig_1(impart, z, -3.0, 3.0)

    return Π0_D_rdf(T, mu, m, param, Phi, Phibar) - repart - repart_med
end

function Πq0_D_analytic_rdf(T, mu, z, m, param, Phi, Phibar)
    if imag(z) == 0
        return realpart_D_normal_q0_rdf(T, mu, real(z), m, param, Phi, Phibar) - 1im * imagpart_D_normal_q0_rdf(T, mu, real(z), m, param, Phi, Phibar)
    end
    integrand_vac(ep) = 4 * sqrt(ep^2 - m^2) * ep * ((-1) / (z +2 * mu + 2 * ep) + (1) / (z + 2 * mu - 2 * ep)) / π^2
    integrand_med(ep) = 4 * sqrt(ep^2 - m^2) * ep * ((2 * fbar_polyakov(T, mu, ep, Phi, Phibar)) / (z + 2 * mu + 2 * ep) + (-2 * f_polyakov(T, mu, ep, Phi, Phibar)) / (z + 2 * mu - 2 * ep)) / π^2
    imagpart = (imag(z) >= 0.0) ? 0.0 : 1im * imagpart_D_normal_q0_C_rdf(T, mu, z, m, param, Phi, Phibar)
    # int_vac = integrand_vac(m, 3.0) - 2 * integrand_vac(sqrt(m^2 + param.Λ^2), 3.0) + integrand_vac(sqrt(m^2 + 2 * param.Λ^2), 3.0)
    return 1 / (2 * param.GD) + integrate(integrand_vac, m, sqrt(param.Λ^2 + m^2)) + integrate(integrand_med, m, 2.0) - 2imagpart
end

function full_real_part_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    function integrand(p, x)
        Em = sqrt(p^2 + q^2/4 + m^2 - 2p*q*x)
        Ep = sqrt(p^2 + q^2/4 + m^2 + 2p*q*x)
        tpm(xi, xip) = 1 - xi*xip*(p^2 + q^2/4 + m^2)/(Em*Ep)
        denom(xi, xip) = f_polyakov(T, mu, xi*Em, Phi, Phibar) - f_polyakov(T, mu, xip*Ep, Phi, Phibar)
        numerator(xi, xip) = ω + xi*Em - xip*Ep
        return p^2 * sum(xi -> sum(xip -> tpm(xi, xip)*denom(xi, xip)/numerator(xi, xip), [-1, 1]), [-1, 1])/(2π^2)
    end
    return 1 / (2 * param.GD) + integrate(x->integrand(x[1], x[2]), [0.0, -1.0], [param.Λ, 1.0])
end

function find_mass_D_q0_rdf(T, mu, m, param, Phi, Phibar, guess=[0.7, 0.01])
    rep(z) = realpart_D_normal_q0_rdf(T, mu, z, m, param, Phi, Phibar)
    if m > mu && rep(0.0) * rep(2(m) - 2mu) < 0.0
        return [bisection(rep, 0.0, 2(m)-2mu) , 0.0]
    end
    function ff!(F, x)
        term = Πq0_D_analytic_rdf(T, mu, x[1] - 1im * x[2] / 2, m, param, Phi, Phibar)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return mcpsolve(ff!, [0.0, 0.0], [1.0, 1.0], guess).zero
end

function spectral_function_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)
    impart = imagpart_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)
    if abs(impart) <= 1e-4
        return 0.0
    end
    repart = realpart_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)

    return impart / (π * (repart^2 + impart^2))
end

function spectral_function_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    impart = imagpart_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    if abs(impart) <= 1e-4
        return 0.0
    end
    repart = realpart_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)

    return impart / (π * (repart^2 + impart^2))
end

function _has_diquark_bound_state_q0_rdf(T, mu, m, param, Phi, Phibar)
    rep(ω) = realpart_D_normal_q0_rdf(T, mu, ω, m, param, Phi, Phibar)
    return _has_bound_state(rep, 2 * m)
end
function _has_diquark_bound_state_rdf(T, mu, q, m, param, Phi, Phibar)
    rep(ω) = realpart_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    return _has_bound_state(rep, 2 * sqrt(m^2 + q^2/4))
end

function wave_function_renormalization_diquark_q0_rdf(T, mu, m, ed, param, Phi, Phibar, tol=1e-6)
    return tol / abs(realpart_D_normal_q0_rdf(T, mu, ed + tol, m, param, Phi, Phibar) - realpart_D_normal_q0_rdf(T, mu, ed, m, param, Phi, Phibar))
end

function wave_function_renormalization_diquark_rdf(T, mu, q, m, ed, param, Phi, Phibar, tol=1e-6)
    return tol / abs(realpart_D_normal_rdf(T, mu, ed + tol, q, m, param, Phi, Phibar) - realpart_D_normal_rdf(T, mu, ed, q, m, param, Phi, Phibar))
end

function _f_sum_diquark_rdf(T, mu, m, param, Phi, Phibar)
    M = _has_diquark_bound_state_q0_rdf(T, mu, m, param, Phi, Phibar)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2 * M * wave_function_renormalization_diquark_q0_rdf(T, mu, m, M, param, Phi, Phibar)
    end

    _fsum += integrate(o -> 2 * o * spectral_function_D_normal_q0_rdf(T, mu, o, m, param, Phi, Phibar), 0.0, max(10.0*T, 2.0))
    return _fsum
end

function _f_sum_diquark_zeroth_rdf(T, mu, m, param, Phi, Phibar)
    M = _has_diquark_bound_state_q0_rdf(T, mu, m, param, Phi, Phibar)
    _fsum = 0.0

    if M > 0.0
        _fsum += wave_function_renormalization_diquark_q0_rdf(T, mu, m, M, param, Phi, Phibar)
    end

    _fsum += integrate(o -> spectral_function_D_normal_q0_rdf(T, mu, o, m, param, Phi, Phibar), 0.0, max(10.0*T, 5.0))
    return _fsum
end

function _f_sum_diquark_rdf(T, mu, q, m, param, Phi, Phibar)
    M = _has_diquark_bound_state_rdf(T, mu, q, m, param, Phi, Phibar)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2 * M * wave_function_renormalization_diquark_rdf(T, mu, q, m, M, param, Phi, Phibar)
    end

    _fsum += integrate(o -> 2 * o * spectral_function_D_normal_rdf(T, mu, o, q, m, param, Phi, Phibar), 0.0, 2.0)
    return _fsum
end
