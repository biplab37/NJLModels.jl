function _integ_s_cor_rdf(spec, T, mu, k, p0, m, param, Phi, Phibar)
    term1 = abs(p0 - En(k, m)) < 1e-7 ? 0.0 : spec(p0 - En(k, m), k) * (1 - f_polyakov(T, mu, En(k, m), Phi, Phibar) + numberB(T, -mu, p0 - En(k, m)))
    term2 = abs(p0 + En(k, m)) < 1e-7 ? 0.0 : -spec(p0 + En(k, m), k) * (fbar_polyakov(T, mu, En(k, m), Phi, Phibar) + numberB(T, -mu, p0 + En(k, m)))
    return k^2 * m * (term1 + term2) / (π * En(k, m)) ## -1 factor..check calc in goodnotes
end

function _imagpart_baryons_spectral_normal_q0_s_cor_rdf(spec, T, mu, p0, m, param::Parameters, Phi, Phibar; tol=1e-5, maxevals=1e7)
    return quadgk(k -> _integ_s_cor_rdf(spec, T, mu, k, p0, m, param, Phi, Phibar), 0.0, param.Λ; atol=tol, rtol=tol, maxevals=maxevals)[1]
end

function _get_interpolated_spectral_rdf(T, mu, m, param, Phi, Phibar)
    order = (500, 100)
    lb_qp = [0.0, 0.0]
    lb_ld = [0.0, 2 * mu]
    ub = [1.0, 2 * param.Λ + 0.1]

    # Landau Damping: Maps [u, q] back to 0 < ω < q - 2μ
    interp_ld = chebinterp(order, lb_ld, ub) do x
        u, q = x[1], x[2]
        ω_lower, ω_upper = -q - 2 * mu, q - 2 * mu
        ω = ω_lower + u * (ω_upper - ω_lower)
        return spectral_function_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    end
    # Quasiparticle: Maps [t, q] back to ω_lower < ω < ω_upper
    interp_qp = chebinterp(order, lb_qp, ub) do x
        t, q = x[1], x[2]
        ω_lower, ω_upper = _get_threshold_cutoff(q, m, mu, param.Λ)

        ω = ω_lower + t * (ω_upper - ω_lower)
        return spectral_function_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    end

    interp_qp2 = chebinterp(order, lb_qp, ub) do x
        t, q = x[1], x[2]
        ω_lower, ω_upper = -sqrt(4 * m^2 + q^2) - 2 * mu, -sqrt(4 * param.Λ^2 + 4 * m^2 + q^2) - 2 * mu

        ω = ω_lower + t * (ω_upper - ω_lower)
        return spectral_function_D_normal_rdf(T, mu, ω, q, m, param, Phi, Phibar)
    end

    function spec(ω, q)
        if q > 2 * mu && -q - 2 * mu < ω < q - 2 * mu
            u = (ω - (-q - 2 * mu)) / (2 * q)
            return interp_ld([u, q])
        end

        ω_lower_qp, ω_upper_qp = _get_threshold_cutoff(q, m, mu, param.Λ)
        if ω_lower_qp < ω < ω_upper_qp
            t = (ω - ω_lower_qp) / (ω_upper_qp - ω_lower_qp)
            return interp_qp([t, q])
        end

        ω_lower_qp2, ω_upper_qp2 = -sqrt(4 * m^2 + q^2) - 2 * mu, -sqrt(4 * param.Λ^2 + 4 * m^2 + q^2) - 2 * mu
        if ω_lower_qp2 > ω > ω_upper_qp2
            t = (ω - ω_lower_qp2) / (ω_upper_qp2 - ω_lower_qp2)
            return interp_qp2([t, q])
        end

        return 0.0
    end
    return spec
end

function _integ_sigma_B_00_s_cor_interp_rdf(spec, T, mu, k, nu, m, param, Phi, Phibar)
    ek = En(k, m)
    stat_factor1 = (-1 + f_polyakov(T, mu, ek, Phi, Phibar) - numberB(T, -mu, nu)) / (ek + nu)
    stat_factor2 = (-fbar_polyakov(T, mu, ek, Phi, Phibar) - numberB(T, -mu, nu)) / (ek - nu)

    return k^2 * m * spec(nu, k) * (stat_factor1 + stat_factor2) / (π^2)
end

function sigma_B_00_s_cor_interp_rdf(spec, T, mu, m, param, Phi, Phibar)
    cutoff = Float64(param.Λ)

    lower_bound = [0.0, 0.0]
    upper_bound = [cutoff, 6 * cutoff]

    return integrate(x -> _integ_sigma_B_00_s_cor_interp_rdf(spec, T, mu, x[1], x[2], m, param, Phi, Phibar), lower_bound, upper_bound)
end

function find_diquark_energy_q_rdf(T, mu, q, m, param, Phi, Phibar)
    EQ = sqrt(q^2 + 4m^2)

    if EQ > mu && realpart_D_normal_rdf(T, mu, 0.0, q, m, param, Phi, Phibar) * realpart_D_normal_rdf(T, mu, EQ - 2mu, q, m, param, Phi, Phibar) < 0.0
        return bisection(x -> realpart_D_normal_rdf(T, mu, x, q, m, param, Phi, Phibar), 0.0, EQ - 2mu)
    end

    return 0.0
end
# function imagpart_baryon_q0_rdf(T, mu, ω, mq, mD, param, Phi, Phibar)
#     if ω == 0
#         return 0.0
#     end
#     eq = 0.5 * abs((ω^2 + mq^2 - mD^2) / ω)
#     eD = 0.5 * abs((ω^2 - mq^2 + mD^2) / ω)
#     # if eD < 1e-5
#     #     @show mq, mD, ω
#     # end
#     p = sqrt(abs(mq^4 + (mD^2 - ω^2)^2 - 2mq^2 * (mD^2 + ω^2))) / (2 * abs(ω))
#     factor = -p * mq / (π)
#     term1 = (1 - FD_dist(T, mu, eq, Phi, Phibar) + numberB(T, mu, eD)) / (eq + eD)
#     term2 = -(1 - FD_dist(T, -mu, eq, Phi, Phibar) + numberB(T, -mu, eD)) / (eq + eD)
#     term3 = -(FD_dist(T, -mu, eq, Phi, Phibar) + numberB(T, mu, eD)) * PrincipalValue(eD - eq)
#     term4 = (FD_dist(T, mu, eq, Phi, Phibar) + numberB(T, -mu, eD)) * PrincipalValue(eD - eq)
#     if abs(ω) > sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2)
#         return 0.0
#     end
#     if ω > 0
#         if ω > (mq + mD)
#             return factor * term2
#         elseif ω < (mD - mq)
#             return factor * term4
#         else
#             return 0.0
#         end
#     else
#         if ω < -(mq + mD)
#             return factor * term1
#         elseif ω > -(mD - mq)
#             return factor * term3
#         else
#             return 0.0
#         end
#     end
# end
#TODO:update this function to include the Polyakov loop distribution functions
function imagpart_baryon_q_rdf(T, mu, ω, q, mq, mD, param, Phi, Phibar)
    if q == 0
        return imagpart_baryon_q0_C_rdf(T, mu, ω, mq, mD, param, Phi, Phibar)
    end
    ω = ω + 3 * mu

    expb(eps) = exp(-(eps+2mu) / 2T)

    term1(xp, xm) = -T*(log((1 - expb(ω+xp))/(1 - expb(ω+xm)))) - (_log_Z_phi_plus(T, mu, -0.5*(ω-xp), Phi, Phibar) - _log_Z_phi_plus(T, mu, -0.5*(ω-xm), Phi, Phibar))/3

    term2(xp, xm) = +T*(log((1 - expb(ω-xp))/(1 - expb(ω-xm)))) + (_log_Z_phi_plus(T, mu, -0.5*(ω+xp), Phi, Phibar) - _log_Z_phi_plus(T, mu, -0.5*(ω+xm), Phi, Phibar))/3

    term3(xp, xm) = T*(log((1 - expb(ω + xp))/(1 - expb(ω-xm)))) + (_log_Z_phi_plus(T, mu, -0.5*(ω-xp), Phi, Phibar) - _log_Z_phi_plus(T, mu, -0.5*(ω+xm), Phi, Phibar))/3
    term4(xp, xm) = -T*(log((1 - expb(ω+xp))/(1 - expb(ω+xm)))) - (_log_Z_phi_plus(T, mu, -0.5*(ω-xp), Phi, Phibar) - _log_Z_phi_plus(T, mu, -0.5*(ω-xm), Phi, Phibar))/3
    term5(xp, xm) = +T*(log((1 - expb(ω-xp))/(1 - expb(ω-xm)))) + (_log_Z_phi_plus(T, mu, -0.5*(ω+xp), Phi, Phibar) - _log_Z_phi_plus(T, mu, -0.5*(ω+xm), Phi, Phibar))/3

    cutoff = max(sqrt(param.Λ^2 + mq^2 + q^2/4) + sqrt( param.Λ^2 + q^2/4 +  mD^2), sqrt(q^2 + 4*param.Λ^2 + (mq + mD)^2))
    factor = -4 * mq / (8 * π * q)

    s = ω^2 - q^2

    if abs(ω) >= cutoff || s == 0
        return 0.0
    end

    if mD > mq
        m1 = mD
        m2 = mq
    else
        m1 = mq
        m2 = mD
    end

    if s > 0
        if (m1 - m2)^2 <= s <= (m1 + m2)^2
            return 0.0
        end
        xm = x_pm(ω, q, m1, m2, -1)
        xp = x_pm(ω, q, m1, m2, +1)
        if s > (m1 + m2)^2
            if ω >= 0
                return -factor * term1(xp, xm)
            else
                return factor * term2(xp, xm)
            end
        else
            if ω >= 0
                return factor * term4(xp, xm) 
            else
                return factor * term5(xp, xm)
            end
        end
    else
        xm = x_pm(ω, q, m1, m2, -1)
        xp = x_pm(ω, q, m1, m2, +1)
        return factor * term3(xp, xm)
    end
end

function wave_function_renormalization_diquark_rdf(T, mu, q, m, ed, param, Phi, Phibar, tol=1e-6)
    return tol / abs(realpart_D_normal_rdf(T, mu, ed + tol, q, m, param, Phi, Phibar) - realpart_D_normal_rdf(T, mu, ed, q, m, param, Phi, Phibar))
end

function wave_function_renormalization_diquark_q0_rdf_1(T, mu, m, ed, param, Phi, Phibar, tol=1e-6)
    return 1/abs(UsefulFunctions._first_derivative_smooth(x -> realpart_D_normal_q0_rdf(T, mu, x, m, param, Phi, Phibar), ed, dx=tol))
end

function imagpart_baryons_spectral_normal_q0_s_bound_rdf(T, mu, p0, m, param, Phi, Phibar)
    md = find_diquark_energy_q_rdf(T, mu, 0.0, m, param, Phi, Phibar)
    if md == 0.0 || abs(p0)<1e-5
        return 0.0
    end
    Zk = wave_function_renormalization_diquark_q0_rdf_1(T, mu, m, md, param, Phi, Phibar)
    ed = 0.5 * abs((p0^2 - m^2 + md^2) / p0)

    return imagpart_baryon_q0_rdf(T, mu, p0, m, md, param, Phi, Phibar) * Zk * (2 * ed)
end

function _imagpart_baryons_spectral_normal_total_q0_interp_rdf(spec, T, mu, p0, m, param, Phi, Phibar)
    return _imagpart_baryons_spectral_normal_q0_s_cor_rdf(spec, T, mu, p0, m, param, Phi, Phibar) + imagpart_baryons_spectral_normal_q0_s_bound_rdf(T, mu, p0, m, param, Phi, Phibar)
end

function _realpart_baryons_spectral_normal_q0_s_interp_rdf(spec, T, mu, p0, m, param, Phi, Phibar)
    return realpart_kramers_kronig(x -> (_imagpart_baryons_spectral_normal_total_q0_interp_rdf(spec, T, mu, x, m, param, Phi, Phibar)), p0, 3 * sqrt(param.Λ^2 + m^2))
end

function _integ_sigma_B_00_s_bound_rdf(T, mu, k, m, param, Phi, Phibar)
    ed = find_diquark_energy_q_rdf(T, mu, k, m, param, Phi, Phibar)

    if ed == 0.0
        return 0.0
    end

    ek = En(k, m)
    Zk = wave_function_renormalization_diquark_rdf(T, mu, k, m, ed, param, Phi, Phibar)

    stat_factor1 = -(1 - FD_dist(T, mu, ek) + numberB(T, mu, ed)) / (ek + ed)
    stat_factor2 = -(FD_dist(T, mu, ek) + numberB(T, mu, ed)) / (ek - ed)

    return 2*k^2 * m * Zk * (stat_factor1 + stat_factor2) / (π^2 * ek)
end

function sigma_B_00_s_bound_rdf(T, mu, m, param, Phi, Phibar)
    cutoff = Float64(param.Λ)
    return integrate(k -> _integ_sigma_B_00_s_bound_rdf(T, mu, k, m, param, Phi, Phibar), 0.0, cutoff)
end

function p00_baryons_interp_rdf(spec, T, mu, m, mD, param, Phi, Phibar)
    1 / coupling_B_rdf(T, mu, m, mD, param, Phi, Phibar) - sigma_B_00_s_bound_rdf(T, mu, m, param, Phi, Phibar) - sigma_B_00_s_cor_interp_rdf(spec, T, mu, m, param, Phi, Phibar)
end
function p00_baryons_rdf(T, mu, m, mD, param, Phi, Phibar)
    1 / coupling_B_rdf(T, mu, m, mD, param, Phi, Phibar) - sigma_B_00_s_bound_rdf(T, mu, m, param, Phi, Phibar) - sigma_B_00_s_cor_rdf(T, mu, m, param, Phi, Phibar)
end

function _get_imag_baryon_total_q0_interp_rdf(T, mu, m, param, Phi, Phibar)
    func_cor(x) = imagpart_baryons_spectral_normal_q0_s_cor_rdf(spec, T, mu, x, m, param, Phi, Phibar)

    threshold_cor = 3 * (m - mu)
    cutoff_cor = 3 * (sqrt(param.Λ^2 + m^2) - mu)

    interp_cor = chebinterp(func_cor, 500, threshold_cor, cutoff_cor)

    cor_interp(x) = (threshold_cor < x < cutoff_cor) ? interp_cor(x) : 0.0

    func_bound(x) = imagpart_baryons_spectral_normal_q0_s_bound_rdf(T, mu, x, m, param, Phi, Phibar)

    mD = find_mass_D(T, mu, param)[1]
    threshold_bound = m + mD - 3 * mu
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + mD^2) - 3 * mu

    interp_bound = chebinterp(func_bound, 1000, threshold_bound, cutoff_bound)

    bound_interp(x) = (threshold_bound < x < cutoff_bound) ? interp_bound(x) : 0.0

    return x -> cor_interp(x) + bound_interp(x)
end

function _get_interpolated_imagpart_baryon_q0_rdf(spec, T, mu, m, param, Phi, Phibar)
    threshold_cor = 3m
    cutoff_cor = 3 * sqrt(param.Λ^2 + m^2)

    imag_cor = chebinterp(x -> _imagpart_baryons_spectral_normal_q0_s_cor_rdf(spec, T, mu, x, m, param, Phi, Phibar), 500, threshold_cor, cutoff_cor)
    function imag_cor_func(omega)
        if omega < 0
            return -imag_cor(-omega)
        end
        if threshold_cor < omega < cutoff_cor
            return imag_cor(omega)
        end
        return 0.0
    end
    md = find_diquark_energy_q_rdf(T, mu, 0.0, m, param, Phi, Phibar)

    threshold_bound_1, threshold_bound_2 = m + md, abs(m - md)
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + md^2)

    imag_bound_1 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_rdf(T, mu, x, m, param, Phi, Phibar), 500, threshold_bound_1, cutoff_bound)
    imag_bound_2 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_rdf(T, mu, x, m, param, Phi, Phibar), 500, 0.0, threshold_bound_2)

    function imag_bound_func(omega)
        if omega < 0
            return -imag_bound_func(-omega)
        end
        if 0.0 < omega < threshold_bound_2
            return imag_bound_2(omega)
        end
        if threshold_bound_1 < omega < cutoff_bound
            return imag_bound_1(omega)
        end
        return 0.0
    end

    return x -> imag_cor_func(x) + imag_bound_func(x)
end

## creating baryon with spectral function with finite momentum.

function _integ_s_cor_q_rdf(spec, T, mu, k, u, q, p0, m, param, Phi, Phibar)
    kq = sqrt(k^2 + q^2 + 2 * k * q * u)
    term1 = abs(p0 - En(k, m)) < 1e-7 ? 0.0 : spec(p0 - En(k, m), kq) * (1 - f_polyakov(T, mu, En(k, m), Phi, Phibar) + numberB(T, -mu, p0 - En(k, m)))
    term2 = abs(p0 + En(k, m)) < 1e-7 ? 0.0 : -spec(p0 + En(k, m), kq) * (fbar_polyakov(T, mu, En(k, m), Phi, Phibar) + numberB(T, -mu, p0 + En(k, m)))
    return k^2 * m * (term1 + term2) / (8π^2 * En(k, m)) ## -1 factor..check calc in goodnotes
end

function _imagpart_baryons_spectral_normal_q0_s_cor_q_rdf(spec, T, mu, q, p0, m, param::Parameters, Phi, Phibar; tol=1e-5, maxevals=1e7)
    return integrate(x -> _integ_s_cor_q_rdf(spec, T, mu, x[1], x[2], q, p0, m, param, Phi, Phibar), [0.0, -1.0], [param.Λ, 1.0])
end

function imagpart_baryons_spectral_normal_q0_s_bound_q_rdf(T, mu, p0, q, m, param, Phi, Phibar)
    md = find_diquark_energy_q_rdf(T, mu, q, m, param, Phi, Phibar)
    if md == 0.0
        return 0.0
    end
    Zk = wave_function_renormalization_diquark_rdf(T, mu, q, m, md, param, Phi, Phibar)
    ed = 0.5 * abs((p0^2 - m^2 + md^2) / p0)

    return imagpart_baryon_q_rdf(T, mu, p0, q, m, md, param, Phi, Phibar) * Zk / (2 * ed)
end

function _get_interpolated_imagpart_baryon_q_rdf(spec, T, mu, q, m, param, Phi, Phibar)
    threshold_cor = sqrt(q^2 + 9m^2) - 3 * mu
    cutoff_cor = sqrt(param.Λ^2 + m^2) + 2 * sqrt(param.Λ^2 + q^2 + m^2) - 3 * mu

    imag_cor = chebinterp(x -> _imagpart_baryons_spectral_normal_q0_s_cor_q_rdf(spec, T, mu, q, x, m, param, Phi, Phibar), 500, threshold_cor, cutoff_cor)
    function imag_cor_func(omega)
        if omega < 0
            return -imag_cor_func(-omega)
        end
        if threshold_cor < omega < cutoff_cor
            return imag_cor(omega)
        end
        return 0.0
    end
    md = find_diquark_energy_q_rdf(T, mu, q, m, param, Phi, Phibar)

    threshold_bound_1, threshold_bound_2 = sqrt(q^2 + (m + md)^2) - 2 * mu, sqrt(q^2 + (m - md)^2) - 2 * mu
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + q^2 + md^2) - 2 * mu

    imag_bound_1 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_q_rdf(T, mu, x, q, m, param, Phi, Phibar), 500, threshold_bound_1, cutoff_bound)
    imag_bound_2 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_q_rdf(T, mu, x, q, m, param, Phi, Phibar), 500, 0.0, threshold_bound_2)

    function imag_bound_func(omega)
        if omega < 0
            return -imag_bound_func(-omega)
        end
        if 0.0 < omega < threshold_bound_2
            return imag_bound_2(omega)
        end
        if threshold_bound_1 < omega < cutoff_bound
            return imag_bound_1(omega)
        end
        return 0.0
    end

    return x -> imag_cor_func(x) + imag_bound_func(x)
end

function _get_baryon_mass_spectral(T, mu, m, param, Phi, Phibar)
    spec = _get_interpolated_spectral_rdf(T, mu, m, param, Phi, Phibar)
    imag = _get_interpolated_imagpart_baryon_q0_rdf(spec, T, mu, m, param, Phi, Phibar)
    md = find_diquark_energy_q_rdf(T, mu, 0.0, m, param, Phi, Phibar)
    p00 = p00_baryons_interp_rdf(spec, T, mu, m, md, param, Phi, Phibar)
    realpart(om) = _get_realpart_baryon_spectral_q(imag, imag, om, m, p00, param)
    if realpart(0.0) * realpart(1.0) < 0.0
        return bisection(realpart, 0.0, 1.0)
    end
    return 0.0
end

function show_realpart_baryon_spectral_q0(spec, imag, T, mu, m, param, Phi, Phibar)
    md = find_diquark_energy_q_rdf(T, mu, 0.0, m, param, Phi, Phibar)
    p00 = p00_baryons_interp_rdf(spec, T, mu, m, md, param, Phi, Phibar)
    realpart(om) = _get_realpart_baryon_spectral_q(imag, imag, om, m, p00, param)
    return realpart
end

function distribution_baryon_q_spectral_rdf(T, mu, mq, q, param, Phi, Phibar)
    spec = _get_interpolated_spectral_rdf(T, mu, mq, param, Phi, Phibar)
    imag_part_q = _get_interpolated_imagpart_baryon_q_rdf(spec, T, mu, q, mq, param, Phi, Phibar)
    imagpart0 = _get_interpolated_imagpart_baryon_q0_rdf(spec, T, mu, mq, param, Phi, Phibar)
    p00 = p00_baryons_interp(spec, T, mu, mq, mD, param)
    phase_shift(om) = _get_phase_shift_baryon_spectral_q(imag_part_q, imagpart0, om, mq, p00, param)
    stat_factor(om)::Float64 = FD_dist(T, 3mu, om) * (1 - FD_dist(T, 3mu, om)) + FD_dist(T, -3mu, om) * (1 - FD_dist(T, -3mu, om))

    integrand(om)::Float64 = stat_factor(om) * phase_shift(om)

    return integrate(integrand, 0.0, 2.0)::Float64
end

function distribution_baryon_q_spectral_rdf(spec,T, mu, mq, q, p00,param, Phi, Phibar)
    imag_part_q = _get_interpolated_imagpart_baryon_q_rdf(spec, T, mu, q, mq, param, Phi, Phibar)
    imagpart0 = _get_interpolated_imagpart_baryon_q0_rdf(spec, T, mu, mq, param, Phi, Phibar)
    phase_shift(om) = _get_phase_shift_baryon_spectral_q(imag_part_q, imagpart0, om, mq, p00, param)
    stat_factor(om)::Float64 = FD_dist(T, 3mu, om) * (1 - FD_dist(T, 3mu, om)) + FD_dist(T, -3mu, om) * (1 - FD_dist(T, -3mu, om))

    integrand(om)::Float64 = stat_factor(om) * phase_shift(om)

    return integrate(integrand, 0.0, 2.0)::Float64
end
