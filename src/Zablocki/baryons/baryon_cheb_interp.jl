## first make the interpolated version of the spectral function of the diquark for speeding up the calculations

function _get_threshold_cutoff(q, m, mu, Λ)
    return  sqrt(4 * m^2 + q^2) - 2 * mu, sqrt(4 * Λ^2 + 4 * m^2 + q^2) - 2 * mu
end

function _spectral_diquark_interpolated_k(T, mu, q, m, param)
    func(x) = spectral_function_D_normal(T, mu, x, q, m, param)
    threshold, cutoff = _get_threshold_cutoff(q, m, mu, param.Λ)

    landau_piece = (q > 2mu) ? chebinterp(func, 200, 0, q - 2mu) : 0.0
    qp_piece = chebinterp(func, 500, threshold, cutoff)

    interpolated_spec(x) = (x + 2mu < q) ? landau_piece(x) : ((threshold < x < cutoff) ? qp_piece(x) : 0.0)

    return interpolated_spec
end

function _integ_s_cor(spec, T, mu, k, p0, m, param)
    term1 = abs(p0 - En(k, m)) < 1e-7 ? 0.0 : spec(p0 - En(k, m), k) * (1 - FD_dist(T, mu, En(k, m)) + numberB(T, mu, p0 - En(k, m)))
    term2 = abs(p0 + En(k, m)) < 1e-7 ? 0.0 : -spec(p0 + En(k, m), k) * (FD_dist(T, mu, En(k, m)) + numberB(T, mu, p0 + En(k, m)))
    return k^2 * m * (term1 + term2) / (4π^2 * En(k, m)) ## -1 factor..check calc in goodnotes
end

function _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, p0, m, param::Parameters; tol=1e-5, maxevals=1e7)
    return quadgk(k -> _integ_s_cor(spec, T, mu, k, p0, m, param), 0.0, param.Λ; atol=tol, rtol=tol, maxevals=maxevals)[1]
end

function _spectral_diquark_interpolated(T, mu, m, param)
    func(x) = spectral_function_D_normal(T, mu, x[1], x[2], m, param)
    _, cutoff = _get_threshold_cutoff(param.Λ, m, mu, param.Λ)

    interpolated_spec = chebinterp(func, (800, 100), [0.0, 0.0], [cutoff, param.Λ])

    return (om, q) -> interpolated_spec([om, q])
end

function _get_interpolated_spectral(T, mu, m, param)
    order = (500, 100)
    lb_qp = [0.0, 0.0]
    lb_ld = [0.0, 2 * mu]
    ub = [1.0, 2 * param.Λ + 0.1]

    # Landau Damping: Maps [u, q] back to 0 < ω < q - 2μ
    interp_ld = chebinterp(order, lb_ld, ub) do x
        u, q = x[1], x[2]
        ω_upper = q - 2 * mu
        return spectral_function_D_normal(T, mu, ω_upper * u, q, m, param)
    end
    # Quasiparticle: Maps [t, q] back to ω_lower < ω < ω_upper
    interp_qp = chebinterp(order, lb_qp, ub) do x
        t, q = x[1], x[2]
        ω_lower, ω_upper = _get_threshold_cutoff(q, m, mu, param.Λ)

        ω = ω_lower + t * (ω_upper - ω_lower)
        return spectral_function_D_normal(T, mu, ω, q, m, param)
    end

    function spec(ω, q)
        if ω < 0
            return -spec(-ω, q)
        end
        if q > 2 * mu && 0.0 < ω < q - 2 * mu
            u = ω / (q - 2 * mu)
            return interp_ld([u, q])
        end

        ω_lower_qp, ω_upper_qp = _get_threshold_cutoff(q, m, mu, param.Λ)
        if ω_lower_qp < ω < ω_upper_qp
            t = (ω - ω_lower_qp) / (ω_upper_qp - ω_lower_qp)
            return interp_qp([t, q])
        end

        return 0.0
    end
    return spec
end
function _integ_sigma_B_00_s_cor_interp(spec, T, mu, k, nu, m, param)
    ek = En(k, m)
    stat_factor1 = (-1 + FD_dist(T, mu, ek) - numberB(T, mu, nu)) / (ek + nu)
    stat_factor2 = (-FD_dist(T, mu, ek) - numberB(T, mu, nu)) / (ek - nu)

    return k^2 * m * spec(nu, k) * (stat_factor1 + stat_factor2) / (2π^2)
end

function sigma_B_00_s_cor_interp(spec, T, mu, m, param)
    cutoff = Float64(param.Λ)

    lower_bound = [0.0, 0.0]
    upper_bound = [cutoff, 6 * cutoff]

    return integrate(x -> _integ_sigma_B_00_s_cor_interp(spec, T, mu, x[1], x[2], m, param), lower_bound, upper_bound)
end

function _imagpart_baryons_spectral_normal_total_q0_interp(spec, T, mu, p0, m, param)
    return _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, p0, m, param) + imagpart_baryons_spectral_normal_q0_s_bound(T, mu, p0, m, param)
end

function _realpart_baryons_spectral_normal_q0_s_interp(spec, T, mu, p0, m, param)
    return realpart_kramers_kronig(x -> (_imagpart_baryons_spectral_normal_total_q0_interp(spec, T, mu, x, m, param)), p0, 3 * sqrt(param.Λ^2 + m^2))
end

function p00_baryons_interp(spec, T, mu, m, mD, param)
    1 / coupling_B(T, mu, m, mD, param) - sigma_B_00_s_bound(T, mu, m, param) - sigma_B_00_s_cor_interp(spec, T, mu, m, param)
end
function p00_baryons(T, mu, m, mD, param)
    1 / coupling_B(T, mu, m, mD, param) - sigma_B_00_s_bound(T, mu, m, param) - sigma_B_00_s_cor(T, mu, m, param)
end

function _get_imag_baryon_total_q0_interp(T, mu, m, param)
    func_cor(x) = imagpart_baryons_spectral_normal_q0_s_cor(T, mu, x, m, param)

    threshold_cor = 3 * (m - mu)
    cutoff_cor = 3 * (sqrt(param.Λ^2 + m^2) - mu)

    interp_cor = chebinterp(func_cor, 500, threshold_cor, cutoff_cor)

    cor_interp(x) = (threshold_cor < x < cutoff_cor) ? interp_cor(x) : 0.0

    func_bound(x) = imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param)

    mD = find_mass_D(T, mu, param)[1]
    threshold_bound = m + mD - 3 * mu
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + mD^2) - 3 * mu

    interp_bound = chebinterp(func_bound, 1000, threshold_bound, cutoff_bound)

    bound_interp(x) = (threshold_bound < x < cutoff_bound) ? interp_bound(x) : 0.0

    return x -> cor_interp(x) + bound_interp(x)
end

function _get_interpolated_imagpart_baryon_q0(spec, T, mu, m, param)
    threshold_cor = 3m
    cutoff_cor = 3 * sqrt(param.Λ^2 + m^2)

    imag_cor = chebinterp(x -> _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, x, m, param), 500, threshold_cor, cutoff_cor)
    function imag_cor_func(omega)
        if omega < 0
            return -imag_cor_func(-omega)
        end
        if threshold_cor < omega < cutoff_cor
            return imag_cor(omega)
        end
        return 0.0
    end
    md = find_diquark_energy_q(T, mu, 0.0, m, param)

    threshold_bound_1, threshold_bound_2 = m + md, abs(m - md)
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + md^2)

    imag_bound_1 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param), 500, threshold_bound_1, cutoff_bound)
    imag_bound_2 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param), 500, 0.0, threshold_bound_2)

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

function _get_realpart_baryon_spectral(imag_part, p0, m, p00, param)
    return p00 - realpart_kramers_kronig(imag_part, p0, 3 * sqrt(param.Λ^2 + m^2))
end

function _get_phase_shift_baryon_spectral(imag_part, p0, m, p00, param)
    real_part = _get_realpart_baryon_spectral(imag_part, p0, m, p00, param)
    imag_part = imag_part(p0)

    return atan(imag_part, real_part)
end

