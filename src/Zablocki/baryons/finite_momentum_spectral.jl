## creating baryon with spectral function with finite momentum.

function _integ_s_cor_q(spec, T, mu, k, u, q, p0, m, param)
    kq = sqrt(k^2 + q^2 + 2 * k * q * u)
    term1 = abs(p0 - En(k, m)) < 1e-7 ? 0.0 : spec(p0 - En(k, m), kq) * (1 - FD_dist(T, mu, En(k, m)) + numberB(T, mu, p0 - En(k, m)))
    term2 = abs(p0 + En(k, m)) < 1e-7 ? 0.0 : -spec(p0 + En(k, m), kq) * (FD_dist(T, mu, En(k, m)) + numberB(T, mu, p0 + En(k, m)))
    return k^2 * m * (term1 + term2) / (8π^2 * En(k, m)) ## -1 factor..check calc in goodnotes
end

function _imagpart_baryons_spectral_normal_q0_s_cor_q(spec, T, mu, q, p0, m, param::Parameters; tol=1e-5, maxevals=1e7)
    return integrate(x -> _integ_s_cor_q(spec, T, mu, x[1], x[2], q, p0, m, param), [0.0, -1.0], [param.Λ, 1.0])
end

function imagpart_baryons_spectral_normal_q0_s_bound_q(T, mu, p0, q, m, param)
    md = find_diquark_energy_q(T, mu, q, m, param)
    if md == 0.0
        return 0.0
    end
    Zk = wave_function_renormalization_diquark(T, mu, q, m, md, param)
    ed = 0.5 * abs((p0^2 - m^2 + md^2) / p0)

    return imagpart_baryon_q(T, mu, p0, q, m, md, param) * Zk / (2 * ed)
end

function _get_interpolated_imagpart_baryon_q(spec, T, mu, q, m, param)
    threshold_cor = sqrt(q^2 + 9m^2) - 3 * mu
    cutoff_cor = sqrt(param.Λ^2 + m^2) + 2 * sqrt(param.Λ^2 + q^2 + m^2) - 3 * mu

    imag_cor = chebinterp(x -> _imagpart_baryons_spectral_normal_q0_s_cor_q(spec, T, mu, q, x, m, param), 500, threshold_cor, cutoff_cor)
    function imag_cor_func(omega)
        if omega < 0
            return -imag_cor_func(-omega)
        end
        if threshold_cor < omega < cutoff_cor
            return imag_cor(omega)
        end
        return 0.0
    end
    md = find_diquark_energy_q(T, mu, q, m, param)

    threshold_bound_1, threshold_bound_2 = sqrt(q^2 + (m + md)^2) - 2 * mu, sqrt(q^2 + (m - md)^2) - 2 * mu
    cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + q^2 + md^2) - 2 * mu

    imag_bound_1 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_q(T, mu, x, q, m, param), 500, threshold_bound_1, cutoff_bound)
    imag_bound_2 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound_q(T, mu, x, q, m, param), 500, 0.0, threshold_bound_2)

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

function _get_realpart_baryon_spectral_q(imag_part_q, imagpart0, p0, m, p00, param)
    return p00 - realpart_kramers_kronig_q(imag_part_q, imagpart0, p0, 3 * sqrt(param.Λ^2 + m^2))
end

function _get_phase_shift_baryon_spectral_q(imag_part_q, imagpart0, p0, m, p00, param)
    real_part = _get_realpart_baryon_spectral_q(imag_part_q, imagpart0, p0, m, p00, param)
    imag_part = imag_part_q(p0)

    return atan(imag_part, real_part)
end

function distribution_baryon_q_spectral(T, mu, q, param)
    mq::Float64, ome::Float64 = massgap(T, mu, param).zero
    mD::Float64 = find_diquark_energy_q(T, mu, q, mq, param)
    mus::Float64 = mu + ome

    spec = _get_interpolated_spectral(T, mu, mq, param)
    imag_part_q = _get_interpolated_imagpart_baryon_q(spec, T, mu, q, mq, param)
    imagpart0 = _get_interpolated_imagpart_baryon_q0(spec, T, mu, mq, param)
    p00 = p00_baryons_interp(spec, T, mu, mq, mD, param)
    phase_shift(om) = _get_phase_shift_baryon_spectral_q(imag_part_q, imagpart0, om, mq, p00, param)
    stat_factor(om)::Float64 = FD_dist(T, 3mu, om) * (1 - FD_dist(T, 3mu, om)) + FD_dist(T, -3mu, om) * (1 - FD_dist(T, -3mu, om))

    integrand(om)::Float64 = stat_factor(om) * phase_shift(om)

    return integrate(integrand, 0.0, 2.0)::Float64
end