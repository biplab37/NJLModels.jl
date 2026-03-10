## This file contains code related to computing spectral functions

function spectral_function_meson_q(T, mu, ω, q, Nm, param::Parameters)
    gap = massgap(T, mu, param).zero

    return spectral_function_meson_q(T, mu, ω, q, gap, Nm, param)
end

function spectral_function_meson_q(T, mu, ω, q, gap, Nm, param::Parameters)
    impart = imagpart_meson_q(ω, T, mu, q, gap, Nm, param)
    repart = realpart_meson_q(ω, T, mu, q, gap, Nm, param)

    return impart / (π * (impart^2 + repart^2))
end

spectral_function_pi_q(T, mu, ω, q, gap, param::Parameters) = spectral_function_meson_q(T, mu, ω, q, gap, 0, param)
spectral_function_sigma_q(T, mu, ω, q, gap, param::Parameters) = spectral_function_meson_q(T, mu, ω, q, gap, 1, param)

function _has_bound_state(rep::Function, threshold)
    if rep(0.0)*rep(threshold) > 0
        return 0.0
    end
    return UsefulFunctions.bisection(rep, 0.0, threshold)
end

function _has_meson_bound_state_q0(T, mu, gap, Nm, param)
    rep(ω) = realpart_meson_q(ω, T, mu, 0.0, gap, Nm, param)
    return _has_bound_state(rep, 2*gap[1])
end

function _has_meson_bound_state(T, mu, q, gap, Nm, param)
    rep(ω) = realpart_meson_q(ω, T, mu, q, gap, Nm, param)
    return _has_bound_state(rep, 2*gap[1])
end

function _has_diquark_bound_state_q0(T, mu, m, param)
    rep(ω) = realpart_D_normal(T, mu, ω, 0.0, m, param)
    return _has_bound_state(rep, 2*m)
end

function _has_diquark_bound_state(T, mu, q, m, param)
    rep(ω) = realpart_D_normal(T, mu, ω, q, m, param)
    return _has_bound_state(rep, 2*sqrt(m^2 + q^2))
end

function _wave_function_renormalization_meson(T, mu, q, gap, Nm, param)
    M = _has_meson_bound_state(T, mu, q, gap, Nm, param)
    if M == 0.0
        return 0.0
    end

    rep(ω) = realpart_meson_q(ω, T, mu, q, gap, Nm, param)
    return 1 / abs(UsefulFunctions._derivative(rep, M))
end

function _f_sum_meson(T, mu, q, gap, Nm, param)
    M = _has_meson_bound_state(T, mu, q, gap, Nm, param)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2*M*_wave_function_renormalization_meson(T, mu, q, gap, Nm, param)
    end

    _fsum += integrate(o->2*o*spectral_function_meson_q(T, mu, o, q, gap, Nm, param), 0.0, 4*param.Λ)
    return _fsum
end

function _f_sum_diquark(T, mu, q, m, param)
    M = _has_diquark_bound_state(T, mu, q, m, param)
    _fsum = 0.0

    if M > 0.0
        _fsum += 2*M*wave_function_renormalization_diquark(T, mu, q, m, M, param)
    end

    _fsum += integrate(o->2*o*spectral_function_D_normal(T, mu, o, q, m, param), 0.0, 4*param.Λ)
    return _fsum
end
