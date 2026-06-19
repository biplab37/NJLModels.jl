## This file contains code related to calculation of Baryons self energy with full diquark spectral function

function imagpart_baryon_q0_spectral(T, mu, ω, mq, param)
    spec1(x) = spectral_function_D_normal(T, mu, ω - x, sqrt(x^2 - mq^2), param)
    spec2(x) = spectral_function_D_normal(T, mu, ω + x, sqrt(x^2 - mq^2), param)

    integrand(x) = mq * sqrt(x^2 - mq^2) * ((1 - FD_dist(T, mu, x) + numberB(T, mu, ω - x)) * spec1(x) + (FD_dist(T, mu, x) + numberB(T, mu, ω + x)) * spec2(x))

    cutoff = sqrt(4 * param.Λ^2 + mq^2)

    return integrate(integrand, mq, cutoff)
end

function imagpart_diquark_interpol(T, mu, q, mq, param)
    cutoff = 2*sqrt(param.Λ^2 + mq^2)
    ωs = range(-cutoff, cutoff, length=200)
    imparts = [imagpart_D_normal(T, mu, ω, q, mq, param) for ω in ωs]
    f(x) = (abs(x) > cutoff) ? 0.0 : interp(imparts)(x/(2cutoff) + 0.5)
    return f
end

function spectral_function_D_normal_fast(T, mu, q, mq,  param)
    cutoff = 2*sqrt(param.Λ^2 + mq^2 + q^2/4)
    ωs = range(-cutoff, cutoff, length=200)
    imparts = [imagpart_D_normal(T, mu, ω, q, mq, param) for ω in ωs]
    impart_func = linear_interpolation(ωs, imparts)
    impart_interp(x) = (abs(x) >= cutoff) ? 0.0 : impart_func(x)
    ωs2 = 2*ωs
    pi0 = Π0_D(T, mu, mq, param)
    repart_vals = [pi0 - realpart_kramers_kronig(x->impart_interp(x), ω_val, cutoff) for ω_val in ωs2]
    repart_func = linear_interpolation(ωs2, repart_vals)
    repart_interp(x) = (abs(x) > 2*cutoff) ? pi0 : repart_func(x)

    return ω -> (abs(ω) > 2*cutoff) ? 0.0 : impart_interp(ω) / (π*(repart_interp(ω)^2 + impart_interp(ω)^2))
end

function spectral_function_D_normal_interpolated(T, mu, mq, param; orange=range(0.0, 2.0, length=200), qrange=range(0, 2.0, length=200))
    data_spectral = zeros(length(orange), length(qrange))

    for i in eachindex(qrange)
        func = spectral_function_D_normal_fast(T, mu, qrange[i], mq, param)
        data_spectral[:, i] = map(func, orange)
    end

    return linear_interpolation((orange, qrange), data_spectral)
end

function imagpart_baryon_q0_spectral_interpolated(T, mu, ω, mq, param)
    spec_func = spectral_function_D_normal_interpolated(T, mu, mq, param)

    spec1(x) = sign(ω -x) * spec_func(abs(ω - x), sqrt(x^2 - mq^2))
    spec2(x) = (ω + x >= 2.0) ? 0.0 : spec_func(ω + x, sqrt(x^2 - mq^2))

    integrand(x) = (abs(ω - x)<1e-3) ? (FD_dist(T, mu, x) + numberB(T, mu, ω + x)) * spec2(x) : mq * sqrt(x^2 - mq^2) * ((1 - FD_dist(T, mu, x) + numberB(T, mu, ω - x)) * spec1(x) + (FD_dist(T, mu, x) + numberB(T, mu, ω + x)) * spec2(x))

    cutoff = sqrt(4 * param.Λ^2 + mq^2)

    return integrate(integrand, mq, cutoff)
end

function realpart_omega0(T, mu, mq, param)
    integrand(ep, nu) = (nu==0) ? 0.0 : mq*sqrt(ep^2 - mq^2)*spectral_function_D_normal(T, mu, nu, sqrt(ep^2 - mq^2), param) * ((1 - FD_dist(T, mu, ep) + numberB(T, mu, nu))*PrincipalValue(nu + ep) + (FD_dist(T, mu, ep) + numberB(T, mu, nu))*PrincipalValue(nu - ep))/π^2
    cutoff_nu = 2 * sqrt(param.Λ^2 + mq^2)
    return integrate(x->integrand(x...), [mq, -cutoff_nu], [sqrt(4 * param.Λ^2 + mq^2), cutoff_nu])
end

function realpart_baryons_q0_spectral(T, mu, ω, mq, param)
    mD = Zablocki.find_mass_D(T, mu, param)[1]
    imaginary(x) = imagpart_baryon_q0_spectral(T, mu, x, mq, param)
    cutoff = 2 * sqrt(param.Λ^2 + mq^2)
    repart = realpart_kramers_kronig(imaginary, ω, cutoff)
    return 1/(coupling_B(T, mu, mq, mD, param)) - repart
end

function realpart_omega0_interpolated(T, mu, mq, param)
    spectral_func = spectral_function_D_normal_interpolated(T, mu, mq, param)
    integrand(ep, nu) = (nu==0) ? 0.0 : mq*sqrt(ep^2 - mq^2)*spectral_func(nu, sqrt(ep^2 - mq^2)) * ((1 - FD_dist(T, mu, ep) + numberB(T, mu, nu))*PrincipalValue(nu + ep) + (FD_dist(T, mu, ep) + numberB(T, mu, nu))*PrincipalValue(nu - ep))/π^2
    cutoff_nu = 2 * sqrt(param.Λ^2 + mq^2)
    return integrate(x->integrand(x...), [mq, 0.0], [sqrt(4 * param.Λ^2 + mq^2), cutoff_nu])
end

function realpart_baryon_q0_spectral_interpolated(T, mu, ω, mq, param)
    mD = Zablocki.find_mass_D(T, mu, param)[1]
    imaginary(x) = imagpart_baryon_q0_spectral_interpolated(T, mu, x, mq, param)
    cutoff = 2 * sqrt(param.Λ^2 + mq^2)
    repart = realpart_kramers_kronig(imaginary, ω, cutoff, maxevals=100)
    return 1/(coupling_B(T, mu, mq, mD, param)) - repart
end
