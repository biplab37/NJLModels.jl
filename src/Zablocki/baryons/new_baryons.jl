# calculating baryon with full diquark spectral function. Once this works should rename the file

function integ_s_cor(T, mu, k, p0, m, param)
    term1 = abs(p0 - En(k, m))<1e-5 ? 0.0 : spectral_function_D_normal(T, mu, p0 - En(k, m), k, m, param) * (1 - numberF(T, mu, En(k, m)) + numberB(T, mu, p0 - En(k, m)))
    term2 = abs(p0 + En(k, m))<1e-5 ? 0.0 : -spectral_function_D_normal(T, mu, p0 + En(k, m), k, m, param) * (numberF(T, mu, En(k, m)) + numberB(T, mu, p0 + En(k, m)))
    return k^2 * m * (term1 + term2) / (4π^2 * En(k, m)) ## -1 factor..check calc in goodnotes
end

function imagpart_baryons_spectral_normal_q0_s_cor(T, mu, p0, m, param::Parameters; tol=1e-3, maxevals=1e5)
    return quadgk(k -> integ_s_cor(T, mu, k, p0, m, param), 0.0, param.Λ; atol=tol, rtol=tol, maxevals=maxevals)[1]
end

function _integ_sigma_B_00_s_cor(T, mu, k, nu, m, param)
    ek = En(k, m)
    stat_factor1 = (-1 + numberF(T, mu, ek) - numberB(T, mu, nu)) / (ek + nu)
    stat_factor2 = (-numberF(T, mu, ek) - numberB(T, mu, nu)) / (ek - nu)

    return k^2 * m * spectral_function_D_normal(T, mu, nu, k, m, param) * (stat_factor1 + stat_factor2) / (2π^2)
end

function sigma_B_00_s_cor(T, mu, m, param)
    cutoff = Float64(param.Λ)

    lower_bound = [0.0, 0.0]
    upper_bound = [cutoff, 6 * cutoff]

    return integrate(x -> _integ_sigma_B_00_s_cor(T, mu, x[1], x[2], m, param), lower_bound, upper_bound)
end

function _integ_sigma_B_00_s_bound(T, mu, k, m, param)
    ed = find_diquark_energy_q(T, mu, k, m, param)

    if ed == 0.0
        return 0.0
    end
    
    ek = En(k, m)
    Zk = wave_function_renormalization_diquark(T, mu, k, m, ed, param)

    stat_factor1 = -(1 - numberF(T, mu, ek) + numberB(T, mu, ed)) / (ek + ed)
    stat_factor2 = -(numberF(T, mu, ek) + numberB(T, mu, ed)) / (ek - ed)

    return k^2 * m * Zk * (stat_factor1 + stat_factor2) / (8π^2*ek)
end

function sigma_B_00_s_bound(T, mu, m, param)
    cutoff = Float64(param.Λ)
    return integrate(k -> _integ_sigma_B_00_s_bound(T, mu, k, m, param), 0.0, cutoff)
end

function realpart_baryons_spectral_normal_q0_s(T, mu, p0, m, param)
    return realpart_kramers_kronig(x -> imagpart_baryons_spectral_normal_q0_s_cor(T, mu, x, m, param), p0, 6 * param.Λ)
end

export imagpart_baryons_spectral_normal_q0_s, realpart_baryons_spectral_normal_q0_s
