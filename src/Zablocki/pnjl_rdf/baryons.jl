function kallen(a, b, c)
    return a^2 + b^2 + c^2 - 2 * (a * b + b * c + c * a)
end
function imagpart_baryon_q0_C_rdf(T, mu, ω, mq, mD, param, Phi, Phibar)
    z = ω + 3 * mu
    if abs(z) < 1e-5
        return 0.0
    end
    func(om) = 2 * mq * sqrt(kallen(om^2, mD^2, mq^2)) * (numberB(T, -mu, (om^2 + mD^2 - mq^2) / (2 * om)) + f_polyakov(T, mu, -(om^2 - mD^2 + mq^2) / (2 * om), Phi, Phibar)) / (4π * om^2)

    if imag(z) != 0.0
        return func(z)
    end

    if z==0 || (mD - mq)^2 <= z^2 <= (mq + mD)^2 || z^2 >= (sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2))^2
        return 0.0
    end

    return func(z)
end

function imagpart_baryon_q0_rdf(T, mu, ω, mq, mD, param, Phi, Phibar)
    if ω == 0
        return 0.0
    end
    ω = ω + 3 * mu
    eq = 0.5 * abs((ω^2 + mq^2 - mD^2) / ω)
    eD = 0.5 * abs((ω^2 - mq^2 + mD^2) / ω)
    # if eD < 1e-5
    #     @show mq, mD, ω
    # end
    p = sqrt(abs(mq^4 + (mD^2 - ω^2)^2 - 2mq^2 * (mD^2 + ω^2))) / (2 * abs(ω))
    factor = -p * mq / (π)
    term1 = (1 - f_polyakov(T, mu, eq, Phi, Phibar) + numberB(T, mu, eD)) / (eq + eD)
    term2 = -(1 - fbar_polyakov(T, mu, eq, Phi, Phibar) + numberB(T, -mu, eD)) / (eq + eD)
    term3 = (fbar_polyakov(T, mu, eq, Phi, Phibar) + numberB(T, mu, eD)) * PrincipalValue(eD - eq)
    term4 = -(f_polyakov(T, mu, eq, Phi, Phibar) + numberB(T, -mu, eD)) * PrincipalValue(eD - eq)
    if abs(ω) > sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2)
        return 0.0
    end
    if ω > 0
        if ω > (mq + mD)
            return factor * term2
        elseif ω < (mD - mq)
            return factor * term4
        else
            return 0.0
        end
    else
        if ω < -(mq + mD)
            return factor * term1
        elseif ω > -(mD - mq)
            return factor * term3
        else
            return 0.0
        end
    end
end
function coupling_B_rdf(T, mu, mq, mD, param, Phi, Phibar)
    pol_D(ω) = realpart_D_normal_q0_rdf(T, mu, ω, mq, param, Phi, Phibar)
    der = UsefulFunctions._derivative(pol_D, mD)
    # return 16 * mD / (mq * abs(der))
    return  72.67151584450362
end


function Π0_B_rdf(T, mu, mq, mD, param, Phi, Phibar)
    factor = mq / (π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    integrand1(p) = p^2 * ((1 - f_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, mu, eD(p)))*PrincipalValue(3*mu + eq(p) + eD(p)) + (1 - fbar_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, -mu, eD(p)))*PrincipalValue(-3*mu + eq(p) + eD(p))) / (eq(p) * eD(p))
    integrand2(p) = p^2 * ((fbar_polyakov(T, -mu, eq(p), Phi, Phibar) + numberB(T, mu, eD(p)))*PrincipalValue(3*mu - eq(p) + eD(p)) + (f_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, -mu, eD(p)))*PrincipalValue(-3*mu - eq(p) + eD(p))) / (eq(p) * eD(p))

    return 1 / (coupling_B_rdf(T, mu, mq, mD, param, Phi, Phibar)) - factor * integrate(p -> integrand1(p) + integrand2(p), 0.0, param.Λ)
end

function realpart_baryon_q0_rdf(T, mu, ω, mq, mD, param, Phi, Phibar)
    factor = mq / (π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    term1(p) = ( - f_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, mu, eD(p))) * PrincipalValue(ω + 3*mu+ eq(p) + eD(p))
    term2(p) = -( - fbar_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, -mu, eD(p))) * PrincipalValue(ω + 3*mu - eq(p) - eD(p))
    term3(p) = (fbar_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, mu, eD(p))) * PrincipalValue(ω + 3*mu - eq(p) + eD(p))
    term4(p) = -(f_polyakov(T, mu, eq(p), Phi, Phibar) + numberB(T, -mu, eD(p))) * PrincipalValue(ω + 3*mu + eq(p) - eD(p))

    integrand_vac(p) = p^2*(PrincipalValue(ω + 3*mu + eq(p) + eD(p)) - PrincipalValue(ω + 3*mu - eq(p) - eD(p)) )/(eq(p) * eD(p))
    integrand(p) = p^2 * (term1(p) + term2(p) + term3(p) + term4(p)) / (eq(p) * eD(p))

    return 1 / coupling_B_rdf(T, mu, mq, mD, param, Phi, Phibar) - factor * integrate(integrand_vac, 0.0, param.Λ) - factor * integrate(integrand, 0.0, 2.0)
end

function realpart_baryon_q_rdf(T, mu, ω, q, mq, mD, param, Phi, Phibar)
    impart(x) = imagpart_baryon_q_rdf(T, mu, x, q, mq, mD, param, Phi, Phibar)

    cutoff = max(sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2)), sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2))
    repart_dependent = realpart_kramers_kronig_1(impart, ω, -cutoff-3mu, cutoff-3mu)

    return Π0_B_rdf(T, mu, mq, mD, param, Phi, Phibar) - repart_dependent
end

function Πq0_baryon_analytic_rdf(T, mu, ω, m, mD, param, Phi, Phibar)
    if imag(ω) == 0
        return realpart_baryon_q0_rdf(T, mu, real(ω), m, mD, param, Phi, Phibar) - 1im * imagpart_baryon_q0_C_rdf(T, mu, real(ω), m, mD, param, Phi, Phibar)
    end
    repart = realpart_baryon_q0_rdf(T, mu, ω, m, mD, param, Phi, Phibar)
    impart = (imag(ω) >= 0.0) ? 0.0 : 1im * imagpart_baryon_q0_C_rdf(T, mu, ω, m, mD, param, Phi, Phibar)

    return repart - 2impart
end

function find_baryon_mass_q0_rdf(T, mu, mq, mD, param, Phi, Phibar, guess=[0.98, 0.0])
    rep(z) = realpart_baryon_q0_rdf(T, mu, z, mq, mD, param, Phi, Phibar)
    if mq + mD > 3*mu && rep(0.0) * rep(mq + mD) < 0.0
        return [bisection(rep, 0.0, (mq + mD)), 0.0]
    end
    function ff!(F, x)
        term = Πq0_baryon_analytic_rdf(T, mu, x[1] - 1im * x[2] / 2, mq, mD, param, Phi, Phibar)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return mcpsolve(ff!, [0.0, 0.0], [2.0, 2.0], guess).zero
    # return [0.0, 0.0]

end

function phase_shift_baryon_q_rdf(T, mus, ω, q, mq, mD, param, Phi, Phibar)
    impart = imagpart_baryon_q_rdf(T, mus, ω, q, mq, mD, param, Phi, Phibar)
    repart = realpart_baryon_q_rdf(T, mus, ω, q, mq, mD, param, Phi, Phibar)

    if abs(impart) <= 1e-4 && repart < 0.0
        if ω<-3mus
            return -π
        else
            return π
        end
    end
    return atan(impart, repart)
end

function distribution_baryon_q_rdf(T, mu, q, mq, mD, param, Phi, Phibar)
    stat_factor(om)::Float64 = FD_dist(T, 3mu, om) * (1 - FD_dist(T, 3mu, om)) + FD_dist(T, -3mu, om) * (1 - FD_dist(T, -3mu, om))

    integrand(om)::Float64 = stat_factor(om) * phase_shift_baryon_q_rdf(T, mu, om, q, mq, mD, param, Phi, Phibar)
    # return phase_shift_baryon_q_rdf(T, mu, 3*mu, q, mq, mD, param, Phi, Phibar) + phase_shift_baryon_q_rdf(T, mu, -3*mu, q, mq, mD, param, Phi, Phibar) 
    return integrate(integrand, 3mu-0.1, 3mu+0.1)::Float64
end

# function _integ_s_cor(spec, T, mu, k, p0, m, param, Phi, Phibar)
#     term1 = abs(p0 - En(k, m)) < 1e-7 ? 0.0 : spec(p0 - En(k, m), k) * (1 - FD_dist(T, mu, En(k, m), Phi, Phibar) + numberB(T, mu, p0 - En(k, m)))
#     term2 = abs(p0 + En(k, m)) < 1e-7 ? 0.0 : -spec(p0 + En(k, m), k) * (FD_dist(T, mu, En(k, m), Phi, Phibar) + numberB(T, mu, p0 + En(k, m)))
#     return k^2 * m * (term1 + term2) / (4π^2 * En(k, m)) ## -1 factor..check calc in goodnotes
# end

# function _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, p0, m, param::Parameters, Phi, Phibar; tol=1e-5, maxevals=1e7)
#     return quadgk(k -> _integ_s_cor(spec, T, mu, k, p0, m, param, Phi, Phibar), 0.0, param.Λ; atol=tol, rtol=tol, maxevals=maxevals)[1]
# end

# function _get_interpolated_spectral(T, mu, m, param, Phi, Phibar)
#     order = (500, 100)
#     lb_qp = [0.0, 0.0]
#     lb_ld = [0.0, 2 * mu]
#     ub = [1.0, 2 * param.Λ + 0.1]

#     # Landau Damping: Maps [u, q] back to 0 < ω < q - 2μ
#     interp_ld = chebinterp(order, lb_ld, ub) do x
#         u, q = x[1], x[2]
#         ω_upper = q - 2 * mu
#         return spectral_function_D_normal(T, mu, ω_upper * u, q, m, param, Phi, Phibar)
#     end
#     # Quasiparticle: Maps [t, q] back to ω_lower < ω < ω_upper
#     interp_qp = chebinterp(order, lb_qp, ub) do x
#         t, q = x[1], x[2]
#         ω_lower, ω_upper = _get_threshold_cutoff(q, m, mu, param.Λ)

#         ω = ω_lower + t * (ω_upper - ω_lower)
#         return spectral_function_D_normal(T, mu, ω, q, m, param, Phi, Phibar)
#     end

#     function spec(ω, q)
#         if ω < 0
#             return -spec(-ω, q)
#         end
#         if q > 2 * mu && 0.0 < ω < q - 2 * mu
#             u = ω / (q - 2 * mu)
#             return interp_ld([u, q])
#         end

#         ω_lower_qp, ω_upper_qp = _get_threshold_cutoff(q, m, mu, param.Λ)
#         if ω_lower_qp < ω < ω_upper_qp
#             t = (ω - ω_lower_qp) / (ω_upper_qp - ω_lower_qp)
#             return interp_qp([t, q])
#         end

#         return 0.0
#     end
#     return spec
# end

# function _integ_sigma_B_00_s_cor_interp(spec, T, mu, k, nu, m, param, Phi, Phibar)
#     ek = En(k, m)
#     stat_factor1 = (-1 + FD_dist(T, mu, ek, Phi, Phibar) - numberB(T, mu, nu)) / (ek + nu)
#     stat_factor2 = (-FD_dist(T, mu, ek, Phi, Phibar) - numberB(T, mu, nu)) / (ek - nu)

#     return k^2 * m * spec(nu, k) * (stat_factor1 + stat_factor2) / (2π^2)
# end

# function sigma_B_00_s_cor_interp(spec, T, mu, m, param, Phi, Phibar)
#     cutoff = Float64(param.Λ)

#     lower_bound = [0.0, 0.0]
#     upper_bound = [cutoff, 6 * cutoff]

#     return integrate(x -> _integ_sigma_B_00_s_cor_interp(spec, T, mu, x[1], x[2], m, param, Phi, Phibar), lower_bound, upper_bound)
# end

# function find_diquark_energy_q(T, mu, q, m, param, Phi, Phibar)
#     EQ = sqrt(q^2 + 4m^2)

#     if EQ > mu && realpart_D_normal(T, mu, 0.0, q, m, param, Phi, Phibar) * realpart_D_normal(T, mu, EQ - 2mu, q, m, param, Phi, Phibar) < 0.0
#         return bisection(x -> realpart_D_normal(T, mu, x, q, m, param, Phi, Phibar), 0.0, EQ - 2mu)
#     end

#     return 0.0
# end

# function imagpart_baryon_q(T, mu, om, q, mq, mD, param, Phi, Phibar)
#     ω = om + 3mu
#     if q == 0
#         return imagpart_baryon_q0(T, mu, ω, mq, mD, param, Phi, Phibar)
#     end
#     expb(eps) = exp(eps / (2T))

#     term1(x) = -2 * T * log((expb(x) * (1 - expb(-(-ω - mu + x)))) / (1 + expb(-(-ω - mu - x))))
#     term2(x) = 2 * T * log((expb(x) * (1 - expb(-(ω + mu + x)))) / (1 + expb(-(ω + mu - x))))

#     term3(x) = 2 * T * log((1 - expb(-(-ω - mu + x))) / (1 + expb(-(ω + mu + x))))
#     term4(x) = -2 * T * log((1 - expb(-(ω + mu + x))) / (1 + expb(-(-ω - mu + x))))

#     cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
#     factor = 2 * mq / (8 * π * q)

#     s = ω^2 - q^2

#     if abs(ω) >= cutoff || s == 0
#         return 0.0
#     end

#     if mD > mq
#         m1 = mD
#         m2 = mq
#     else
#         m1 = mq
#         m2 = mD
#     end

#     if s > 0
#         if (m1 - m2)^2 <= s <= (m1 + m2)^2
#             return 0.0
#         end
#         x1 = x_pm(abs(ω), q, m1, m2, -1)
#         x2 = x_pm(abs(ω), q, m1, m2, +1)
#         if s > (m1 + m2)^2
#             if ω >= 0
#                 return factor * (term2(x2) - term2(x1))
#             else
#                 return factor * (term1(x2) - term1(x1))
#             end
#         else
#             if ω >= 0
#                 return factor * (term4(x2) - term4(x1))
#             else
#                 return factor * (term3(x2) - term3(x1))
#             end
#         end
#     else
#         x2 = x_pm(ω, q, m1, m2, +1)
#         return factor * (term3(x2) + term4(x2))
#     end
# end


# function imagpart_baryons_spectral_normal_q0_s_bound(T, mu, p0, m, param, Phi, Phibar)
#     md = find_diquark_energy_q(T, mu, 0.0, m, param, Phi, Phibar)
#     if md == 0.0
#         return 0.0
#     end
#     Zk = wave_function_renormalization_diquark(T, mu, 0.0, m, md, param, Phi, Phibar)
#     ed = 0.5 * abs((p0^2 - m^2 + md^2) / p0)

#     return imagpart_baryon_q0(T, mu, p0, m, md, param, Phi, Phibar) * Zk / (2 * ed)
# end

# function _imagpart_baryons_spectral_normal_total_q0_interp(spec, T, mu, p0, m, param, Phi, Phibar)
#     return _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, p0, m, param, Phi, Phibar) + imagpart_baryons_spectral_normal_q0_s_bound(T, mu, p0, m, param, Phi, Phibar)
# end

# function _realpart_baryons_spectral_normal_q0_s_interp(spec, T, mu, p0, m, param, Phi, Phibar)
#     return realpart_kramers_kronig(x -> (_imagpart_baryons_spectral_normal_total_q0_interp(spec, T, mu, x, m, param, Phi, Phibar)), p0, 3 * sqrt(param.Λ^2 + m^2))
# end

# function p00_baryons_interp(spec, T, mu, m, mD, param, Phi, Phibar)
#     1 / coupling_B(T, mu, m, mD, param) - sigma_B_00_s_bound(T, mu, m, param) - sigma_B_00_s_cor_interp(spec, T, mu, m, param)
# end
# function p00_baryons(T, mu, m, mD, param)
#     1 / coupling_B(T, mu, m, mD, param) - sigma_B_00_s_bound(T, mu, m, param) - sigma_B_00_s_cor(T, mu, m, param)
# end

# function _get_imag_baryon_total_q0_interp(T, mu, m, param)
#     func_cor(x) = imagpart_baryons_spectral_normal_q0_s_cor(T, mu, x, m, param)

#     threshold_cor = 3 * (m - mu)
#     cutoff_cor = 3 * (sqrt(param.Λ^2 + m^2) - mu)

#     interp_cor = chebinterp(func_cor, 500, threshold_cor, cutoff_cor)

#     cor_interp(x) = (threshold_cor < x < cutoff_cor) ? interp_cor(x) : 0.0

#     func_bound(x) = imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param)

#     mD = find_mass_D(T, mu, param)[1]
#     threshold_bound = m + mD - 3 * mu
#     cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + mD^2) - 3 * mu

#     interp_bound = chebinterp(func_bound, 1000, threshold_bound, cutoff_bound)

#     bound_interp(x) = (threshold_bound < x < cutoff_bound) ? interp_bound(x) : 0.0

#     return x -> cor_interp(x) + bound_interp(x)
# end

# function _get_interpolated_imagpart_baryon_q0(spec, T, mu, m, param)
#     threshold_cor = 3m
#     cutoff_cor = 3 * sqrt(param.Λ^2 + m^2)

#     imag_cor = chebinterp(x -> _imagpart_baryons_spectral_normal_q0_s_cor(spec, T, mu, x, m, param), 500, threshold_cor, cutoff_cor)
#     function imag_cor_func(omega)
#         if omega < 0
#             return -imag_cor_func(-omega)
#         end
#         if threshold_cor < omega < cutoff_cor
#             return imag_cor(omega)
#         end
#         return 0.0
#     end
#     md = find_diquark_energy_q(T, mu, 0.0, m, param)

#     threshold_bound_1, threshold_bound_2 = m + md, abs(m - md)
#     cutoff_bound = sqrt(param.Λ^2 + m^2) + sqrt(param.Λ^2 + md^2)

#     imag_bound_1 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param), 500, threshold_bound_1, cutoff_bound)
#     imag_bound_2 = chebinterp(x -> imagpart_baryons_spectral_normal_q0_s_bound(T, mu, x, m, param), 500, 0.0, threshold_bound_2)

#     function imag_bound_func(omega)
#         if omega < 0
#             return -imag_bound_func(-omega)
#         end
#         if 0.0 < omega < threshold_bound_2
#             return imag_bound_2(omega)
#         end
#         if threshold_bound_1 < omega < cutoff_bound
#             return imag_bound_1(omega)
#         end
#         return 0.0
#     end

#     return x -> imag_cor_func(x) + imag_bound_func(x)
# end

# function _get_realpart_baryon_spectral(imag_part, p0, m, p00, param)
#     return p00 - realpart_kramers_kronig(imag_part, p0, 3 * sqrt(param.Λ^2 + m^2))
# end

# function _get_phase_shift_baryon_spectral(imag_part, p0, m, p00, param)
#     real_part = _get_realpart_baryon_spectral(imag_part, p0, m, p00, param)
#     imag_part = imag_part(p0)

#     return atan(imag_part, real_part)
# end