function imagpart_baryon_q1(T, mu, ω, q, mq, mD, param)
    if q == 0
        return imagpart_baryon_q0(T, mu, ω, mq, mD, param)
    end

    expb(eps) = exp(eps / T)

    term1(x) = 2 * T * log((expb(x / 2) + expb(mu - abs(ω) / 2)) / (1 - expb(mu - abs(ω) / 2 + x / 2)))

    term3(y) = 2 * T * log((1 - expb(y + 2 * mu - ω)) / (expb(y) - expb(2 * mu - ω)))

    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    factor = mq / (8 * π * q)

    if ω^2 - q^2 > (mq + mD)^2 && abs(ω) < cutoff
        x1 = x_pm(abs(ω), q, mq, mD, -1)
        x2 = x_pm(abs(ω), q, mq, mD, +1)
        return factor * sign(ω) * (term1(x2) - term1(x1))
    end
    # for now ignoring the landau damping terms
    # if ω^2 - q^2 < 0
    #     return term3(cutoff) - term3(x1)
    # end
    # if sqrt(q^2 + (mD - mq)^2) > ω > q
    #     return term3(x2) - term3(x1)
    # end
    return 0.0
end

function imagpart_baryon_q(T, mu, om, q, mq, mD, param)
    ω = om + 3mu
    if q == 0
        return imagpart_baryon_q0(T, mu, ω, mq, mD, param)
    end
    expb(eps) = exp(eps / (2T))

    term1(x) = -2 * T * log((expb(x) * (1 - expb(-(-ω - mu + x)))) / (1 + expb(-(-ω - mu - x))))
    term2(x) = 2 * T * log((expb(x) * (1 - expb(-(ω + mu + x)))) / (1 + expb(-(ω + mu - x))))

    term3(x) = 2 * T * log((1 - expb(-(-ω - mu + x))) / (1 + expb(-(ω + mu + x))))
    term4(x) = -2 * T * log((1 - expb(-(ω + mu + x))) / (1 + expb(-(-ω - mu + x))))

    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    factor = 2 * mq / (8 * π * q)

    s = ω^2 - q^2

    if abs(ω) >= cutoff || s == 0
        return 0.0
    end

    if s > 0
        if (mD - mq)^2 <= s <= (mq + mD)^2
            return 0.0
        end
        x1 = x_pm(abs(ω), q, mD, mq, -1)
        x2 = x_pm(abs(ω), q, mD, mq, +1)
        if s > (mq + mD)^2
            if ω >= 0
                return factor * (term2(x2) - term2(x1))
            else
                return factor * (term1(x2) - term1(x1))
            end
        else
            if ω >= 0
                return factor * (term4(x2) - term4(x1))
            else
                return factor * (term3(x2) - term3(x1))
            end
        end
    else
        x2 = x_pm(ω, q, mD, mq, +1)
        return factor * (term3(x2) + term4(x2))
    end
end

function imagpart_baryon_q_asymmetric(T, mu, om, q, mq, mD, param)
    ω = om + 3mu
    if q == 0
        return imagpart_baryon_q0(T, mu, ω, mq, mD, param)
    end
    expb(eps) = exp(eps / (2T))

    term1(x) = -2 * T * log((expb(x) * (1 - expb(-(-ω - mu + x)))) / (1 + expb(-(-ω - mu - x))))
    term2(x) = 2 * T * log((expb(x) * (1 - expb(-(ω + mu + x)))) / (1 + expb(-(ω + mu - x))))

    term3(x) = 2 * T * log((1 - expb(-(-ω - mu + x))) / (1 + expb(-(ω + mu + x))))
    term4(x) = -2 * T * log((1 - expb(-(ω + mu + x))) / (1 + expb(-(-ω - mu + x))))

    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    factor = 2 * mq / (8 * π * q)

    s = ω^2 - q^2

    if abs(ω) >= cutoff || s == 0
        return 0.0
    end

    if s > 0
        if (mD - mq)^2 <= s <= (mq + mD)^2
            return 0.0
        end
        x1 = x_pm(abs(ω), q, mD, mq, -1)
        x2 = min(x_pm(abs(ω), q, mD, mq, +1), 2 * sqrt((param.Λ + q)^2 + mD^2) - abs(ω))
        if s > (mq + mD)^2
            if ω >= 0
                return factor * (term2(x2) - term2(x1))
            else
                return factor * (term1(x2) - term1(x1))
            end
        else
            if ω >= 0
                return factor * (term4(x2) - term4(x1))
            else
                return factor * (term3(x2) - term3(x1))
            end
        end
    else
        x2 = x_pm(ω, q, mD, mq, +1)
        return factor * (term3(x2) + term4(x2))
    end
end

function imagpart_baryon_q0(T, mu, ω, mq, mD, param)
    if ω == 0
        return 0.0
    end
    eq = 0.5 * abs((ω^2 + mq^2 - mD^2) / ω)
    eD = 0.5 * abs((ω^2 - mq^2 + mD^2) / ω)
    if eD < 1e-5
        @show mq, mD, ω
    end
    p = sqrt(abs(mq^4 + (mD^2 - ω^2)^2 - 2mq^2 * (mD^2 + ω^2))) / (2 * abs(ω))
    factor = -p * mq / (π)
    term1 = (1 - numberF(T, mu, eq) + numberB(T, mu, eD)) / (eq + eD)
    term2 = -(1 - numberF(T, -mu, eq) + numberB(T, -mu, eD)) / (eq + eD)
    term3 = -(numberF(T, -mu, eq) + numberB(T, mu, eD)) * PrincipalValue(eD - eq)
    term4 = (numberF(T, mu, eq) + numberB(T, -mu, eD)) * PrincipalValue(eD - eq)
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

function coupling_B(T, mu, mq, mD, param)
    pol_D(ω) = realpart_D_normal_dependent_part(T, mu, ω, 0.0, mq, param)
    der = UsefulFunctions._derivative(pol_D, mD)
    # return 16 * mD / (mq * abs(der))
    return 82.284
end

function Π0_B(T, mu, mq, mD, param)
    factor = mq / (π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    integrand1(p) = p^2 * ((1 - numberF(T, mu, eq(p)) + numberB(T, mu, eD(p))) + (1 - numberF(T, -mu, eq(p)) + numberB(T, -mu, eD(p)))) / (eq(p) * eD(p) * (eq(p) + eD(p)))
    integrand2(p) = p^2 * (-(numberF(T, -mu, eq(p)) + numberB(T, mu, eD(p))) - (numberF(T, mu, eq(p)) + numberB(T, -mu, eD(p)))) * PrincipalValue(eD(p) - eq(p)) / (eq(p) * eD(p))

    return 1 / (coupling_B(T, mu, mq, mD, param)) - factor * integrate(p -> integrand1(p) + integrand2(p), 0.0, param.Λ)
end

function full_real_part_baryon(T, mu, ω, mq, mD, param)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    term1(p) = (1 - numberF(T, mu, eq(p)) + numberB(T, mu, eD(p))) * PrincipalValue(ω + 3 * mu + eq(p) + eD(p))
    term2(p) = -(1 - numberF(T, -mu, eq(p)) + numberB(T, -mu, eD(p))) * PrincipalValue(ω + 3 * mu - eq(p) - eD(p))
    term3(p) = -(numberF(T, -mu, eq(p)) + numberB(T, mu, eD(p))) * PrincipalValue(ω + 3 * mu - eq(p) + eD(p))
    term4(p) = (numberF(T, mu, eq(p)) + numberB(T, -mu, eD(p))) * PrincipalValue(ω + 3 * mu + eq(p) - eD(p))

    factor(p) = mq * p^2 / (π^2 * eq(p) * eD(p))

    integrand(p) = factor(p) * (term1(p) + term2(p) + term3(p) + term4(p))

    return 1 / coupling_B(T, mu, mq, mD, param) - integrate(integrand, 0.0, param.Λ)
end

function realpart_baryon_q(T, mu, ω, q, mq, mD, param)
    impart(x, y) = imagpart_baryon_q(T, mu, x, y, mq, mD, param)

    cutoff = max(sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2)), sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2))
    repart_dependent = realpart_kramers_kronig_q_1(impart, ω, q, -cutoff - 3mu, cutoff - 3mu)

    return Π0_B(T, mu, mq, mD, param) - repart_dependent
end

function realpart_baryon_KK(T, mu, ω, q, mq, mD, param)
    impart(x) = imagpart_baryon_q0(T, mu, x, mq, mD, param)

    cutoff = max(sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2)), sqrt(param.Λ^2 + mq^2) + sqrt(param.Λ^2 + mD^2), 2.0)
    repart_dependent = realpart_kramers_kronig_mod(impart, ω, cutoff)

    return Π0_B(T, mu, mq, mD, param) - repart_dependent
end

function realpart_baryon_q0(T, mu, ω, mq, mD, param)
    factor = mq / (π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    term1(p) = (1 - numberF(T, mu, eq(p)) + numberB(T, mu, eD(p))) * PrincipalValue(ω + eq(p) + eD(p))
    term2(p) = (1 - numberF(T, -mu, eq(p)) + numberB(T, -mu, eD(p))) * PrincipalValue(ω - eq(p) - eD(p))
    term3(p) = (numberF(T, -mu, eq(p)) + numberB(T, mu, eD(p))) * PrincipalValue(ω - eq(p) + eD(p))
    term4(p) = (numberF(T, mu, eq(p)) + numberB(T, -mu, eD(p))) * PrincipalValue(ω + eq(p) - eD(p))

    integrand(p) = p^2 * (term1(p) - term2(p) - term3(p) + term4(p)) / (eq(p) * eD(p))

    return 1 / coupling_B(T, mu, mq, mD, param) - factor * integrate(integrand, 0.0, param.Λ)
end

function mass_diquark_normal(T, mu, param)
    return find_mass_D(T, mu, param)[1]
end

function phase_shift_baryon_q(T, mu, ω, q, param)
    mq, ome = massgap(T, mu, param).zero
    mD = mass_diquark_normal(T, mu, param)
    mus = mu + ome

    impart = imagpart_baryon_q(T, mus, ω, q, mq, mD, param)
    repart = realpart_baryon_q(T, mus, ω, q, mq, mD, param)

    return atan(impart, repart)
end

function baryon_mass(T, mu_0, param)
    m, ome = massgap(T, mu_0, param).zero
    mu = mu_0 + ome
    mD = find_mass_D(T, mu, param)[1]
    rep(ω) = realpart_baryon_q(T, mu, ω, 0.0, m, mD, param)
    # if m > mu && rep(0.0) * rep(m + mD) < 0.0
    #     return fzero(rep, 0.4)
    # end
    return UsefulFunctions.bisection(rep, 0.0, m + mD)
    # function ff!(F, x)
    #     term = Πq0_D_analytic(T, mu, x[1] - 1im * x[2] / 2, m, param)
    #     F[1] = real(term)
    #     F[2] = imag(term)
    # end
    # return nlsolve(ff!, [2 * (m - mu), 0.1]).zero
    # return 0.0
end

function realpart_baryon_q0_1(T, mu, ω, mq, mD, param)
    factor = mq / (π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    term1(p) = (1 - numberF(T, mu, eq(p)) + numberB(T, mu, eD(p))) / (ω + eq(p) + eD(p))
    term2(p) = (1 - numberF(T, -mu, eq(p)) + numberB(T, -mu, eD(p))) / (ω - eq(p) - eD(p))
    term3(p) = (numberF(T, -mu, eq(p)) + numberB(T, mu, eD(p))) / (ω - eq(p) + eD(p))
    term4(p) = (numberF(T, mu, eq(p)) + numberB(T, -mu, eD(p))) / (ω + eq(p) - eD(p))

    integrand(p) = p^2 * (term1(p) - term2(p) - term3(p) + term4(p)) / (eq(p) * eD(p))

    return 1 / coupling_B(T, mu, mq, mD, param) - factor * integrate(integrand, 0.0, param.Λ)
end

function distribution_baryon_q(T, mu, q, param)
    stat_factor(om) = numberF(T, 3mu, om) * (1 - numberF(T, 3mu, om)) + numberF(T, -3mu, om) * (1 - numberF(T, -3mu, om))

    integrand(om) = stat_factor(om) * phase_shift_baryon_q(T, mu, om, q, param)

    return integrate(integrand, 0.0, 2.0)
end
