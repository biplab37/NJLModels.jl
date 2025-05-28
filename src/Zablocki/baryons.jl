function x_pm(ω, q, m1, m2, sign)
    s = ω^2 - q^2
    return (m1^2 - m2^2) * ω / s + sign * q * sqrt((1 - (m1 - m2)^2 / s) * (1 - (m1 + m2)^2 / s))
end

function imagpart_baryon_q(T, mu, ω, q, mq, mD, param)
    if q == 0
        q = 0.01
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

function coupling_B(T, mu, mq, mD, param)
    pol_D(ω) = realpart_D_normal(T, mu, ω, 0.0, mq, param)
    der = UsefulFunctions._derivative(pol_D, mD)
    return -16 * mD / (mq * der)
    # return 90
end

function Π0_B(T, mu, mq, mD, param)
    factor = mq / (2π^2)
    eq(p) = sqrt(p^2 + mq^2)
    eD(p) = sqrt(p^2 + mD^2)

    integrand(p) = p^2 * ((1 - numberF(T, mu, eq(p)) - numberF(T, -mu, eq(p))) + (1 + numberB(T, mu, eD(p)) + numberB(T, -mu, eD(p)))) / (eq(p) * eD(p) * (eq(p) + eD(p)))

    return 1 / (coupling_B(T, mu, mq, mD, param)) + factor * integrate(integrand, 0.0, param.Λ)
end

function realpart_baryon_q(T, mu, ω, q, mq, mD, param)
    impart(x, y) = imagpart_baryon_q(T, mu, x, y, mq, mD, param)

    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    repart_dependent = realpart_kramers_kronig_q_1(impart, ω, q, -cutoff, cutoff)

    return Π0_B(T, mu, mq, mD, param) - repart_dependent
end

function mass_diquark_normal(T, mu, m, param)
    return find_mass_D(T, mu, m, param)[1]
end

function phase_shift_baryon_q(T, mu, ω, q, param)
    mq, ome = massgap(T, mu, param).zero
    mD = mass_diquark_normal(T, mu, mq, param)
    mus = mu + ome

    impart = imagpart_baryon_q(T, mus, ω, q, mq, mD, param)
    repart = realpart_baryon_q(T, mus, ω, q, mq, mD, param)

    return atan(impart, repart)
end
