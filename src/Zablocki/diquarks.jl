# Normal phase

ω_D(ω, μ, sign) = ω + 2 * sign * μ

s_D(ω, q, mu, sign) = ω_D(ω, mu, sign)^2 - q^2

E_D(sign_up, sign_down, ω, q, mu, m) = 0.5 * (ω_D(ω, mu, sign_down) + sign_up * q * sqrt(1 - (4 * m^2 / s_D(ω, q, mu, sign_down))))

function J_D_pair(sign, T, mu, ω, q, m)
    NN(sign_up, sign_fermi, sign_energy) = numberF(T, -sign * sign_fermi * mu, sign_energy * E_D(sign_up, sign, ω, q, mu, m))

    return T * log((NN(-1, +1, -1) * NN(-1, -1, +1)) / (NN(+1, +1, -1) * NN(+1, -1, +1)))
end

function J_D_Landau(sign, T, mu, ω, q, m)
    NN(sign_up, sign_fermi, sign_energy) = numberF(T, -sign * sign_fermi * mu, sign_energy * E_D(sign_up, sign, ω, q, mu, m))

    return 2 * T * log(NN(-1, -1, +1) / NN(+1, +1, -1))
end

function imagpart_D_normal_q0(T, mu, ω, m, param)
    factor = (3 - 1) * 2 / (8π)
    positive_term, negative_term = 0.0, 0.0
    if ω * ((ω / 4) + mu) > m^2 - mu^2
        positive_term = (ω + 2mu)^2 * sqrt(1 - 4 * m^2 / (ω + 2mu)^2)
    end
    if ω * ((ω / 4) - mu) > m^2 - mu^2
        negative_term = (ω - 2mu)^2 * sqrt(1 - 4 * m^2 / (ω - 2mu)^2)
    end

    return factor * (1 - 2 * numberF(T, 0.0, ω / 2)) * (positive_term + negative_term)
end

function imagpart_D_normal(T, mu, ω, q, m, param)
    if ω^2 > 4 * (param.Λ^2 + m^2 + (q^2 / 4))
        return 0.0
    end

    if q == 0.0
        return imagpart_D_normal_q0(T, mu, ω, m, param)
    end

    factor = 2 / (4 * π * q)

    # postive branch
    sp = s_D(ω, q, mu, +1)

    if sp > 4m^2
        positive = sp * J_D_pair(+1, T, mu, ω, q, m)
    elseif sp < 0
        positive = sp * J_D_Landau(+1, T, mu, ω, q, m)
    else
        positive = 0.0
    end

    #negative branch
    sm = s_D(ω, q, mu, -1)

    if sm > 4m^2
        negative = sm * J_D_pair(-1, T, mu, ω, q, m)
    elseif sm < 0
        negative = sm * J_D_Landau(-1, T, mu, ω, q, m)
    else
        negative = 0.0
    end

    return factor * (positive + negative)
end

function Π0_D(T, mu, m, param)
    factor = (3 - 1) * 2 / (2 * π^2)
    function integrand_pi0(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * numberF(T, -mu, ep)) * PrincipalValue(ep + mu)) + (1 - 2 * numberF(T, mu, ep)) * PrincipalValue(ep - mu))
    end

    return 1 / (2 * param.GD) + factor * integrate(integrand_pi0, m, sqrt(param.Λ^2 + m^2))
end

function realpart_D_normal(T, mu, ω, q, m, param)
    impart(x, y) = imagpart_D_normal(T, mu, x, y, m, param)

    cutoff = 2 * sqrt(param.Λ^2 + m^2 + (q^2 / 4))
    repart_dependent = realpart_kramers_kronig_q(impart, ω, q, cutoff)

    return Π0_D(T, mu, m, param) - repart_dependent
end

function phase_shift_D_normal(T, mu, ω, q, param)
    m, ome = massgap(T, mu, param).zero
    mus = mu + ome

    impart = imagpart_D_normal(T, mus, ω, q, m, param)
    repart = realpart_D_normal(T, mus, ω, q, m, param)

    return atan(impart, repart)
end

function integrand_pressure_D_normal(T, mu, ω, q, param)
    return 2 * numberB(T, 0.0, ω) * phase_shift_D_normal(T, mu, ω, q, param) / (2 * π^2)
end

function pressure_fl_D_normal(T, mu, param)
    integrand_fl_D(x) = integrand_pressure_D_normal(T, mu, x[1], x[2], param)

    return integrate(integrand_fl_D, [0.0, 0.0], [sqrt(5) * param.Λ, param.Λ])
end
