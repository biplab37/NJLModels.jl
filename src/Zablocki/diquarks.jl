# Normal phase

ω_D(ω, μ, sign) = ω + 2 * sign * μ

s_D(ω, q, mu, sign) = ω_D(ω, mu, sign)^2 - q^2

E_D(sign_up, sign_down, ω, q, mu, m) = 0.5 * (ω_D(ω, mu, sign_down) + sign_up * q * sqrt(1 - 4 * m^2 / s_D(ω, q, mu, sign_down)))

function J_D_pair(sign, T, mu, ω, q, m)
    NN(sign_up, sign_fermi, sign_energy) = numberF(T, -sign * sign_fermi * mu, sign_energy * E_D(sign_up, sign, ω, q, mu, m))

    return T * log((NN(-1, +1, -1) * NN(-1, -1, +1)) / (NN(+1, +1, -1) * NN(+1, -1, +1)))
end

function J_D_Landau(sign, T, mu, ω, q, m)
    NN(sign_up, sign_fermi, sign_energy) = numberF(T, -sign * sign_fermi * mu, sign_energy * E_D(sign_up, sign, ω, q, mu, m))

    return 2 * T * log(NN(-1, -1, +1) / NN(+1, +1, -1))
end

function imagpart_D_Normal(T, mu, ω, q, m, param)
    factor = 2 / (16 * π * q)

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
