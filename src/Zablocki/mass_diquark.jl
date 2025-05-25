function realpart_diquark_width(ω, Γ, m, T, mu, param)
    factor = 2 / (π^2)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * numberF(T, mu, ep)) * (ω / 2 + ep - mu) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) + ((1 - 2 * numberF(T, -mu, ep)) * (-ω / 2 + ep + mu) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end
    cutoff = sqrt(param.Λ^2 + m^2)

    return 1 / (2 * param.GD) - factor * integrate(integrand, m, 0.71)
end

function imagpart_diquark_width(ω, Γ, m, T, mu, param)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * numberF(T, mu, ep)) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) - ((1 - 2 * numberF(T, -mu, ep)) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end

    return (1 / (2 * π^2)) * Γ * integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I0(m, param)
    f(ep) = 2 * sqrt(ep^2 - m^2) / π^2
    return integrate(f, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I1_real(M, T, mu, m, param)
    f(ep) = -sqrt(ep^2 - m^2) * ep * ((numberF(T, -mu, ep)) / (2ep + (M + 2mu)) + (numberF(T, mu, ep)) / (2ep - (M + 2mu)))
    return 4 * integrate(f, m, sqrt(param.Λ^2 + m^2)) / π^2
end

function diquark_I1_imag(M, T, mu, m, param)
    if (M + 2mu) < 2m
        @info "Mass should be larger than 2m" T mu
        return 0.0
    end
    return sqrt((M + 2mu)^2 - 4 * m^2) * (M + 2 * mu) * numberF(T, 0.0, M / 2)
end

function diquark_I2_real(M, T, mu, m, param)
    f(ep) = -2 * sqrt(ep^2 - m^2) * PrincipalValue((M + 2 * mu)^2 - 4 * ep^2) / π^2
    return integrate(f, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I2_imag(M, T, mu, m, param)
    if (M + 2mu) < 2m
        @info "Mass should be larger than 2m" T mu
        return 0.0
    end
    return sqrt(1 - (4 * m^2 / (M + 2 * mu)^2)) / (4 * π)
end

function approx_diquark_mass(T, mu, param)
    m, ome = massgap(T, mu, param).zero
    mus = mu + ome

    I0 = diquark_I0(m, param)
    imI1(mD) = diquark_I1_imag(mD, T, mus, m, param)
    imI2(mD) = diquark_I2_imag(mD, T, mus, m, param)
    reI1(mD) = diquark_I2_imag(mD, T, mus, m, param)
    reI2(mD) = diquark_I2_imag(mD, T, mus, m, param)

    gamma_m(mD) = (imI2(mD) * (1 / (2 * param.GD) - reI1(mD) - I0) + imI1(mD) * reI2(mD)) / ((mD + 2 * mu) * (imI2(mD)^2 + reI2(mD)^2))
    term1(mD) = sqrt((reI2(mD) * (1 / (2 * param.GD) - reI1(mD) - I0) - imI1(mD) * imI2(mD)) / (imI2(mD)^2 + reI2(mD)^2) + gamma_m(mD)^2 / 4) - 2 * mu
    @show gamma_m(0.7)

    mass_D = bisection(term1, 2(m - mu) + 0.2, 0.8)
    return mass_D, gamma_m(mass_D), imI1(mass_D), imI2(mass_D), reI1(mass_D), reI2(mass_D), I0
end

function diquark_mass(T, mu, param, guess=[2 * massgap(T, mu, param).zero[1], 0.1])
    m, ome = massgap(T, mu, param).zero
    mus = mu + ome
    function ff!(F, x)
        (realpart_diquark_width(x[1], x[2], m, T, mus, param)^2 + imagpart_diquark_width(x[1], x[2], m, T, mus, param)^2)
    end

    return mcpsolve(ff!, [2(m - mu), 0.0], [1.0, 1.0], guess, iterations=5_000)
end
