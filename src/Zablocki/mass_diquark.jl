function realpart_diquark_width(ω, Γ, m, T, mu, param)
    factor = 2 / (π^2)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * FD_dist(T, mu, ep)) * (ω / 2 + ep - mu) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) + ((1 - 2 * FD_dist(T, -mu, ep)) * (-ω / 2 + ep + mu) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end
    cutoff = sqrt(param.Λ^2 + m^2)

    return 1 / (2 * param.GD) - factor * integrate(integrand, m, cutoff)
end

function imagpart_diquark_width(ω, Γ, m, T, mu, param)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * FD_dist(T, mu, ep)) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) - ((1 - 2 * FD_dist(T, -mu, ep)) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end

    return -(1 / (2 * π^2)) * Γ * integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I0(m, param)
    f(ep) = 2 * sqrt(ep^2 - m^2) / π^2
    return integrate(f, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I1_real(M, T, mu, m, param)
    f(ep) = -sqrt(ep^2 - m^2) * ep * ((FD_dist(T, -mu, ep)) / (2ep + (M + 2mu)) + (FD_dist(T, mu, ep)) / (2ep - (M + 2mu)))
    return 4 * integrate(f, m, sqrt(param.Λ^2 + m^2)) / π^2
end

function diquark_I1_imag(M, T, mu, m)
    if (M + 2mu) < 2m
        @info "Mass should be larger than 2m" T mu
        return 0.0
    end
    return sqrt((M + 2mu)^2 - 4 * m^2) * (M + 2 * mu) * FD_dist(T, 0.0, M / 2)
end

function diquark_I2_real(M, T, mu, m, param)
    f(ep) = -2 * sqrt(ep^2 - m^2) * PrincipalValue((M + 2 * mu)^2 - 4 * ep^2) / π^2
    return integrate(f, m, sqrt(param.Λ^2 + m^2))
end

function diquark_I2_imag(M, T, mu, m)
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
    imI1(mD) = diquark_I1_imag(mD, T, mus, m)
    imI2(mD) = diquark_I2_imag(mD, T, mus, m)
    reI1(mD) = diquark_I2_real(mD, T, mus, m, param)
    reI2(mD) = diquark_I2_real(mD, T, mus, m, param)

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
        F[1] = realpart_diquark_width(x[1], x[2], m, T, mus, param)^2
        F[2] = imagpart_diquark_width(x[1], x[2], m, T, mus, param)^2
    end

    return mcpsolve(ff!, [2(m - mu), 0.0], [1.0, 1.0], guess, iterations=5_000)
end

function find_mass_D(T, mu_0, param)
    m = massgap_m(T, mu_0, param)
    mu = mu_0
    rep(ω) = realpart_D_normal(T, mu, ω - 2mu, 0.0, m, param)
    if m > mu && rep(0.0) * rep(2(m)) < 0.0
        return bisection(rep, 0.0, 2(m)) - 2mu, 0.0
    end
    function ff!(F, x)
        term = Πq0_D_analytic(T, mu, x[1] - 1im * x[2] / 2 - 2mu, m, param)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return mcpsolve(ff!, [0.0, 0.0], [2.0, 1.0], [max(0.1, 0.6), 0.1]).zero - [2mu, 0.0]
end
function find_mass_D(T, mu, m, param)
    rep(ω) = realpart_D_normal(T, mu, ω - 2mu, 0.0, m, param)
    if m > mu && rep(0.0) * rep(2(m)) < 0.0
        return bisection(rep, 0.0, 2(m)) - 2mu, 0.0
    end
    function ff!(F, x)
        term = Πq0_D_analytic(T, mu, x[1] - 1im * x[2] / 2 - 2mu, m, param)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return mcpsolve(ff!, [0.0, 0.0], [2.0, 1.0], [max(0.1, 0.6), 0.1]).zero - [2mu, 0.0]
end

function find_mass_D_q(T, mu_0, q, param)
    m, ome = massgap(T, mu_0, param).zero
    mu = mu_0 + ome
    rep(ω) = realpart_D_normal(T, mu, ω, q, m, param)
    EQ = sqrt(q^2 + 4m^2)
    if EQ > mu && rep(0.0) * rep(EQ - 2mu) < 0.0
        return bisection(rep, 0.0, EQ - 2mu), 0.0
    end
    function ff!(F, x)
        term = Πq0_D_analytic(T, mu, x[1] - 1im * x[2] / 2, m, param)
        F[1] = real(term)
        F[2] = imag(term)
    end
    return nlsolve(ff!, [2 * (EQ - mu), 0.1]).zero
end

function find_diquark_energy_q(T, mu, q, m, param)
    EQ = sqrt(q^2 + 4m^2)

    if EQ > mu && realpart_D_normal(T, mu, 0.0, q, m, param) * realpart_D_normal(T, mu, EQ - 2mu, q, m, param) < 0.0
        return bisection(x -> realpart_D_normal(T, mu, x, q, m, param), 0.0, EQ - 2mu)
    end

    return 0.0
end
