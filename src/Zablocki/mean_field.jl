function mean_field_pressure(T, mu, param)
    m1, ω1 = massgap(T, mu, param).zero
    m2, ω2 = massgap(0.01, 0.0, param).zero

    scalar_part(m, ω) = (m - param.m0)^2 / (4 * param.Gs) - ω^2 / (4 * param.Gv)
    factor = 2 * 2 * 3 / (2 * π^2)
    En(p, m) = sqrt(p^2 + m^2)
    divergentpart(p, m) = p^2 * En(p, m)
    mediumpart(p, m, ω) = T * p^2 * (log(1 + exp(-(En(p, m) - mu - ω) / T)) + log(1 + exp(-(En(p, m) + mu + ω) / T)))

    div_part = integrate(q -> divergentpart(q, m1), 0.0, param.Λ)
    div_part_zero = integrate(q -> divergentpart(q, m2), 0.0, param.Λ)
    medium_term = integrate(q -> mediumpart(q, m1, ω1), 0.0, 10.0)

    return -scalar_part(m1, ω1) + scalar_part(m2, ω2) + factor * (div_part + medium_term - div_part_zero)
end

function mean_field_pressure_m0(T, mu, param)
    m1, ω1 = massgap(T, mu, param).zero
    m2, ω2 = massgap(0.01, 0.0, param).zero

    scalar_part(m, ω) = (m - param.m0)^2 / (4 * param.Gs) - ω^2 / (4 * param.Gv)
    factor = 2 * 2 * 3 / (2 * π^2)
    En(p, m) = sqrt(p^2 + m^2)
    divergentpart(p, m) = p^2 * En(p, m)
    mediumpart(p, m, ω) = T * p^2 * (log(1 + exp(-(En(p, m) - mu - ω) / T)) + log(1 + exp(-(En(p, m) + mu + ω) / T)))

    div_part = integrate(q -> divergentpart(q, m1), 0.0, param.Λ)
    div_part_zero = integrate(q -> divergentpart(q, m2), 0.0, param.Λ)
    medium_term = integrate(q -> mediumpart(q, m1, ω1), 0.0, 10.0)

    return -scalar_part(m1, ω1) + scalar_part(m2, ω2) + factor * (div_part + medium_term - div_part_zero)
end

function mean_field_pressure_scalar(T, mu, param)
    m1, ω1 = massgap(T, mu, param).zero
    scalar_part(m, ω) = (m - param.m0)^2 / (4 * param.Gs) - ω^2 / (4 * param.Gv)
    return -scalar_part(m1, ω1)
end

function mean_field_pressure_quark(T, mu, param)
    m1, ω1 = massgap(T, mu, param).zero
    m2, ω2 = massgap(0.01, 0.0, param).zero

    factor = 2 * 2 * 3 / (2 * π^2)
    En(p, m) = sqrt(p^2 + m^2)
    divergentpart(p, m) = p^2 * En(p, m)
    mediumpart(p, m, ω) = T * p^2 * (log(1 + exp(-(En(p, m) - mu - ω) / T)) + log(1 + exp(-(En(p, m) + mu + ω) / T)))

    div_part = integrate(q -> divergentpart(q, m1), 0.0, param.Λ)
    div_part_zero = integrate(q -> divergentpart(q, m2), 0.0, param.Λ)
    medium_term = integrate(q -> mediumpart(q, m1, ω1), 0.0, 10.0)

    return factor * [(div_part - div_part_zero), medium_term]
end

function mean_field_pressure_quark_nonrenormalized(T, mu, param, m1, ω1)
    factor = 2 * 2 * 3 / (2 * π^2)
    En(p, m) = sqrt(p^2 + m^2)
    divergentpart(p, m) = p^2 * En(p, m)
    mediumpart(p, m, ω) = T * p^2 * (log(1 + exp(-(En(p, m) - mu - ω) / T)) + log(1 + exp(-(En(p, m) + mu + ω) / T)))

    div_part = integrate(q -> divergentpart(q, m1), 0.0, param.Λ)
    medium_term = integrate(q -> mediumpart(q, m1, ω1), 0.0, 10.0)

    return factor * (medium_term)
end

function mean_field_density(T, mu, param)
    Pres(μ) = mean_field_pressure(T, μ, param)
    return UsefulFunctions._derivative(Pres, mu)
end

function mean_field_density1(T, mu, param)
    m, ω = massgap(T, mu, param).zero
    factor = 6 / π^2
    integrand(p) = p^2 * (numberF(T, mu + ω, En(p, m)) - numberF(T, -mu - ω, En(p, m)))
    return factor * integrate(integrand, 0.0, param.Λ)
end
