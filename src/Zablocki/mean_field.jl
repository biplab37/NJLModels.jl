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