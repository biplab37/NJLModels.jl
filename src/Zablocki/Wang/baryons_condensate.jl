
function quark_propagator_without_condensate11(p0, k0, k, mu, m)
    return ((p0 - k0 - mu) + m) * PrincipalValue((p0 - k0 + mu)^2 - (k^2 + m^2), 0.05)
end
function quark_propagator_without_condensate22(p0, k0, k, mu, m)
    return (-(p0 - k0 - mu) + m) / ((p0 - k0 + mu)^2 - (k^2 + m^2))
end

function baryon_self_energy_without_condensate(T, mu, p0, param)
    m, ome = massgap(T, mu, param).zero
    mu += ome
    integrand(k0, k) = k^2 * quark_propagator_without_condensate11(p0, k0, k, mu, m) * diquark_propagator_without_condensate(T, mu, k0, k, m, param)
    return integrate(x -> integrand(x[1], x[2]), [-2.0, 0], [0.0, param.Λ], rtol=1e-1)
end
