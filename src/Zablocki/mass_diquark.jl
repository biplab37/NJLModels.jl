function realpart_diquark_width(ω, Γ, m, T, mu, param)
    factor = 2 / (π^2)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * numberF(T, mu, ep)) * (ω / 2 + ep - mu) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) + ((1 - 2 * numberF(T, -mu, ep)) * (-ω / 2 + ep + mu) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end

    return 1 / (2 * param.GD) - factor * integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function imagpart_diquark_width(ω, Γ, m, T, mu, param)
    function integrand(ep)
        return sqrt(ep^2 - m^2) * ep * (((1 - 2 * numberF(T, mu, ep)) / ((ω / 2 + ep - mu)^2 + Γ^2 / 16)) - ((1 - 2 * numberF(T, -mu, ep)) / ((-ω / 2 + ep + mu)^2 + Γ^2 / 16)))
    end

    return (1 / (2 * π^2)) * Γ * integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end
