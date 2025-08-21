# This file contains code to calculate the Baryon self energy from the quark diquark pair.
# We use the full spectral function calculation of the diquark propagator and use the
# matsubara sum analytically to avoid numerical issues.

function identity_part_im_self_B(T, mu, m, p0, param)
    func(Ek) = m * sqrt(Ek^2 - m^2) * (spectral_function_D_normal(T, mu, p0 + mu - Ek, sqrt(Ek^2 - m^2), param) * (numberB(T, -mu, p0 - Ek) + numberF(T, mu, Ek) - 1.0) + spectral_function_D_normal(T, mu, p0 + mu + Ek, sqrt(Ek^2 - m^2), param) * (numberF(T, -mu, Ek) - numberB(T, -mu, p0 + Ek))) / (4π)

    return integrate(func, m, sqrt(param.Λ^2 + m^2))
end

function gamma0_part_im_self_B(T, mu, m, p0, param)
    func(Ek) = Ek * sqrt(Ek^2 - m^2) * (spectral_function_D_normal(T, mu, p0 + mu - Ek, sqrt(Ek^2 - m^2), param) * (numberB(T, -mu, p0 - Ek) + numberF(T, mu, Ek) - 1.0) + spectral_function_D_normal(T, mu, p0 + mu + Ek, sqrt(Ek^2 - m^2), param) * (-numberF(T, -mu, Ek) + numberB(T, -mu, p0 + Ek))) / (4π)

    return integrate(func, m, sqrt(param.Λ^2 + m^2))
end

function total_part_im_self_B(T, mu, m, p0, param)
    func(Ek) = Ek * sqrt(Ek^2 - m^2) * (spectral_function_D_normal(T, mu, p0 + mu - Ek, sqrt(Ek^2 - m^2), param) * (numberB(T, -mu, p0 - Ek) * (1 + m / Ek) + numberF(T, mu, Ek) - 1.0) + spectral_function_D_normal(T, mu, p0 + mu + Ek, sqrt(Ek^2 - m^2), param) * (-numberF(T, -mu, Ek) + numberB(T, -mu, p0 + Ek)) * (1 - m / Ek)) / (4π)

    return integrate(func, m, sqrt(param.Λ^2 + m^2))
end

function realpart_self_B(T, mu, m, p0, param)
    func(x) = total_part_im_self_B(T, mu, m, x, param)

    return realpart_kramers_kronig_1(func, p0, -2.0, 2.0, rtol=1e-1, maxevals=1000)
end

function imagpart_baryons(T, p0, m, mD, ΓD, param)#At zero chemical potential
    rho(om, k) = spectral_function_D_normal(T, 0.0, om, k, param)
    k(en) = sqrt(en^2 - m^2)
    f(x) = numberF(T, 0.0, x)
    g(x) = numberB(T, 0.0, x)
    func(en) = m * k(en) * (rho(p0 - en, k(en)) * (1 - f(en) + g(p0 - en)) - rho(p0 + en, k(en)) * (f(en) + g(p0 + en))) / π
    impa = (ΓD != 0) ? 0.0 : imagpart_baryon_q0(T, 0.0, p0, m, mD, param)
    return 0.25 * integrate(func, m, sqrt(param.Λ^2 + m^2)) + impa
end

function realpart_baryons(T, p0, m, mD, ΓD, param)
    func(x) = imagpart_baryons(T, x, m, mD, ΓD, param)

    return realpart_kramers_kronig_1(func, p0, -4.0, 4.0, rtol=1e-2, maxevals=500)
end
