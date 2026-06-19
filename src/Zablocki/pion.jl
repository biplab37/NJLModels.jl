# function Π0_meson(T, mu, m, NM, param::Parameters; Nc=3, Nf=2)
#     pauli(p) = 1 - FD_dist(T, mu, En(p, m)) - FD_dist(T, -mu, En(p, m))
#     return 1 / (2 * param.Gs) - Nc * Nf * (integrate(p -> (1 - m^2 / En(p, m)^2)^(NM - 0.5) * p^2 * pauli(p) / En(p, m), 0, param.Λ)) / (1 * π^2)
# end
# function imagpart_meson_q0(T, mu, ω::Real, m, NM, param::Parameters; Nc=3, Nf=2)
#     # NM=1/2 for pion and 3/2 ofr sigma
#     if ω^2 < 4 * m^2 || ω^2 > 4 * (param.Λ^2 + m^2)
#         return 0.0
#     end
#     return imagpart_meson_q0_C(T, mu, ω, m, NM, Nc=Nc, Nf=Nf)
# end

# function imagpart_meson_q0_C(T, mu, ω, m, NM; Nc=3, Nf=2)
#     # NM=1/2 for pion and 3/2 for sigma
#     factor = Nc * Nf / 8π
#     pauli_term = 1 - FD_dist(T, mu, ω / 2) - FD_dist(T, -mu, ω / 2)
#     return factor * ω^2 * ((1 - (4 * m^2 / ω^2))^NM) * pauli_term
# end

# function Π0_meson_analytic(T, mu, ω, m, NM, param::Parameters; Nc=3, Nf=2)
#     pauli(p) = 1 - FD_dist(T, mu, En(p, m)) - FD_dist(T, -mu, En(p, m))
#     ep(p) = En(p, m)
#     repart = 1 / (2 * param.Gs) - 2 * Nc * Nf * (integrate(p -> (1 - m^2 / ep(p)^2)^(NM - 0.5) * p^2 * pauli(p) * (1 / (ep(p) + ω / 2) + 1 / (ep(p) - ω / 2)), 0, param.Λ)) / (4 * π^2)
#     if imag(ω) > 0
#         return repart
#     elseif imag(ω) < 0
#         return repart - 2im * imagpart_meson_q0_C(T, mu, ω, m, NM, Nc=Nc, Nf=Nf)
#     else
#         return realpart_meson_q0(T, mu, ω, m, NM, param, Nc=Nc, Nf=Nf)
#     end
# end

function find_meson_mass(T, mu, NM, m, param::Parameters, Phi=0.0, Phibar=0.0, guess=[0.2, 0.01])
    repart(ω) = Π0_meson_analytic(T, mu, ω, m, NM, param, Phi, Phibar)

    if repart(0.0) * repart(2 * m) < 0 # there is a bound state with 0<M<2m
        return [bisection(repart, 2 * m, 0.0), 0.0]
    end

    # return [fzero(repart, 0.1), 0.0]
    # println()
    function ff!(F, x)
        term = repart(x[1] - 0.5im * x[2]) # Look for pole in the second Reimann sheet.
        F[1] = real(term)
        F[2] = imag(term)
    end
    return nlsolve(ff!, guess).zero
end
find_pion_mass(T, mu, m, param, Phi=0.0, Phibar=0.0, guess=[0.1, 0.01]) = find_meson_mass(T, mu, 0.5, m, param, Phi, Phibar, guess)
find_sigma_mass(T, mu, m, param, Phi=0.0, Phibar=0.0, guess=[0.2, 0.01]) = find_meson_mass(T, mu, 1.5, m, param, Phi, Phibar, guess)