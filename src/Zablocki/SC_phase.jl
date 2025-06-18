# This file contains code for 2SC phase at zero momentum.

function polarization_matrix_2SC(T, mu, ω, gap, param)
    pol_matrix = zeros(8, 8)

    pol_pi_pi = polarization_pi_pi(T, mu, ω, gap, param)
    pol_sigma_sigma = polarization_sigma_sigma(T, mu, ω, gap, param)
    pol_sigma_del2 = polarization_sigma_del2(T, mu, ω, gap, param)
    pol_del2_del2 = polarization_del2_del2(T, mu, ω, gap, param)
    pol_del2_star_del2 = polarization_del2_star_del2(T, mu, ω, gap, param)
    pol_del5_del5 = polarization_del5_del5(T, mu, ω, gap, param)

    pol_matrix[1, 1] = pol_pi_pi
    pol_matrix[2, 2] = pol_sigma_sigma
    pol_matrix[5, 5], pol_matrix[6, 6], pol_matrix[7, 7], pol_matrix[8, 8] = repeat([pol_del5_del5], 4)
end
