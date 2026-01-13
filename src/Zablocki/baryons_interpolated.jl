# usage of interpolation and other approximate improvement towards the calculation of Baryon properties

function baryon_phase_shift_interpolated(T, mu, q, mq, mD, param; o_length=200)
    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    ωrange = range(0, cutoff, length=o_length)

    imparts_data = [imagpart_baryon_q(T, mu, ω, q, mq, mD, param) for ω in ωrange]
    impart_interp = linear_interpolation(ωrange, imparts_data)
    impart_func(x) = (x >= cutoff) ? 0.0 : impart_interp(x)

    pi0 = Π0_B(T, mu, mq, mD, param)

    orange = 2 .* ωrange
    reparts_data = [pi0 - realpart_kramers_kronig(x -> impart_func(x), om, cutoff) for om in orange]
    repart_interp = linear_interpolation(orange, reparts_data)
    repart_func(x) = (x >= 2 * cutoff) ? pi0 : repart_interp(x)

    return x -> atan(impart_func(x), repart_func(x))
end

function integrand_number_density(T, mu, q, mq, mD, param)
    phase_shift(x, μ) = baryon_phase_shift_interpolated(T, μ, q, mq, mD, param)(x)

    integrand(x) = numberF(T, mu, x) * UsefulFunctions._derivative(μ -> phase_shift(x, μ), mu)

    return integrand
end

function distribution_function_q(T, mu, q, mq, mD, param)
    cutoff = sqrt(q^2 + 4 * param.Λ^2 + 2(mq^2 + mD^2))
    integrand = integrand_number_density(T, mu, q, mq, mD, param)

    return integrate(integrand, 0.0, cutoff)
end
