#### This file contains code related to imaginary part calculation of the polarization function

## Zero External Momenta

# Equal masses
function I_equal_mass_q0(F::Function, ω, m, xi)
    if ω^2 <= 4m^2 || xi * ω < 0
        return 0.0
    end

    return -sqrt(1 - 4m^2 / ω^2) * F(abs(ω) / 2, abs(ω) / 2) / (4π)
end

# Different Masses, NOTE: m1>m2
function I1_different_mass_q0(F::Function, ω, m1, m2, xi)
    if xi * ω < 0 || abs(ω) < (m1 + m2)
        return 0.0
    end
    factor = -sqrt((1 - (m1 + m2)^2 / (ω^2))(1 - (m1 - m2)^2 / (ω^2))) / (4π)
    return factor * F((ω^2 + m1^2 - m2^2) / (2abs(ω)), (ω^2 - m1^2 + m2^2) / (2ω))
end

function I2_different_mass_q0(F::Function, ω, m1, m2, xi)
    if xi * ω < 0 || abs(ω) > (m1 - m2)
        return 0.0
    end
    factor = sqrt((1 - (m1 + m2)^2 / (ω^2))(1 - (m1 - m2)^2 / (ω^2))) / (4π)
    return factor * F((ω^2 + m1^2 - m2^2) / (2abs(ω)), -(ω^2 - m1^2 + m2^2) / (2ω))
end

## Finite Momenta

# Equal masses
function I1_equal_masses_q(F::Function, ω::Number, q::Number, m::Number, xi)
    if xi * ω < 0 || ω^2 - q^2 <= 4m^2
        return 0.0
    end

    integrand(E2) = F((abs(ω) + E2) / 2, (abs(ω) - E2) / 2)

    y = q * sqrt(1 - 4m^2 / (ω^2 - q^2))

    return -integrate(integrand, -y, y) / (8π * q)
end

function I2_equal_masses_q(F::Function, ω, q, m, xi, max_integration_limit=10.0)
    if ω^2 - q^2 >= 0.0
        return 0.0
    end

    y = q * sqrt(1 - 4m^2 / (ω^2 - q^2))

    integrand(E1) = F((E1 + xi * ω) / 2, (E1 - xi * ω) / 2)

    return -integrate(integrand, y, max_integration_limit) / (8π * q)
end

# Different Masses
function x_pm(ω, q, m1, m2, sign)
    if m1 < m2
        error("Check the order of masses!")
    end
    s = ω^2 - q^2
    return (m1^2 - m2^2) * ω / s + sign * q * sqrt((1 - (m1 - m2)^2 / s) * (1 - (m1 + m2)^2 / s))
end

function I1_different_masses_q(F::Function, ω, q, m1, m2, xi)
    if xi * ω < 0 || ω^2 - q^2 <= (m1 + m2)^2
        return 0.0
    end

    xp = x_pm(ω, q, m1, m2, +1)
    xm = x_pm(ω, q, m1, m2, -1)

    integrand(E2) = F((abs(ω) + E2) / 2, (abs(ω) - E2) / 2)

    return -integrate(integrate, xm, xp) / (8π * q)
end

function I2_different_masses_q_part1(F::Function, ω, q, m1, m2, xi, max_integration_limit=10.0)
    if ω^2 - q^2 >= 0.0
        return 0.0
    end

    xp = x_pm(ω, q, m1, m2, +1)

    integrand(E1) = F((E1 + xi * ω) / 2, (E1 - xi * ω) / 2)

    return -integrate(integrate, xp, max_integration_limit) / (8π * q)
end

function I2_different_masses_q_part2(F::Function, ω, q, m1, m2, xi)
    e2max = sqrt(q^2 + (m1 - m2)^2)

    if xi > 0
        if (ω - q) * (e2max - ω) <= 0
            return 0.0
        end
    else
        if (ω + q) * (e2max + ω) >= 0
            return 0.0
        end
    end

    xp = x_pm(ω, q, m1, m2, +1)
    xm = x_pm(ω, q, m1, m2, -1)

    integrand(E1) = F((E1 + xi * ω) / 2, (E1 - xi * ω) / 2)

    return integrate(integrate, xm, xp) / (8π * q)
end
