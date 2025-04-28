function integrate(func, a, b; rtol=1e-3, maxevals=1e4)
    return quadgk(func, a, b, rtol=rtol, maxevals=maxevals)[1]
end

function integrate(func, start::Vector, finish::Vector; rtol=1e-3, maxevals=100000)
    return hcubature(func, start, finish, reltol=rtol, maxevals=maxevals)[1]
end

function fzero(f, guess)
    sol = nlsolve(x -> f(x...), [guess])
    if sol.f_converged
        return sol.zero[1]
    else
        return 0.0
    end
end

"""

    realpart_kramers_kronig_q(imagpart::Function, ω, q, cutoff)

This function returns the real part from imaginary part of a function f
using the Kramers Kronig relation. Note that the imaginary part have to
a function of two variables ω and q. Also the function only returns the
ω, q dependent part.
"""
function realpart_kramers_kronig_q(imagpart::Function, ω, q, cutoff)
    integrand(ν) = 2 * ν * (imagpart(ν, q) * PrincipalValue(ν^2 - ω^2) - imagpart(ν, 0.0) * PrincipalValue(ν^2)) / π
    return integrate(integrand, 0.0, cutoff)
end
