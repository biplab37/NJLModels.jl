function integrate(func, a, b; rtol=1e-3, maxevals=1e4)
    return quadgk(func, a, b, rtol=rtol, maxevals=maxevals)[1]
end

function fzero(f, guess)
    sol = nlsolve(x -> f(x...), [guess])
    if sol.f_converged
        return sol.zero[1]
    else
        return 0.0
    end
end
