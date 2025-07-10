## Typical calculations with the gap equations, given in terms of the type Gaps

function solve_gap_equations(variables, gaps::Gaps, param::Parameters; guess=gaps.initial_guess)
    eqns(x) = gaps.gap_eqn(x, variables...)
    return nlsolve(eqns, guess).zero
end

function phase_diagram(variables_region, gaps::Gaps, param::Parameters)
    ##TODO: finish this function
    return nothing
end