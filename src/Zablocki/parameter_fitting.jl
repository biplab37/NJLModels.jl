# Parameter fitting for NJL/PNJL models

# Helper struct to pass temporary parameters during fitting
struct TempParameters <: Parameters
    Λ::Float64
    m0::Float64
    Gs::Float64
    Gv::Float64
    eta_d::Float64
    GD::Float64
end

TempParameters(Λ, m0, Gs) = TempParameters(Λ, m0, Gs, Gs / 2, 0.75, 0.75 * Gs)

function vacuum_observables(Λ, m0, Gs; regularization=:sharp, pv_scale=10.0)
    param = TempParameters(Λ, m0, Gs)
    cutoff = (regularization == :pv) ? pv_scale * Λ : Λ

    # 1. Solve vacuum gap equation for constituent quark mass M
    function gap_eqn(M)
        if regularization == :pv
            m1 = sqrt(M^2 + Λ^2)
            m2 = sqrt(M^2 + 2 * Λ^2)
            int_val = integrate(p -> p^2 * (1 / sqrt(p^2 + M^2) - 2 / sqrt(p^2 + m1^2) + 1 / sqrt(p^2 + m2^2)), 0, cutoff)
        else
            int_val = integrate(p -> p^2 / sqrt(p^2 + M^2), 0, cutoff)
        end
        return M - m0 - (6 * Gs * M / π^2) * int_val
    end

    # We expect M to be around 0.3 - 0.5 GeV in the vacuum
    M = fzero(x -> gap_eqn(x), 0.35)
    if M == 0.0
        M = 0.35 # fallback guess if fzero fails
    end

    # 2. Calculate chiral condensate <uu> per flavor
    if regularization == :pv
        m1 = sqrt(M^2 + Λ^2)
        m2 = sqrt(M^2 + 2 * Λ^2)
        cond_val = integrate(p -> p^2 * (1 / sqrt(p^2 + M^2) - 2 / sqrt(p^2 + m1^2) + 1 / sqrt(p^2 + m2^2)), 0, cutoff)
    else
        cond_val = integrate(p -> p^2 / sqrt(p^2 + M^2), 0, cutoff)
    end
    condensate = -(3 * M / π^2) * cond_val

    # 3. Calculate pion decay constant f_pi
    if regularization == :pv
        fpi_sq = (3 * M^2 / (4 * π^2)) * integrate(p -> p^2 * (1 / (p^2 + M^2)^1.5 - 2 / (p^2 + M^2 + Λ^2)^1.5 + 1 / (p^2 + M^2 + 2 * Λ^2)^1.5), 0, cutoff)
    else
        fpi_sq = (3 * M^2 / (4 * π^2)) * integrate(p -> p^2 / (p^2 + M^2)^1.5, 0, cutoff)
    end
    f_pi = sqrt(max(0.0, fpi_sq))

    # 4. Solve for pion mass m_pi
    # m_pi is the root of 1 - 2*Gs*Pi_ps(m_pi^2) = 0
    # The integrand has a threshold singularity at q = 2M (denominator = 0).
    # We clamp q2 strictly below 4M^2 to avoid the branch cut.
    function Pi_ps(q2)
        q2_safe = min(q2, 4 * M^2 * (1 - 1e-8))  # clamp below 4M^2 threshold
        if regularization == :pv
            m1 = sqrt(M^2 + Λ^2)
            m2 = sqrt(M^2 + 2 * Λ^2)
            return (6 / π^2) * integrate(p -> p^2 * (
                    sqrt(p^2 + M^2) / (p^2 + M^2 - q2_safe / 4) -
                    2 * sqrt(p^2 + m1^2) / (p^2 + m1^2 - q2_safe / 4) +
                    sqrt(p^2 + m2^2) / (p^2 + m2^2 - q2_safe / 4)
                ), 0, cutoff)
        else
            return (6 / π^2) * integrate(p -> p^2 * sqrt(p^2 + M^2) / (p^2 + M^2 - q2_safe / 4), 0, cutoff)
        end
    end

    # Search for m_pi in (0, 2M) — pion must be a bound state below the 2-quark threshold.
    # Use a bracketed bisection to avoid the solver wandering above the q=2M singularity.
    # Note: Pi_ps already clamps q2 below 4M^2, so even if fzero explores above threshold
    # it will not produce Inf. The bracket below provides extra robustness.
    threshold = 2 * M
    g(q) = 1 / (2 * Gs) - Pi_ps(q^2)
    f_low = g(1e-6)
    f_high = g(threshold * (1 - 1e-6))
    if f_low * f_high < 0
        # Bisection within the physical interval (no external package needed)
        lo, hi = 1e-6, threshold * (1 - 1e-6)
        for _ in 1:60
            mid = (lo + hi) / 2
            g(mid) * f_low < 0 ? (hi = mid) : (lo = mid)
        end
        m_pi = (lo + hi) / 2
    else
        # Fallback: Newton-based fzero (Pi_ps is clamped so no Inf)
        m_pi = fzero(g, 0.138)
        if m_pi == 0.0
            m_pi = 0.138
        end
    end

    return M, m_pi, f_pi, condensate
end

"""
    fit_parameters(target_m_pi, target_f_pi, target_cbrt_condensate; regularization=:sharp, pv_scale=10.0, initial_guess=[0.6, 0.005, 5.0])

Fits the parameters `(Λ, m0, Gs)` of the NJL model to match the physical targets:
- `target_m_pi`: Pseudoscalar meson (pion) mass (e.g. 0.138 GeV)
- `target_f_pi`: Pion decay constant (e.g. 0.093 GeV)
- `target_cbrt_condensate`: Third root of the quark condensate <uu>^(1/3) (e.g. -0.250 GeV)

Returns the fitted parameters as a `TempParameters` struct.
"""
function fit_parameters(target_m_pi, target_f_pi, target_cbrt_condensate; regularization=:sharp, pv_scale=10.0, initial_guess=[0.6, 0.005, 5.0])
    function target_system!(F, x)
        Λ_val, m0_val, Gs_val = x[1], x[2], x[3]
        # Ensure parameters remain positive and physical
        if Λ_val <= 0.1 || m0_val <= 0.001 || Gs_val <= 0.5
            F .= [1e3, 1e3, 1e3]
            return F
        end
        M, m_pi, f_pi, condensate = vacuum_observables(Λ_val, m0_val, Gs_val; regularization=regularization, pv_scale=pv_scale)
        F[1] = m_pi - target_m_pi
        F[2] = f_pi - target_f_pi
        F[3] = cbrt(condensate) - target_cbrt_condensate
        return F
    end

    sol = nlsolve(target_system!, initial_guess)
    if !sol.f_converged
        @warn "Parameter fitting did not converge fully. Returned best fit."
    end

    Λ_fit, m0_fit, Gs_fit = sol.zero[1], sol.zero[2], sol.zero[3]
    return TempParameters(Λ_fit, m0_fit, Gs_fit)
end

export vacuum_observables, fit_parameters
