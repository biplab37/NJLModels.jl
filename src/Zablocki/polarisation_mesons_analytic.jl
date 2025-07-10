# This file contains code for the analytic structure of mesons

function polarisation_meson(T, mu, ω, q, gap, Nm, param)
    if imag(ω) == 0
        return realpart_meson_q(ω, T, mu, q, gap, Nm, param)
    end
    if imag(ω) > 0
        return nothing
    end

end
