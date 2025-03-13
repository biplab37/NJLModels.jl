# For the calculation of the Mott momenta we want the ressurection of the bound state. 
# In other word the phase shift should go to π at non-zero momenta for temperature
# greater than Mott temeperature. It is enough to check at the threshold s=4m^2 for 
# the bound state. The imaginary part of the polarization function at the threshold
# is zero so we want the real part to be negative. So to find the Mott momenta we 
# need to look for momenta at which the real part of the polarization function changes
# sign from positive to negative. i.e. to find the zero of realpart as a function of
# the momenta at a specific temperature.

function mott_momenta(T, mu, param::Parameters)
    gap = massgap(T, mu, param).zero
    threshold(q) = sqrt(q^2 + 4 * gap[1]^2)
    re_part(q) = realpart_meson_q(threshold(q), T, mu, q, gap, 0, param)

    return bisection(re_part, 0.0, 1.0)
end