# Using interpolattion for the imaginary part, to calculate the real part
# Here we will specifically use the parts where the function is strictly zero
# I think the interpolation was bad previously on the region where the funtion
# was becoming finite from zero.

# typical region for imagpart being zero is 0 < ω^2 - q^2 < 4*m^2 and abs(ω) >cutoff

function interpolate_imag_part(func::Function, cutoff, threshold)
    
end
