module Zablocki

using UsefulFunctions, NLsolve, QuadGK

include("../custom_types.jl")
include("custom_types_specific.jl")
include("../helper_functions.jl")
include("gapeqn.jl")
include("polarization_function.jl")
include("finite_momentum.jl")
include("mott_momenta.jl")
include("mean_field.jl")
include("fluctuations.jl")
include("diquarks.jl")

end
