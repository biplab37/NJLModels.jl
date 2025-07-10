abstract type Parameters end

Base.@kwdef mutable struct Parameters1 <: Parameters
    Λ = 0.5879
    m0 = 0.005588
    Gs = 2.442 / Λ^2
    Gv = Gs / 2
    eta_d = 0.75
    GD = eta_d * Gs
end

Base.@kwdef mutable struct Parameters2 <: Parameters
    Λ = 0.639
    m0 = 0.0055
    Gs = 2.134 / Λ^2
    Gv = Gs / 2
    eta_d = 0.75
    GD = eta_d * Gs
end

function Base.show(io::IO, ::MIME"text/plain", p::Parameters)
    return dump(p)
end

export Parameters
