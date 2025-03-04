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
    println(io, "Parameters:")
    println(io, "Λ = ", p.Λ)
    println(io, "m0 = ", p.m0)
    println(io, "Gs = ", p.Gs)
    println(io, "Gv = ", p.Gv)
    println(io, "eta_d = ", p.eta_d)
    println(io, "GD = Gs*eta_d = ", p.GD)
end

export Parameters
