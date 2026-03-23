"""
Top level abstraction for controlling the constant parameters
"""
abstract type Parameters end

Base.@kwdef struct Parameters1 <: Parameters
    Λ::Float64 = 0.5879
    m0::Float64 = 0.005588
    Gs::Float64 = 2.442 / Λ^2
    Gv::Float64 = Gs / 2
    eta_d::Float64 = 0.75
    GD::Float64 = eta_d * Gs
end

Base.@kwdef struct Parameters2 <: Parameters
    Λ::Float64 = 0.639
    m0::Float64 = 0.0055
    Gs::Float64 = 2.134 / Λ^2
    Gv::Float64 = Gs / 2
    eta_d::Float64 = 0.75
    GD::Float64 = eta_d * Gs
end

# From Maslov paper
Base.@kwdef struct ParametersMaslov <: Parameters
    Λ::Float64 = 0.651
    m0::Float64 = 0.0055
    Gs::Float64 = 5.04
    Gv::Float64 = Gs / 2
    eta_d::Float64 = 0.75
    GD::Float64 = eta_d * Gs
end

function Base.show(io::IO, ::MIME"text/plain", p::Parameters)
    return dump(p)
end

export Parameters
