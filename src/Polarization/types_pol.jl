
abstract type Factor end

Base.@kwdef mutable struct Factor_Pi <: Factor
    prefactor::Float64
    factor::Function = (ω, q, m) -> (ω^2 - q^2 - 4 * m^2)
end

abstract type Denominator end

Base.@kwdef mutable struct Denominator1
    sign::Int = +1
    diffmass::Bool = false
end

Base.@kwdef mutable struct Denominator2
    sign::Int = +1
    diffmass::Bool = false
end

abstract type Numerator end

Base.@kwdef mutable struct Numerator_Pi_1
    numerator::Function

end

abstract type PolarizationFunction end

Base.@kwdef mutable struct CustomPolarization <: PolarizationFunction
    factor::Factor
    numerator::Numerator
    denominator::Denominator
end

abstract type MeanField end

abstract type Model end

Base.@kwdef mutable struct CustomModel <: Model
    param::Parameters
    meanfield::MeanField
    polarization::PolarizationFunction
end
