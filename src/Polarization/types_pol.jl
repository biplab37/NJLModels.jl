
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

Base.@kwdef mutable struct MeanField_P <: MeanField
    pressure::Function
end

abstract type Gap end

Base.@kwdef mutable struct Gaps <: Gap
    variables::Vector
    gap_eqn::Function
    num_var::Int = length(variables)
    initial_guess::Vector
end

abstract type Model end

Base.@kwdef mutable struct CustomModel <: Model
    param::Parameters
    meanfield::MeanField
    polarization::PolarizationFunction
end
