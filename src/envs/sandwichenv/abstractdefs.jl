abstract type AbstractSandwichEnv end


abstract type AbstractSingleLayerSandwichEnv <: AbstractSandwichEnv end


Base.length(x::AbstractSandwichEnv) = length(x.middle)


