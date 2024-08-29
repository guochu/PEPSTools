abstract type AbstractSandwichEnv end


abstract type AbstractSquareTNSandwichEnv <: AbstractSandwichEnv end


Base.length(x::AbstractSandwichEnv) = length(x.middle)


