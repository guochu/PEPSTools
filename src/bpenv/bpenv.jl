abstract type AbstractBPEnvironment end

Base.size(m::AbstractBPEnvironment) = size(m.peps)


struct BP{M <: MessageNormalizationAlgorithm} <: ImaginaryTimePEPSUpdateAlgorithm
	normalize_alg::M
	msg_tol::Float64 
	msg_maxiter::Int
	damping::Float64
	initguess::Symbol
	seed::Int
	verbosity::Int 
end

function BP(normalize_alg::MessageNormalizationAlgorithm; msg_maxiter::Int=100, msg_tol::Real=1.0e-8, 
														  initguess::Symbol=:unit, damping::Real=0.2, 
														  seed::Int=1354, verbosity::Int=0)
	mixing = convert(Float64, damping)
	(0 <= mixing <= 1) || throw(ArgumentError("damping should be in [0, 1]"))
	return BP(normalize_alg, convert(Float64, msg_tol), msg_maxiter, mixing, initguess, seed, verbosity)	
end
BP(; normalize_alg::MessageNormalizationAlgorithm=FixedNorm(), kwargs...) = BP(normalize_alg; kwargs...)


include("graphinterface.jl")
include("messageinitializers.jl")
include("normalizemessage.jl")
include("updatemessages.jl")
include("doublelayer.jl")
include("classical.jl")