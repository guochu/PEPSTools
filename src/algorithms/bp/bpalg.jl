
struct BP <: ImaginaryTimePEPSUpdateAlgorithm
	msg_tol::Float64 
	msg_maxiter::Int
	verbosity::Int 
end

BP(; msg_maxiter::Int=10, msg_tol::Real=1.0e-8, verbosity::Int=1) = BP(msg_tol, msg_maxiter, verbosity)