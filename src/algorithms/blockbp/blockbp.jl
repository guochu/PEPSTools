
abstract type AbstractBlockBPPEPSUpdateAlgorithm <: AbstractPEPSUpdateAlgorithm end

struct BlockBP{M<:BoundaryMPS} <: AbstractBlockBPPEPSUpdateAlgorithm
	block_size::Tuple{Int, Int} 
	msg_tol::Float64 
	msg_maxiter::Int
	msg_D::Int 
	verbosity::Int 
	update_alg::M 
end

BlockBP(; block_size::Tuple{Int, Int}=(3,3), msg_maxiter::Int=10, msg_tol::Real=1.0e-5, verbosity::Int=1, 
	update_alg::BoundaryMPS=BoundaryMPS(D2=3, D1=2*D2^2, verbosity=verbosity), msg_D=(update_alg.D2)^2) = BlockBP(
	block_size, msg_tol, msg_maxiter, msg_D, verbosity, update_alg)

# We may want to use different (smaller) bond dimension to compute messages to reduce the complexity
function get_msg_mult_alg(x::AbstractBlockBPPEPSUpdateAlgorithm)
	mult_alg = get_mult_alg(x.update_alg)
	return changeD(mult_alg, D=x.msg_D) 
end 


compute_messages!(blk::AbstractBlockBPPartitionPEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm) = compute_messages!(blk, 
	get_msg_mult_alg(alg), maxiter=alg.msg_maxiter, tol=alg.msg_tol, verbosity=alg.verbosity )
