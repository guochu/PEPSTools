

struct ClassicalBPEnv{T <: Number} <: AbstractBPEnvironment
	peps::SquareTN{T}
	messages::SquareLatticeBondMessages{T}
end


bp_environments(state::SquareTN, nitr::Int, tol::Real=-1; verbosity::Int=1) = bp_environments(state, unit_c_bondmessages(state), nitr, tol, verbosity=verbosity)

function bp_environments(state::SquareTN, init_msgs::SquareLatticeBondMessages, nitr::Int, tol::Real=-1; verbosity::Int=1)
	(size(state) === size(init_msgs)) || throw(ArgumentError("graph mismatch"))
	
	msgs = fixedpoint_messages(state, init_msgs, FixedNorm(), nitr, tol, verbosity=verbosity)
	return ClassicalBPEnv(state, msgs)
end