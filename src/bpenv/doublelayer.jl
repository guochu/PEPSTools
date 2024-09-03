

struct DoubleLayerBPEnv{T <: Number} <: AbstractBPEnvironment
	peps::PEPS{T}
	messages::SquareLatticeBondMessages{T}
end


bp_environments(state::PEPS, nitr::Int, tol::Real=-1; verbosity::Int=1) = bp_environments(state, unit_q_bondmessages(state), nitr, tol, verbosity=verbosity)

function bp_environments(state::PEPS, init_msgs::SquareLatticeBondMessages, nitr::Int, tol::Real=-1; verbosity::Int=1)
	(size(state) === size(init_msgs)) || throw(ArgumentError("graph mismatch"))
	
	msgs = fixedpoint_messages(state, init_msgs, FixedNorm(), nitr, tol, verbosity=verbosity)
	return DoubleLayerBPEnv(state, msgs)
end
