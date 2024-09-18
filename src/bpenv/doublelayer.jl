

struct DoubleLayerBPEnv{T <: Number} <: AbstractBPEnvironment
	peps::PEPS{T}
	messages::SquareLatticeBondMessages{T}
end


environments(state::PEPS, alg::BP) = environments(state, _init_q_bondmessages(state, alg.initguess, alg.seed), alg)

function environments(state::PEPS, init_msgs::SquareLatticeBondMessages, alg::BP)
	(size(state) === size(init_msgs)) || throw(ArgumentError("graph mismatch"))
	
	msgs, converged = fixedpoint_messages(state, init_msgs, alg)
	return DoubleLayerBPEnv(state, msgs)
end
