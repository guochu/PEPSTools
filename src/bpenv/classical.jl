

struct ClassicalBPEnv{T <: Number} <: AbstractBPEnvironment
	peps::SquareTN{T}
	messages::SquareLatticeBondMessages{T}
end


environments(state::SquareTN, nitr::Int, alg::BP) = environments(state, _init_c_bondmessages(state, alg.initguess, alg.seed), alg)

function environments(state::SquareTN, init_msgs::SquareLatticeBondMessages, alg::BP)
	(size(state) === size(init_msgs)) || throw(ArgumentError("graph mismatch"))
	
	msgs, converged = fixedpoint_messages(state, init_msgs, alg)
	return ClassicalBPEnv(state, msgs)
end