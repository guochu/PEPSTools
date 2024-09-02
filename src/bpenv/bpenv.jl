abstract type AbstractBPEnvironment end

include("bondmessages.jl")
include("updatemessages.jl")



# struct ExpectationEnvironment{T <: Number}
# 	tn::Vector{Tuple{Array{T, N}, Array{T, N}} where N}
# 	messages::Dict{Pair{Int, Int}, Vector{T}}
# 	graph::SimpleGraph
# end

# environments(state::GraphState, nitr::Int, tol::Real=-1; verbosity::Int=1) = environments(state, unit_qmessage(state), nitr, tol, verbosity=verbosity)

# function environments(state::GraphState, init_msgs::GraphMessage, nitr::Int, tol::Real=-1; verbosity::Int=1)
# 	(state.graph === init_msgs.graph) || throw(ArgumentError("graph mismatch"))
# 	tn = dot_tn(state, state)
# 	if tol > 0
# 		msgs = fixedpoint_messages(tn, init_msgs, FixedNorm(), nitr, tol, verbosity=verbosity)
# 	else
# 		msgs = fixedpoint_messages(tn, init_msgs, FixedNorm(), nitr, verbosity=verbosity)
# 	end
	
# 	return ExpectationEnvironment(tn.data, msgs.data, msgs.graph)
# end

# graphmessage(env::ExpectationEnvironment) = GraphMessage(env.messages, env.graph)