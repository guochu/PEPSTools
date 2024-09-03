# variational BP
function energy(h::SquareLatticeOperator, state::PEPS, env::DoubleLayerBPEnv)
	(size(h) === size(state) === size(env)) || throw(DimensionMismatch("size mismatch"))
	
	energy = zero(scalartype(state))
	messages = env.messages

	# # onebody observables
	# for node in 1:length(state)
	# 	msgs_in = get_in_messages(messages, node)
	# 	energy += node_energy(state[node], msgs_in, h.onebody[node])
	# end

	# two body observables
	for _edge in edges(messages)
		src, dst = _edge
		i1 = findfirst(x->x==dst, neighbors(messages, src))
		i2 = findfirst(x->x==src, neighbors(messages, dst))
		msgs_src = get_in_messages(messages, src)
		msgs_dst = get_in_messages(messages, dst)
		hj = h[(_edge[1], _edge[2])]
		if !isnothing(hj)
			tmp = bond_energy(state[src], state[dst], msgs_src, msgs_dst, i1=>i2, hj)
			energy += tmp
		end
	end
	return energy
end

function expectation(h::SquareLatticeOperator, state::PEPS, env::DoubleLayerBPEnv)
	(size(h) === size(state) === size(env)) || throw(DimensionMismatch("size mismatch"))

	energies = SquareLatticeBonds(zeros(scalartype(state), size(state)), zeros(scalartype(state), size(state)))
	
	messages = env.messages

	# # onebody observables
	# for node in 1:length(state)
	# 	msgs_in = get_in_messages(messages, node)
	# 	energy += node_energy(state[node], msgs_in, h.onebody[node])
	# end

	# two body observables
	for _edge in edges(messages)
		src, dst = _edge
		i1 = findfirst(x->x==dst, neighbors(messages, src))
		i2 = findfirst(x->x==src, neighbors(messages, dst))
		msgs_src = get_in_messages(messages, src)
		msgs_dst = get_in_messages(messages, dst)
		__edge = (_edge[1], _edge[2])
		hj = h[__edge]
		if !isnothing(hj)
			tmp = bond_energy(state[src], state[dst], msgs_src, msgs_dst, i1=>i2, hj)
			energies[__edge] += tmp
		end
	end
	return energies
end


function node_energy(t::AbstractArray, msgs::Vector, h::AbstractMatrix)
	out = compute_physical(t, msgs)
	num = mapreduce(*, +, out, h)
	return num / tr(out)
end

function bond_energy(a::AbstractArray, b::AbstractArray, msgs_a::Vector, msgs_b::Vector, edge::Pair{Int, Int}, op::AbstractArray{<:Number, 4})
	neighbour1 = compute_out_message_and_physical(a, msgs_a, edge[1])
	neighbour2 = compute_out_message_and_physical(b, msgs_b, edge[2])
	@tensor num = neighbour1[1,2,3] * op[2,3,4,5] * neighbour2[1,4,5]
	@tensor den = neighbour1[1,2,2] * neighbour2[1,3,3]
	return num / den
end

