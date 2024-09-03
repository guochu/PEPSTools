
Zygote.@adjoint energy(h::SquareLatticeOperator, state::PEPS, env::DoubleLayerBPEnv) = begin
	(size(h) === size(state) === size(env)) || throw(DimensionMismatch("size mismatch"))
	
	energy = zero(real(scalartype(state)))
	messages = env.messages

	# # onebody observables
	# for node in 1:length(state)
	# 	msgs_in = get_in_messages(messages, node)
	# 	energy += node_energy(state[node], msgs_in, h.onebody[node])
	# end

	# two body observables
	tn_new = [zero(item) for item in state.data]
	for _edge in edges(messages)
		src, dst = _edge
		i1 = findfirst(x->x==dst, neighbors(messages, src))
		i2 = findfirst(x->x==src, neighbors(messages, dst))
		msgs_src = get_in_messages(messages, src)
		msgs_dst = get_in_messages(messages, dst)
		hj = h[(_edge[1], _edge[2])]
		if !isnothing(hj)
			tmp, back = Zygote.pullback(bond_energy, state[src], state[dst], msgs_src, msgs_dst, i1=>i2, hj)
			state_src_back, state_dst_back = back(one(tmp))
			tn_new[src] .+= state_src_back
			tn_new[dst] .+= state_dst_back
			energy += real(tmp)
		end
	end

	return energy, z -> begin
		for item in tn_new
			LinearAlgebra.lmul!(z, item)
		end
		return nothing, PEPS(tn_new), nothing
	end
end

