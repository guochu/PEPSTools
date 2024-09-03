
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


Zygote.@adjoint compute_out_message_and_physical(a::Array{T, M}, msg_in::Vector, i::Int) where {T, M} = compute_out_message_and_physical(a, msg_in, i), z -> begin
	N = M-1
	d = size(a, M)
	a_back = zeros(T, size(a))

	for k1 in 1:d
		a′ = selectdim(a_back, M, k1)
		indices = CartesianIndices(a′)
		kron_indices = map((x, y)->LinearIndices((x, y)), size(a′), size(a′))
		for k2 in 1:d
			b′ = selectdim(a, M, k2)
			c′ = view(z, :, k1, k2)
			for (idxa, aj) in zip(indices, a′)
				for (idxb, bj) in zip(indices, b′)
					v1 = bj
					v2 = bj
					# v = conj(c′[kron_indices[i][idxa[i], idxb[i]]]) * bj
					for j in 1:i-1
						v1 *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
						v2 *= conj(msg_in[j][kron_indices[j][idxb[j], idxa[j]]])
					end
					for j in i+1:N
						v1 *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
						v2 *= conj(msg_in[j][kron_indices[j][idxb[j], idxa[j]]])
					end
					cj = c′[kron_indices[i][idxa[i], idxb[i]]]
					a′[idxa] += conj(cj) * v1 + cj * v2
				end
			end
		end				
	end
	return a_back, nothing, nothing
end	
