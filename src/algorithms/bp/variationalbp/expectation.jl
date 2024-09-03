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


function compute_physical(a::AbstractArray{T, M}, msg_in::Vector) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	# index = ntuple(i->size(a, i), N)
	d = size(a, M)
	out = zeros(T, d, d)
	for k1 in 1:d, k2 in 1:d
		a′ = selectdim(a, M, k1)
		b′ = selectdim(a, M, k2)
		indices = CartesianIndices(a′)
		kron_indices = map((x, y)->LinearIndices((x, y)), size(a′), size(b′))
		for (idxa, aj) in zip(indices, a′), (idxb, bj) in zip(indices, b′)
			v = conj(aj) * bj
			for j in 1:N
				v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
			end
			out[k1, k2] += v
		end
	end
	return out
end

# contract all the dimensions except physical and i-th virtual
function compute_out_message_and_physical_v2(a::Array{T, M}, msg_in::Vector, i::Int) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	d = size(a, M)
	# out = zeros(T, length(msg_in[i]), d, d)
	out = Zygote.Buffer(a, length(msg_in[i]), d, d)
	for i in 1:length(out)
		out[i] = 0
	end

	for k1 in 1:d, k2 in 1:d
		a′ = selectdim(a, M, k1)
		b′ = selectdim(a, M, k2)
		indices = CartesianIndices(a′)
		kron_indices = map((x, y)->LinearIndices((x, y)), size(a′), size(b′))
		for (idxa, aj) in zip(indices, a′), (idxb, bj) in zip(indices, b′)
			v = conj(aj) * bj
			for j in 1:i-1
				v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
			end
			for j in i+1:N
				v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
			end
			out[kron_indices[i][idxa[i], idxb[i]], k1, k2] += v
		end				
	end
	return copy(out)
end

# slightly more efficient than the above implementation
function compute_out_message_and_physical(a::Array{T, M}, msg_in::Vector, i::Int) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	d = size(a, M)
	out = zeros(T, length(msg_in[i]), d, d)

	for k2 in 1:d 
		b′ = selectdim(a, M, k2)
		indices = CartesianIndices(b′)
		kron_indices = map((x, y)->LinearIndices((x, y)), size(b′), size(b′))
		for k1 in 1:d 
			a′ = selectdim(a, M, k1)
			out′ = view(out, :, k1, k2)
			for (idxb, bj) in zip(indices, b′)
				for (idxa, aj) in zip(indices, a′)
					v = conj(aj) * bj
					for j in 1:i-1
						v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
					end
					for j in i+1:N
						v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
					end
					out′[kron_indices[i][idxa[i], idxb[i]]] += v
				end	
			end
		end			
	end
	return out
end
