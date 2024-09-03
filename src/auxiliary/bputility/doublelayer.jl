
dl_compute_out_messages(t::Array, msg_in::Vector) = [dl_compute_out_message(t, msg_in, i) for i in 1:length(msg_in)]


function dl_compute_out_message(t::AbstractArray{T, M}, msg_in::Vector, i::Int) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	a = t
	b = t
	d = size(a, M)
	out = zeros(T, length(msg_in[i]))

	for k in 1:d
		a′ = selectdim(a, M, k)
		b′ = selectdim(b, M, k)
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
			out[kron_indices[i][idxa[i], idxb[i]]] += v
		end				
	end
	return out
end


function dl_contract_node(t::AbstractArray{T, M}, msg_in::Vector) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	a = t
	b = t
	# index = ntuple(i->size(a, i), N)
	d = size(a, M)
	out = zero(T)
	for i in 1:d
		a′ = selectdim(a, M, i)
		b′ = selectdim(b, M, i)
		indices = CartesianIndices(a′)
		kron_indices = map((x, y)->LinearIndices((x, y)), size(a′), size(b′))
		for (idxa, aj) in zip(indices, a′), (idxb, bj) in zip(indices, b′)
			v = conj(aj) * bj
			for j in 1:N
				v *= msg_in[j][kron_indices[j][idxa[j], idxb[j]]]
			end
			out += v
		end
	end
	return out
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
