
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