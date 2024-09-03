
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
