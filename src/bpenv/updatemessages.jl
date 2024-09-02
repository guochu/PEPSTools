

# function fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm, nitr::Int; verbosity::Int=1)
# 	err = 1.
# 	for i in 1:nitr
# 		msg_2 = update_messages(tn, msg)
# 		msg_2 = normalize!(msg_2, alg)
# 		err = average_distance2(msg_2, msg)
# 		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
# 		msg = msg_2
# 	end
# 	(verbosity > 0) && println("final BP error is ", err)
# 	return msg
# end

function fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm, nitr::Int, tol::Real; verbosity::Int=1)
	err = 1.
	i = 0
	while i < nitr
		msg_2 = update_messages(tn, msg)
		msg_2 = normalize!(msg_2, alg)
		err = distance2(msg_2, msg)
		msg = msg_2
		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
		i += 1
		if err < tol
			(verbosity > 0) && println("BP converges in $i iterations, error=", err)
			return msg
		end
	end
	if (verbosity > 0)
		(i == nitr) && println("BP fails to converge in $(nitr) iterations, error=", err)
	end
	return msg
end


function update_messages(tn::SingleLayerTN, msg::SquareLatticeBondMessages)
	msg_data_new = typeof(msg.data)()
	workspace = Vector{scalartype(msg)}()
	for node in 1:length(tn.data)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		# msg_out = compute_out_messages(aj_data, msg_in)
		msg_out = compute_out_messages_2(aj_data, msg_in, workspace)
		for (j, n) in enumerate(neighbors(tn.graph, node))
			msg_data_new[node=>n] = msg_out[j]
		end
	end
	return GraphMessage(msg.graph, msg_data_new)
end
function update_messages(tn::DoubleLayerTN, msg::SquareLatticeBondMessages)
	msg_data_new = typeof(msg.data)()
	for node in 1:length(tn.data)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		msg_out = compute_out_messages(aj_data, msg_in)
		for (j, n) in enumerate(neighbors(tn.graph, node))
			msg_data_new[node=>n] = msg_out[j]
		end
	end
	return GraphMessage(msg.graph, msg_data_new)
end

get_in_messages(msg::GraphMessage, node::Int) = [msg[n=>node] for n in neighbors(msg.graph, node)]


function compute_out_messages_ascend!(msg_out::AbstractVector, t::AbstractArray{T, 1}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T<:Number}
	@assert length(msg_in) == length(msg_out)  == 1
	axpy!(true, StridedView(t), msg_out[1])
end
function compute_out_messages_ascend!(msg_out::AbstractVector, t::AbstractArray{T, N}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T<:Number, N}
	@assert length(msg_in) == N
	t_size = size(t)
	t_size_tail = tail(t_size)
	t_size_tail_prod = prod(t_size_tail)
	@assert length(workspace) >= t_size_tail_prod
	tmp = reshape(view(workspace, 1:t_size_tail_prod), 1, t_size_tail_prod)

	t′ = reshape(mul!(tmp, transpose(msg_in[1]), reshape(t, t_size[1], t_size_tail_prod)), t_size_tail)
	return compute_out_messages_ascend!(view(msg_out, 2:N), t′, view(msg_in, 2:N), view(workspace, (t_size_tail_prod+1):length(workspace)))
end

function _compute_out_messages_2!(msg_out::AbstractVector, t::AbstractArray{T, 1}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T}
	@assert length(msg_in) == length(msg_out) == 1
	axpy!(true, StridedView(t), msg_out[1])
	return msg_out
end

function _compute_out_messages_2!(msg_out::AbstractVector, t::AbstractArray{T, N}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T, N}
	@boundscheck begin
		(length(msg_in) == length(msg_out) == N) || throw(ArgumentError("message size mismatch with node rank"))
		for (i, (a, b)) in enumerate(zip(msg_in, msg_out))
			(length(a) == length(b) === size(t, i)) || throw(ArgumentError("$i-th message size mismatch"))
		end
	end
	compute_out_messages_ascend!(msg_out, t, msg_in, workspace)
	t_size = size(t)
	t_size_front = front(t_size)
	t_size_front_prod = prod(t_size_front)
	@assert length(workspace) >= t_size_front_prod

	tmp = view(workspace, 1:t_size_front_prod)
	t′ = reshape(mul!(tmp, reshape(t, t_size_front_prod, t_size[N]), msg_in[N]), t_size_front)
	_compute_out_messages_2!(view(msg_out, 1:N-1), t′, view(msg_in, 1:N-1), view(workspace, (t_size_front_prod+1):length(workspace) ))
	return msg_out
end
function compute_out_messages_2(t::Array{T, N}, msg_in::Vector, workspace::Vector{T}=Vector{T}()) where {T, N}
	L = com_2_workspace(t)
	if length(workspace) < L 
		resize!(workspace, L)
	end
	return _compute_out_messages_2!(zero.(msg_in), t, msg_in, workspace)
end

function com_2_workspace(t::AbstractArray{T, N}) where {T, N}
	t_size = size(t)
	L1 = t_size[end]
	tmp = L1
	for j in N-1:-1:2
		tmp *= t_size[j] 
		L1 += tmp
	end
	L2 = t_size[1]
	tmp = L2
	for j in 2:N-1
		tmp *= t_size[j]
		L2 += tmp
	end
	return max(L1, L2)
end


compute_out_messages(t::Tuple, msg_in::Vector) = [compute_out_message(t, msg_in, i) for i in 1:length(msg_in)]
compute_out_messages(t::Array{T, N}, msg_in::Vector) where {T, N} = [compute_out_message(t, msg_in, i) for i in 1:N]

function compute_out_message_debug(t::Array{T, N}, msg_in::Vector, i::Int) where {T, N}
	@assert length(msg_in) == N
	out = Zygote.Buffer(t, size(t, i))
	for i in 1:length(out)
		out[i] = 0
	end
	for (index, tj) in zip(CartesianIndices(t), t)
		v = tj
		for j in 1:i-1
			v *= msg_in[j][index[j]]
		end
		for j in i+1:N
			v *= msg_in[j][index[j]]
		end
		out[index[i]] += v
	end
	return copy(out)
end

function compute_out_message(t::AbstractArray{T, N}, msg_in::AbstractVector, i::Int) where {T, N}
	@assert length(msg_in) == N
	out = zeros(T, size(t, i))
	for (index, tj) in zip(CartesianIndices(t), t)
		v = tj
		for j in 1:i-1
			v *= msg_in[j][index[j]]
		end
		for j in i+1:N
			v *= msg_in[j][index[j]]
		end
		out[index[i]] += v
	end
	return out
end

function compute_out_message_2(t::AbstractArray{T, N}, msg_in::AbstractVector, i::Int) where {T, N}
	@assert length(msg_in) == N
	out = zeros(T, size(t, i))
	msg_in′ = [m for (n, m) in enumerate(msg_in) if n != i]
	for j in 1:size(t, i)
		t′ = StridedView(selectdim(t, i, j))
		out[j] = contract_node_2(t′, msg_in′)
	end
	return out
end

function compute_out_message(t::Tuple{<:AbstractArray{T, M}, <:AbstractArray{T, M}}, msg_in::Vector, i::Int) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	a, b = t
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

function contract_node(t::AbstractArray{T, N}, msg_in::Vector) where {T, N}
	@assert length(msg_in) == N
	out = zero(T)
	for (index, tj) in zip(CartesianIndices(t), t)
		v = tj
		for j in 1:N
			v *= msg_in[j][index[j]]
		end
		out += v
	end
	return out
end

function contract_node_2_util(t::AbstractArray{T, 1}, msg_in::AbstractVector, v::T) where {T}
	@assert length(msg_in) == 1
	return v * sum(((x, y),)->x*y, zip(t, msg_in[1]))
end

function contract_node_2_util(t::AbstractArray{T, N}, msg_in::AbstractVector, v::T) where {T, N}
	@assert length(msg_in) == N
	out = zero(T)
	for j in 1:size(t, N)
		t′ = StridedView(selectdim(t, N, j))
		# v = v * msg_in[N][j]
		out += contract_node_2_util(t′, view(msg_in, 1:N-1), v * msg_in[N][j])
	end
	return out
end
contract_node_2(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N} = contract_node_2_util(t, msg_in, one(T))

function contract_node(t::Tuple{<:AbstractArray{T, M}, <:AbstractArray{T, M}}, msg_in::Vector) where {T, M}
	N = M-1
	@assert length(msg_in) == N
	a, b = t
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

