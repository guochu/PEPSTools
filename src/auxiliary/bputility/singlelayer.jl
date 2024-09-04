
# single layer algorithm
function sl_compute_out_messages_ascend!(msg_out::AbstractVector, t::AbstractArray{T, 1}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T<:Number}
	@assert length(msg_in) == length(msg_out)  == 1
	axpy!(true, StridedView(t), msg_out[1])
end
function sl_compute_out_messages_ascend!(msg_out::AbstractVector, t::AbstractArray{T, N}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T<:Number, N}
	@assert length(msg_in) == N
	t_size = size(t)
	t_size_tail = tail(t_size)
	t_size_tail_prod = prod(t_size_tail)
	@assert length(workspace) >= t_size_tail_prod
	tmp = reshape(view(workspace, 1:t_size_tail_prod), 1, t_size_tail_prod)

	t′ = reshape(mul!(tmp, transpose(msg_in[1]), reshape(t, t_size[1], t_size_tail_prod)), t_size_tail)
	return sl_compute_out_messages_ascend!(view(msg_out, 2:N), t′, view(msg_in, 2:N), view(workspace, (t_size_tail_prod+1):length(workspace)))
end

function _sl_compute_out_messages!(msg_out::AbstractVector, t::AbstractArray{T, 1}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T}
	@assert length(msg_in) == length(msg_out) == 1
	axpy!(true, StridedView(t), msg_out[1])
	return msg_out
end

function _sl_compute_out_messages!(msg_out::AbstractVector, t::AbstractArray{T, N}, msg_in::AbstractVector, workspace::AbstractVector{T}) where {T, N}
	@boundscheck begin
		(length(msg_in) == length(msg_out) == N) || throw(ArgumentError("message size mismatch with node rank"))
		for (i, (a, b)) in enumerate(zip(msg_in, msg_out))
			(length(a) == length(b) === size(t, i)) || throw(ArgumentError("$i-th message size mismatch"))
		end
	end
	sl_compute_out_messages_ascend!(msg_out, t, msg_in, workspace)
	t_size = size(t)
	t_size_front = front(t_size)
	t_size_front_prod = prod(t_size_front)
	@assert length(workspace) >= t_size_front_prod

	tmp = view(workspace, 1:t_size_front_prod)
	t′ = reshape(mul!(tmp, reshape(t, t_size_front_prod, t_size[N]), msg_in[N]), t_size_front)
	_sl_compute_out_messages!(view(msg_out, 1:N-1), t′, view(msg_in, 1:N-1), view(workspace, (t_size_front_prod+1):length(workspace) ))
	return msg_out
end
function sl_compute_out_messages(t::Array{T, N}, msg_in::Vector, workspace::Vector{T}=Vector{T}()) where {T, N}
	L = com_2_workspace(t)
	if length(workspace) < L 
		resize!(workspace, L)
	end
	return _sl_compute_out_messages!(zero.(msg_in), t, msg_in, workspace)
end

function sl_compute_out_messages_v1_ascend(t::AbstractArray{T, 1}, msg_in::AbstractVector) where {T}
	@assert length(msg_in) == 1
	return t
end

function sl_compute_out_messages_v1_ascend(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N}
	@assert length(msg_in) == N
	t_size = size(t)
	t_size_tail = tail(t_size)
	t_size_tail_prod = prod(t_size_tail)
	t2 = reshape(transpose(msg_in[1]) * reshape(t, t_size[1], t_size_tail_prod), t_size_tail)
	return sl_compute_out_messages_v1_ascend(t2, view(msg_in, 2:N))
end

function sl_compute_out_messages_v1_util(t::AbstractArray{T, 1}, msg_in::AbstractVector) where {T}
	@assert length(msg_in) == 1
	return [t]
end
function sl_compute_out_messages_v1_util(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N}
	@boundscheck begin
		(length(msg_in) == N) || throw(ArgumentError("message size mismatch with node rank"))
		for (i, a) in enumerate(msg_in)
			(length(a) == size(t, i)) || throw(ArgumentError("$i-th message size mismatch"))
		end
	end	
	msg_out = sl_compute_out_messages_v1_ascend(t, msg_in)

	t_size = size(t)
	t_size_front = front(t_size)
	t_size_front_prod = prod(t_size_front)

	t2 = reshape(reshape(t, t_size_front_prod, t_size[N]) * msg_in[N], t_size_front)
    return [sl_compute_out_messages_v1_util(t2, view(msg_in, 1:N-1))..., msg_out]
end
sl_compute_out_messages_v1(t::Array{T, N}, msg_in::Vector) where {T, N} = sl_compute_out_messages_v1_util(t, msg_in)

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


sl_compute_out_messages_v2(t::Array{T, N}, msg_in::Vector) where {T, N} = [sl_compute_out_message(t, msg_in, i) for i in 1:N]

function sl_compute_out_message(t::AbstractArray{T, N}, msg_in::AbstractVector, i::Int) where {T, N}
	@assert length(msg_in) == N
	@boundscheck begin
		for k in 1:N
			(size(t, k) == length(msg_in[k])) || throw(ArgumentError("dimension $k size mismatch"))
		end
	end
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

function sl_compute_out_message_v2(t::AbstractArray{T, N}, msg_in::AbstractVector, i::Int) where {T, N}
	@assert length(msg_in) == N
	out = zeros(T, size(t, i))
	msg_in′ = [m for (n, m) in enumerate(msg_in) if n != i]
	for j in 1:size(t, i)
		t′ = StridedView(selectdim(t, i, j))
		out[j] = sl_contract_node_v2(t′, msg_in′)
	end
	return out
end

function sl_compute_out_message_v3(t::Array{T, N}, msg_in::Vector, i::Int) where {T, N}
	@assert length(msg_in) == N
	@boundscheck begin
		for k in 1:N
			(size(t, k) == length(msg_in[k])) || throw(ArgumentError("dimension $k size mismatch"))
		end
	end
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


sl_contract_node(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N} = sl_contract_node_util(t, msg_in)

function sl_contract_node_util(t::AbstractArray{T, 1}, msg_in::AbstractVector) where {T}
	@assert length(msg_in) == 1
	return transpose(t) * msg_in[1]
end
function sl_contract_node_util(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N}
	t_size = size(t)
	t_size_front = front(t_size)
	t_size_front_prod = prod(t_size_front)
	t2 = reshape(reshape(t, t_size_front_prod, t_size[N]) * msg_in[N], t_size_front)
	return sl_contract_node_util(t2, view(msg_in, 1:N-1))
end


function sl_contract_node_v2(t::AbstractArray{T, N}, msg_in::Vector) where {T, N}
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


# sl_contract_node_v3(t::AbstractArray{T, N}, msg_in::AbstractVector) where {T, N} = sl_contract_node_v3_util(t, msg_in, one(T))
# function sl_contract_node_v3_util(t::AbstractArray{T, 1}, msg_in::AbstractVector, v::T) where {T}
# 	@assert length(msg_in) == 1
# 	# return v * sum(((x, y),)->x*y, zip(t, msg_in[1]))
# 	return v * (transpose(t) * msg_in[1])
# end

# function sl_contract_node_v3_util(t::AbstractArray{T, N}, msg_in::AbstractVector, v::T) where {T, N}
# 	@assert length(msg_in) == N
# 	out = zero(T)
# 	for j in 1:size(t, N)
# 		t′ = StridedView(selectdim(t, N, j))
# 		# v = v * msg_in[N][j]
# 		out += sl_contract_node_v3_util(t′, view(msg_in, 1:N-1), v * msg_in[N][j])
# 	end
# 	return out
# end




