# Auto differentiation utilities for single layer
Zygote.@adjoint SquareTN(data::AbstractMatrix{Array{T, 4}}) where {T<:Number} = SquareTN(data), z -> (z,)


Zygote.@adjoint sl_compute_out_message(t::Array{T, N}, msg_in::Vector, i::Int) where {T, N} = begin
	@assert length(msg_in) == N

	t_back = similar(t)
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
		t_back[index] = v
	end

	return out, z -> begin

		# z_conj = conj(z)
		
		# msg_in_back = [zeros(T, size(item)) for item in msg_in]
		msg_in_back = zero.(msg_in) 

		for (index, tj) in zip(CartesianIndices(t), t)
			v = t_back[index] * conj(z[index[i]])

			# assuming that there is no zeros...
			for j in 1:i-1
				msg_in_back[j][index[j]] += v / msg_in[j][index[j]]
			end	
			for j in i+1:N
				msg_in_back[j][index[j]] += v / msg_in[j][index[j]]
			end						

			t_back[index] = conj(v / tj)

		end

		conj!.(msg_in_back)

		return t_back, msg_in_back, nothing
		
	end	
end

Zygote.@adjoint sl_compute_out_messages(t::Array{T, N}, msg_in::Vector, workspace::Vector{T}) where {T, N} = sl_compute_out_messages(t, msg_in, workspace), z_msg -> begin
	# L = com_2_workspace(t)
	# if length(workspace) < L 
	# 	resize!(workspace, L)
	# end
	t_back, msg_back = sl_compute_out_messages_bp(t, msg_in, z_msg, workspace)
	return t_back, msg_back, nothing
end

function sl_compute_out_messages_bp(t::AbstractArray{T, N}, msg_in::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T, N}
	@boundscheck begin
		(length(msg_in) == length(z_msg) == N) || throw(ArgumentError("message size mismatch with node rank"))
		for (i, (a, b)) in enumerate(zip(msg_in, z_msg))
			(length(a) == length(b) === size(t, i)) || throw(ArgumentError("$i-th message size mismatch"))
		end
	end

	t_back = zero(t)
	msg_back = zero.(msg_in)
	conj!.(msg_in)
	sl_compute_out_messages_bp_node!(t_back, msg_in, z_msg)
	conj!.(msg_in)
	sl_compute_out_messages_bp_edge!(msg_back, t, msg_in, conj.(z_msg), workspace)
	return t_back, conj!.(msg_back)
end


function sl_compute_out_messages_bp_node_n!(t_back::AbstractArray{T, 1}, msg_in_conj::AbstractVector, z_msg::AbstractVector, v::T) where {T}
	axpy!(v, z_msg[1], t_back)
end
function sl_compute_out_messages_bp_node_z!(t_back::AbstractArray{T, 1}, msg_in_conj::AbstractVector, z_msg::AbstractVector, v::T) where {T}
	axpy!(v, msg_in_conj[1], t_back)
end

function sl_compute_out_messages_bp_node_n!(t_back::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector, v::T) where {T, N}
	for j in 1:size(t_back, N)
		t′ = selectdim(t_back, N, j)
		msg_in_conj′ = view(msg_in_conj, 1:N-1)
		z_msg′ = view(z_msg, 1:N-1)
		sl_compute_out_messages_bp_node_z!(t′, msg_in_conj′, z_msg′, v * z_msg[N][j])
		sl_compute_out_messages_bp_node_n!(t′, msg_in_conj′, z_msg′, v * msg_in_conj[N][j])
	end
end
function sl_compute_out_messages_bp_node_z!(t_back::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector, v::T) where {T, N}
	for j in 1:size(t_back, N)
		t′ = selectdim(t_back, N, j)
		msg_in_conj′ = view(msg_in_conj, 1:N-1)
		z_msg′ = view(z_msg, 1:N-1)
		sl_compute_out_messages_bp_node_z!(t′, msg_in_conj′, z_msg′, v * msg_in_conj[N][j])
	end
end
sl_compute_out_messages_bp_node!(t_back::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector) where {T, N} = sl_compute_out_messages_bp_node_n!(
								t_back, msg_in_conj, z_msg, one(T))

function sl_compute_out_messages_bp_edge_ascend_z!(msg_back::AbstractVector, t::AbstractArray{T, 1}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T}
	@assert length(msg_back) == length(msg_in_conj) == length(z_msg) == 1
	axpy!(true, t, msg_back[1])
end
function sl_compute_out_messages_bp_edge_ascend_n!(msg_back::AbstractVector, t::AbstractArray{T, 1}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T}
end

# function compute_out_messages_2_bp_edge_ascend_n!(msg_back::AbstractVector, t::AbstractArray{T, 2}, msg_in_conj::AbstractVector, z_msg::AbstractVector) where {T}
# 	mul!(msg_back[1], transpose(t), z_msg[2])
# 	mul!(msg_back[2], t, z_msg[1])
# end

function sl_compute_out_messages_bp_edge_ascend_n!(msg_back::AbstractVector, t::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T, N}
	t_size = size(t)
	t_size_tail = tail(t_size)
	t_size_tail_prod = prod(t_size_tail)
	tt = reshape(t, t_size[1], t_size_tail_prod)
	@assert length(workspace) >= t_size_tail_prod
	msg_back′ = view(msg_back, 2:N)
	msg_in_conj′ = view(msg_in_conj, 2:N)
	z_msg′ = view(z_msg, 2:N)
	if N > 1
		t1 = reshape(view(workspace, 1:t_size_tail_prod), 1, t_size_tail_prod)
		t1 = reshape(mul!(t1, transpose(msg_in_conj[1]), tt), t_size_tail)
		sl_compute_out_messages_bp_edge_ascend_n!(msg_back′, t1, msg_in_conj′, z_msg′, view(workspace, (t_size_tail_prod+1):length(workspace)))		
	end
	t2 = reshape(view(workspace, 1:t_size_tail_prod), 1, t_size_tail_prod)
	t2 = reshape(mul!(t2, transpose(z_msg[1]), tt), t_size_tail)
	sl_compute_out_messages_bp_edge_ascend_z!(msg_back′, t2, msg_in_conj′, z_msg′, view(workspace, (t_size_tail_prod+1):length(workspace)))
end
function sl_compute_out_messages_bp_edge_ascend_z!(msg_back::AbstractVector, t::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T, N}
	t_size = size(t)
	t_size_tail = tail(t_size)
	t_size_tail_prod = prod(t_size_tail)
	tt = reshape(t, t_size[1], t_size_tail_prod)
	@assert length(workspace) >= t_size_tail_prod
	msg_back′ = view(msg_back, 2:N)
	msg_in_conj′ = view(msg_in_conj, 2:N)
	z_msg′ = view(z_msg, 2:N)
	t2 = reshape(view(workspace, 1:t_size_tail_prod), 1, t_size_tail_prod)
	t2 = reshape(mul!(t2, transpose(msg_in_conj[1]), tt), t_size_tail)
	sl_compute_out_messages_bp_edge_ascend_z!(msg_back′, t2, msg_in_conj′, z_msg′, view(workspace, (t_size_tail_prod+1):length(workspace)))
end

function sl_compute_out_messages_bp_edge!(msg_back::AbstractVector, t::AbstractArray{T, 1}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T}
	@assert length(msg_back) == length(msg_in_conj) == length(z_msg) == 1
	return
end

function sl_compute_out_messages_bp_edge!(msg_back::AbstractVector, t::AbstractArray{T, N}, msg_in_conj::AbstractVector, z_msg::AbstractVector, workspace::AbstractVector) where {T, N}
	sl_compute_out_messages_bp_edge_ascend_n!(msg_back, t, msg_in_conj, z_msg, workspace)
	t_size = size(t)
	t_size_front = front(t_size)
	t_size_front_prod = prod(t_size_front)
	t1 = view(workspace, 1:t_size_front_prod)
	t1 = reshape(mul!(t1, reshape(t, t_size_front_prod, t_size[end]), z_msg[end]), t_size_front)
	msg_back′ = view(msg_back, 1:N-1)
	msg_in_conj′ = view(msg_in_conj, 1:N-1)
	z_msg′ = view(z_msg, 1:N-1)
	_sl_compute_out_messages!(msg_back′, t1, msg_in_conj′, view(workspace, (t_size_front_prod+1):length(workspace)))

	t′ = view(workspace, 1:t_size_front_prod)
	t′ = reshape(mul!(t′, reshape(t, t_size_front_prod, t_size[end]), msg_in_conj[end]), t_size_front)
	sl_compute_out_messages_bp_edge!(msg_back′, t′, msg_in_conj′, z_msg′, view(workspace, (t_size_front_prod+1):length(workspace)))
end


Zygote.@adjoint sl_contract_node(t::AbstractArray{T, N}, msg_in::Vector) where {T, N} = begin
	@assert length(msg_in) == N
	out = zero(T)
	t_back = similar(t, size(t))
	for (index, tj) in zip(CartesianIndices(t), t)
		v = tj
		for j in 1:N
			v *= msg_in[j][index[j]]
		end
		out += v
		t_back[index] = v
	end
	return out, z -> begin
		msg_in_back = [zero(item) for item in msg_in]
		conj_z = conj(z)
		for (index, tj) in zip(CartesianIndices(t), t)
			v = conj_z * t_back[index]
			for j in 1:N
				msg_in_back[j][index[j]] += v / msg_in[j][index[j]]
			end	
			t_back[index] = conj(v / tj)
		end
		conj!.(msg_in_back)
		return t_back, msg_in_back
	end
end