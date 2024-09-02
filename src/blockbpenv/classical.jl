


"""
    struct BlockBPPartitionSquareTN{T<:Number}

two dimensional tensor network.
The index convention for site tensor of 2d tn
------------------2-------------
--------------1-------3---------
------------------4-------------
"""
struct BlockBPPartitionSquareTN{T<:Number, _MESSAGE} <: AbstractBlockBPPartitionPEPS{T}
	peps::SquareTN{T}
	partition::SquareLatticePartition
	row_msgs::PeriodicArray{ _MESSAGE, 2}
	col_msgs::PeriodicArray{ _MESSAGE, 2}
end


function BlockBPPartitionSquareTN(peps::SquareTN{T}, partition::SquareLatticePartition) where {T <: Number}
	n_rows, n_cols = nrows(partition), ncols(partition)
	_MESSAGET = Message{MPS{T, real(T)}, MPS{T, real(T)}}
	row_msgs = PeriodicArray{ _MESSAGET, 2 }(undef, n_rows, n_cols)
	col_msgs = PeriodicArray{ _MESSAGET, 2 }(undef, n_rows, n_cols)
	for i in 1:n_rows
		for j in 1:n_cols
			_bk = subblock(peps.data, partition, i, j)
			row_msgs[i, j] = random_double_layer_mps_message(T, [size(item, 2) for item in _bk[1, :]])
			col_msgs[i, j] = random_double_layer_mps_message(T, [size(item, 1) for item in _bk[:, 1]])
		end
	end
	return BlockBPPartitionSquareTN(peps, partition, row_msgs, col_msgs)
end

peps_partition(peps::SquareTN, partition::SquareLatticePartition) = BlockBPPartitionSquareTN(peps, partition)


function random_double_layer_boundary_mps(::Type{T}, ds::AbstractVector{Int}; D::Int) where {T <: Number}
	L = length(ds)
	r = Vector{Array{T, 3}}(undef, L)
	for i in 1:L
		Dl = (i==1) ? 1 : D
		Dr = (i==L) ? 1 : D
		tmp = randn(T, Dl, ds[i], Dr) / D
		r[i] = tmp
	end
	return MPS(r)
end
random_double_layer_boundary_mps(::Type{T}, L::Int; D::Int) where {T <: Number} = random_double_layer_boundary_mps(T, [D for i in 1:L]; D=D)


function _random_double_layer_mps_message(::Type{T}, physpaces::Vector{Int}; Di::Int=maximum(physpaces), Do::Int=Di) where {T <: Number}
	i = random_double_layer_boundary_mps(T, physpaces, D=Di)
	# normalize!(i, iscanonical=false)
	rightorth!(i, alg=Orthogonalize(alg=QR(), normalize=true))

	o = random_double_layer_boundary_mps(T, physpaces, D=Do)
	# normalize!(o, iscanonical=false)
	rightorth!(o, alg=Orthogonalize(alg=QR(), normalize=true))
	return Message(i, o)
	
end

function random_double_layer_mps_message(::Type{T}, physpaces::Vector{Int}; Di::Int=maximum(physpaces), Do::Int=Di) where {T <: Number}
	_physpaces = [round(Int, sqrt(d)) for d in physpaces]
	if _physpaces.^2 == physpaces
		random_mps_message(T, _physpaces; Di=round(Int, sqrt(Di)), Do=round(Int, sqrt(Do)))
	else
		_random_double_layer_mps_message(T, physpaces, Di=Di, Do=Do)
	end
end 
