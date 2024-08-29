




"""
    struct BeliefPEPSBlock{T<:Number}

two dimensional tensor network.
The index convention for site tensor of 2d tn
------------------2-------------
--------------1-------3---------
------------------4-------------
"""
struct BeliefPEPSBlock{T<:Number, _MESSAGE} <: AbstractBlockBPPartitionPEPS{T}
	peps::PEPS{T}
	partition::SquareLatticePartition
	row_msgs::PeriodicArray{ _MESSAGE,2}
	col_msgs::PeriodicArray{ _MESSAGE,2}
end


function BeliefPEPSBlock(peps::PEPS{T}, partition::SquareLatticePartition) where {T <: Number}
	n_rows, n_cols = nrows(partition), ncols(partition)
	_MESSAGET = Message{MPS{T, real(T)}, MPS{T, real(T)}}
	row_msgs = PeriodicArray{ _MESSAGET, 2}(undef, n_rows, n_cols)
	col_msgs = PeriodicArray{ _MESSAGET, 2 }(undef, n_rows, n_cols)
	for i in 1:n_rows
		for j in 1:n_cols
			_bk = subblock(peps.data, partition, i, j)
			row_msgs[i, j] = random_mps_message(T, [size(item, 2) for item in _bk[1, :]])
			col_msgs[i, j] = random_mps_message(T, [size(item, 1) for item in _bk[:, 1]])
		end
	end
	return BeliefPEPSBlock(peps, partition, row_msgs, col_msgs)
end

peps_partition(peps::PEPS, partition::SquareLatticePartition) = BeliefPEPSBlock(peps, partition)