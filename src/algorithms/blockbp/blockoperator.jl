

struct BlockOperator{M}
	data::SquareLatticeBonds{Union{M, Nothing}}
	partition::SquareLatticePartition


function BlockOperator{M}(data::SquareLatticeBonds{Union{M, Nothing}}, partition::SquareLatticePartition) where M
	@assert size(data) == size(partition)
	@assert check_v_block_operator(data.V, partition)
	@assert check_h_block_operator(data.H, partition)
	new{M}(data, partition)
end

end

BlockOperator(data::SquareLatticeBonds{Union{M, Nothing}}, partition::SquareLatticePartition) where M = BlockOperator{M}(data, partition)


Base.size(x::BlockOperator) = size(x.data)
nrows(x::BlockOperator) = nrows(x.partition)
ncols(x::BlockOperator) = ncols(x.partition)

n_nontrivial_terms(x::BlockOperator) = n_nontrivial_terms(x.data) 


function subblock(blk::BlockOperator{M}, i::Int, j::Int) where {M}
	rows = rowindices(blk.partition, i)
	cols = colindices(blk.partition, j)

	subU = _empty_square_lattice(M, length(rows), length(cols))
	for i in 1:size(subU,1)
		for j in 1:size(subU,2)-1
			subU.H[i, j] = blk.data.H[rows[i], cols[j]]
		end
	end
	for i in 1:size(subU,1)-1
		for j in 1:size(subU,2)
			subU.V[i, j] = blk.data.V[rows[i], cols[j]]
		end
	end
	return subU	
end

function subblocks(blk::BlockOperator{M}) where {M}
	r = Matrix{SquareLatticeBonds{Union{M, Nothing}}}(undef, nrows(blk), ncols(blk))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			r[i, j] = subblock(blk, i, j)
		end
	end
	return r
end


function _empty_square_lattice(::Type{M}, m::Int, n::Int) where {M}
	V = Matrix{Union{M, Nothing}}(nothing, m, n)
	H = Matrix{Union{M, Nothing}}(nothing, m, n)
	return SquareLatticeBonds(V, H)
end

function check_v_block_operator(V::PeriodicArray{Union{M, Nothing}, 2}, partition::SquareLatticePartition) where M
	m, n = size(partition)
	rows = Set(vbonds_rowcenterindices(partition))
	cols = Set(vbonds_colcenterindices(partition))
	for i in 1:m
		for j in 1:n
			if !((i in rows) && (j in cols))
				isnothing(V[i, j]) || return false
			end
		end
	end
	return true
end

function check_h_block_operator(H::PeriodicArray{Union{M, Nothing}, 2}, partition::SquareLatticePartition) where M
	m, n = size(partition)
	rows = Set(hbonds_rowcenterindices(partition))
	cols = Set(hbonds_colcenterindices(partition))
	for i in 1:m
		for j in 1:n
			if !((i in rows) && (j in cols))
				isnothing(H[i, j]) || return false
			end
		end
	end
	return true
end

