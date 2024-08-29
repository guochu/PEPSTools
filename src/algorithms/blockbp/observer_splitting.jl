
struct BlockLocalOperator{M}
	data::SquareLatticeSites{M}
	partition::SquareLatticePartition


function BlockLocalOperator{M}(data::SquareLatticeSites{M}, partition::SquareLatticePartition) where M
	@assert size(data) == size(partition)
	@assert check_local_observers(data, partition)
	new{M}(data, partition)
end

end

BlockLocalOperator(data::SquareLatticeSites{M}, partition::SquareLatticePartition) where M = BlockLocalOperator{M}(data, partition)

Base.size(x::BlockLocalOperator) = size(x.data)
nrows(x::BlockLocalOperator) = nrows(x.partition)
ncols(x::BlockLocalOperator) = ncols(x.partition)

n_nontrivial_terms(x::BlockLocalOperator) = n_nontrivial_terms(x.data) 

function subblock(blk::BlockLocalOperator, i::Int, j::Int) 
	@assert (i <= nrows(blk)) && (j <= ncols(blk))
	return SquareLatticeSites(PeriodicArray(blk.data.data[rowindices(blk.partition, i), colindices(blk.partition, j)]))
end

function subblocks(blk::BlockLocalOperator{M}) where {M}
	r = Matrix{SquareLatticeSites{M}}(undef, nrows(blk), ncols(blk))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			r[i, j] = subblock(blk, i, j)
		end
	end
	return r
end

function operator_splitting(U::SquareLatticeSites{M}, pts::Vector{SquareLatticePartition}) where M
	isempty(pts) && error("no partition.")
	# U = default_embeding(U_ori, size(pts[1]))
	m, n = size(pts[1])
	Htable = get_table(U.data)

	blks = BlockLocalOperator{M}[]
	for pt in pts
		Uk = _empty_local_observables(M, m, n)
		fill_observer!(Uk, U, pt, Htable)
		if n_nontrivial_terms(Uk) > 0
			push!(blks, BlockLocalOperator(Uk, pt))
		end
	end

	all(Htable) || error("incompatible splitting.")

	check_observer_splitting(U, blks) || error("incompatible splitting.")
	return blks
end

function fill_observer!(Uk, U, blk, Htable)
	s1, s2 = shifted_starting(size(Uk), size(U))
	rows = hbonds_rowcenterindices(blk)
	cols = vbonds_colcenterindices(blk)

	for i in rows
		for j in cols
			if !(Htable[i-s1, j-s2])
				Uk[i, j] = U[i-s1, j-s2]
				Htable[i-s1, j-s2] = true
			end
		end
	end
end

center_splitting(U::SquareLatticeSites, block_size::Tuple{Int, Int}; center::Tuple{Int, Int}=(1,1)) = operator_splitting(
	U, center_splitting(size(U), block_size, center))

function _empty_local_observables(::Type{M}, m::Int, n::Int) where M
	data = Matrix{Union{M, Nothing}}(nothing, m, n)
	return SquareLatticeSites(PeriodicArray(data))
end

function check_observer_splitting(x::SquareLatticeSites{M}, blks::Vector{BlockLocalOperator{M}}) where {M}
	# return n_nontrivial_terms(x) == sum(n_nontrivial_terms.(blks)) 
	return n_nontrivial_terms(x) == sum(n_nontrivial_terms.(blks)) == sum([sum(n_nontrivial_terms.(subblocks(blk))) for blk in blks])
end


function check_local_observers(V::SquareLatticeSites, partition::SquareLatticePartition) 
	m, n = size(partition)
	rows = Set(hbonds_rowcenterindices(partition))
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