
struct BlockLocalOperator{M}
	data::LocalObservers{M}
	partition::SquareLatticePartition


function BlockLocalOperator{M}(data::LocalObservers{M}, partition::SquareLatticePartition) where M
	@assert size(data) == size(partition)
	@assert check_local_observers(data, partition)
	new{M}(data, partition)
end

end

BlockLocalOperator(data::LocalObservers{M}, partition::SquareLatticePartition) where M = BlockLocalOperator{M}(data, partition)

Base.size(x::BlockLocalOperator) = size(x.data)
nrows(x::BlockLocalOperator) = nrows(x.partition)
ncols(x::BlockLocalOperator) = ncols(x.partition)

n_nontrivial_terms(x::BlockLocalOperator) = n_nontrivial_terms(x.data) 

function subblock(blk::BlockLocalOperator, i::Int, j::Int) 
	@assert (i <= nrows(blk)) && (j <= ncols(blk))
	return LocalObservers(PeriodicArray(blk.data.data[rowindices(blk.partition, i), colindices(blk.partition, j)]))
end

function subblocks(blk::BlockLocalOperator{M}) where {M}
	r = Matrix{LocalObservers{M}}(undef, nrows(blk), ncols(blk))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			r[i, j] = subblock(blk, i, j)
		end
	end
	return r
end

# default_embeding(x::LocalObservers, shape::Tuple{Int, Int}) = LocalObservers(default_embeding(x.data, shape))

# function default_splitting(x_ori::LocalObservers{M}, block_size::Tuple{Int, Int}) where {M}
# 	odd_partition = block_partition(size(x_ori), block_size, (1,1))
# 	x = default_embeding(x_ori, size(odd_partition))

# 	m, n = size(x)
# 	even = _empty_local_observables(M, m, n)
# 	odd = _empty_local_observables(M, m, n)

# 	even_partition = VHshift(odd_partition, 1, 1)
# 	for j in 1:n
# 		col_pos = findfirst(k -> k == j, odd_partition.cols)
# 		if isnothing(col_pos)
# 			for i in 1:m
# 				odd.data[i, j] = x.data[i, j]
# 			end
# 		else
# 			for i in 1:m
# 				even.data[i, j] = x.data[i, j]
# 			end
# 		end
# 	end


# 	blks = [BlockLocalOperator(even, even_partition), BlockLocalOperator(odd, odd_partition)]

# 	check_observer_splitting(x, blks) || error("incompatible splitting.")
# 	return blks
# end

function operator_splitting(U::LocalObservers{M}, pts::Vector{SquareLatticePartition}) where M
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

center_splitting(U::LocalObservers, block_size::Tuple{Int, Int}; center::Tuple{Int, Int}=(1,1)) = operator_splitting(
	U, center_splitting(size(U), block_size, center))


# function center_splitting(U_ori::LocalObservers{M}, block_size::Tuple{Int, Int}, center_size::Tuple{Int, Int}=(1,1)) where M
# 	ori_partition = block_partition(size(U_ori), block_size)
# 	U = default_embeding(U_ori, size(ori_partition))
# 	center_size = min_shape(size(U_ori), center_size)

# 	@assert (block_size[1] >= center_size[1]+2 ) && (block_size[2] >= center_size[2]+2 )
# 	m, n = size(U)
# 	Htable = get_table(U.data)


# 	blks = BlockLocalOperator{M}[]
# 	for k2 in 0:n-1
# 		for k1 in 0:m-1
# 			Uk = _empty_local_observables(M, m, n)
# 			partition = VHshift(ori_partition, k1, k2)
# 			for i in 1:nrows(partition) 
# 				_rows = rowindices(partition, i)
# 				for j in 1:ncols(partition)
# 					_cols = colindices(partition, j)
# 					center = find_center(_rows, _cols, center_size)
# 					if !isnothing(center)
# 						_center_rows, _center_cols = center
# 						get_center!(Uk, U, _center_rows, _center_cols, Htable)
# 					end
# 				end
# 			end
# 			if n_nontrivial_terms(Uk) > 0
# 				push!(blks, BlockLocalOperator(Uk, partition) )
# 			end
# 		end
# 	end

# 	all(Htable) || error("incompatible splitting.")

# 	check_observer_splitting(U, blks) || error("incompatible splitting.")
# 	return blks
# end

# function get_center!(Uk::LocalObservers{M}, U::LocalObservers{M}, rows, cols, Htable) where M
# 	for i in 1:length(rows)
# 		row_i = rows[i]
# 		for j in 1:length(cols)
# 			col_j = cols[j]
# 			if !(Htable[row_i, col_j])
# 				Uk.data[row_i, col_j] = U.data[row_i, col_j]
# 				Htable[row_i, col_j] = true
# 			end
# 		end
# 	end
# end

function _empty_local_observables(::Type{M}, m::Int, n::Int) where M
	data = Matrix{Union{M, Nothing}}(nothing, m, n)
	return LocalObservers(PeriodicArray(data))
end

function check_observer_splitting(x::LocalObservers{M}, blks::Vector{BlockLocalOperator{M}}) where {M}
	# return n_nontrivial_terms(x) == sum(n_nontrivial_terms.(blks)) 
	return n_nontrivial_terms(x) == sum(n_nontrivial_terms.(blks)) == sum([sum(n_nontrivial_terms.(subblocks(blk))) for blk in blks])
end


function check_local_observers(V::LocalObservers, partition::SquareLatticePartition) 
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