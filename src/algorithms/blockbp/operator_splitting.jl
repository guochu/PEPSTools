
# function default_embeding(x::AbstractMatrix{Union{M, Nothing}}, shape::Tuple{Int, Int}) where M
# 	old_shape = size(x)
# 	new_shape = max_shape(old_shape, shape)
# 	x_new = Matrix{Union{M, Nothing}}(nothing, new_shape)
# 	x_new[get_center_positions(new_shape[1], old_shape[1]), get_center_positions(new_shape[2], old_shape[2])] = x
# 	return x_new
# end
# default_embeding(U::SquareLattice, shape::Tuple{Int, Int}) = SquareLattice(default_embeding(U.V, shape), default_embeding(U.H, shape))

starting_partition(shape::Tuple{Int, Int}, block_size::Tuple{Int, Int}, center_size::Tuple{Int, Int}) = lattice_partition(
	shape, block_size, (1, 1), center_size)

function default_splitting(shape::Tuple{Int, Int}, block_size::Tuple{Int, Int})
	odd_partition = starting_partition(shape, block_size, min_shape(shape, block_size))
	even_partition = VHshift(odd_partition, 1, 1)
	pts = [even_partition, odd_partition]
	is_valid_hamiltonian_splitting(shape, pts) || error("invalid hamiltonian splitting.")
	return pts
end

function center_splitting(shape::Tuple{Int, Int}, block_size::Tuple{Int, Int}, center_size::Tuple{Int, Int})
	center_size = min_shape(shape, center_size)
	m, n = shape
	ori_partition = starting_partition(shape, block_size, center_size)
	pts = SquareLatticePartition[]
	for k2 in 0:n-1
		for k1 in 0:m-1
			partition = VHshift(ori_partition, k1, k2)
			push!(pts, partition)
		end
	end
	# is_valid_hamiltonian_splitting(shape, pts) || error("invalid hamiltonian splitting.")
	return pts
end

"""
	operator_splitting(x::SquareLatticeBonds{Union{M, Nothing}}, pts::Vector{SquareLatticePartition})

Spliting the Hamiltonian into block operators

If block size is larger than the hamiltonian size, then the hamiltonian is embeded into 
the center of the block, otherwise the starting point of the hamiltonian is (1,1)
"""
function operator_splitting(x::SquareLatticeBonds{Union{M, Nothing}}, pts::Vector{SquareLatticePartition}) where M
	isempty(pts) && error("no partition.")
	m, n = size(pts[1])
	Htable = get_table(x.H) 
	Vtable = get_table(x.V)

	blks = BlockOperator{M}[]
	for pt in pts
		Uk = _empty_square_lattice(M, m, n)
		fill_block!(Uk, x, pt, Htable, Vtable)

		if n_nontrivial_terms(Uk) > 0
			push!(blks, BlockOperator(Uk, pt))
		end
	end

	(all(Htable) && all(Vtable)) || error("incompatible splitting.")

	# println("number of blocks $(length(blks))")

	check_operator_splitting(x, blks) || error("incompatible splitting.")
	return blks
end

function shifted_starting(new_shape::Tuple{Int, Int}, ori_shape::Tuple{Int, Int})
	s1 = div(new_shape[1] - ori_shape[1], 2)
	s2 = div(new_shape[2] - ori_shape[2], 2)
	return s1, s2	
end

function fill_block!(Uk, U, blk, Htable, Vtable)
	s1, s2 = shifted_starting(size(Uk), size(U))
	vbonds_rows = vbonds_rowcenterindices(blk)
	vbonds_cols = vbonds_colcenterindices(blk)
	hbonds_rows = hbonds_rowcenterindices(blk)
	hbonds_cols = hbonds_colcenterindices(blk)

	for i in hbonds_rows
		for j in hbonds_cols
			if !(Htable[i-s1, j-s2])
				Uk.H[i, j] = U.H[i-s1, j-s2]
				Htable[i-s1, j-s2] = true
			end
		end
	end

	for i in vbonds_rows
		for j in vbonds_cols
			if !(Vtable[i-s1, j-s2])
				Uk.V[i, j] = U.V[i-s1, j-s2]
				Vtable[i-s1, j-s2] = true
			end
		end
	end

end

default_splitting(x::SquareLatticeBonds, block_size::Tuple{Int, Int}) = operator_splitting(x, default_splitting(size(x), block_size)) 
function center_splitting(x::SquareLatticeBonds, block_size::Tuple{Int, Int}; center::Tuple{Int, Int}=(2,2))
	pts = center_splitting(size(x), block_size, center)
	is_valid_hamiltonian_splitting(size(x), pts) || error("invalid hamiltonian splitting.")
	return operator_splitting(x, pts)
end 



function check_operator_splitting(x::SquareLatticeBonds{Union{M, Nothing}}, blks::Vector{BlockOperator{M}}) where {M}
	return n_nontrivial_terms(x) == sum(n_nontrivial_terms.(blks)) == sum([sum(n_nontrivial_terms.(subblocks(blk))) for blk in blks])
end

function get_table(H::AbstractMatrix{Union{M, Nothing}}) where M
	Htable = PeriodicArray{Bool,2}(undef, size(H))
	Htable .= false
	for i in 1:length(H)
		if isnothing(H[i])
			Htable[i] = true
		end
	end
	return Htable
end

function is_valid_hamiltonian_splitting(shape::Tuple{Int, Int}, blks::Vector{SquareLatticePartition})
	isempty(blks) && error("no partition.")
	for i in 2:length(blks)
		(size(blks[i]) == size(blks[1])) || error("partition size mismatch.")
	end
	new_shape = size(blks[1])
	@assert (shape[1] <= new_shape[1]) && (shape[2] <= new_shape[2])
	# m, n = new_shape
	# s1 = div(new_shape[1] - shape[1], 2)
	# s2 = div(new_shape[2] - shape[2], 2)

	vbonds_rows = Set{Int}()
	vbonds_cols = Set{Int}()
	hbonds_rows = Set{Int}()
	hbonds_cols = Set{Int}()
	for blk in blks
		for item in vbonds_rowcenterindices(blk)
			push!(vbonds_rows, item)
		end
		for item in vbonds_colcenterindices(blk)
			push!(vbonds_cols, item)
		end
		for item in hbonds_rowcenterindices(blk)
			push!(hbonds_rows, item)
		end
		for item in hbonds_colcenterindices(blk)
			push!(hbonds_cols, item)
		end
	end
	# println(vbonds_rows)
	# println(vbonds_cols)
	# println(hbonds_rows)
	# println(hbonds_cols)
	# println(get_center_positions(new_shape[1], shape[1]) )
	# println(get_center_positions(new_shape[2], shape[2]))

	for k in get_center_positions(new_shape[1], shape[1]) 
		((k in vbonds_rows) && (k in hbonds_rows)) || return false
	end
	for k in get_center_positions(new_shape[2], shape[2])
		((k in vbonds_cols) && (k in hbonds_cols)) || return false
	end
	return true
end
