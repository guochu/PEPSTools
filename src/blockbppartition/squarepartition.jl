abstract type AbstractPartition end


"""
	struct SquareLatticePartition

Partition does not known system size, a partition can be larger than system size
"""
struct SquareLatticePartition <: AbstractPartition
	shape::Tuple{Int, Int}
	rows::Vector{Int}
	cols::Vector{Int}
	center::Tuple{Int, Int}

function SquareLatticePartition(shape::Tuple{Int, Int}, rows::Vector{Int}, cols::Vector{Int}, center::Tuple{Int, Int})
	@assert issorted(rows) && issorted(cols)
	m, n = shape
	@assert (rows[end] - rows[1] == m) && (cols[end]- cols[1] == n)
	@assert (shape[1] >= center[1]) && (shape[2] >= center[2])
	return new(shape, rows, cols, center)
end

end

# SquareLatticePartition(rows::Vector{Int}, cols::Vector{Int}) = SquareLatticePartition((rows[end]-rows[1], cols[end]-cols[1]), rows, cols)

nrows(m::SquareLatticePartition) = length(m.rows) - 1
ncols(m::SquareLatticePartition) = length(m.cols) - 1

rowindices(x::SquareLatticePartition, i::Int) = x.rows[i]+1:x.rows[i+1]
colindices(x::SquareLatticePartition, j::Int) = x.cols[j]+1:x.cols[j+1]

rowcenterindices(x::SquareLatticePartition, i::Int) = rowcenterenumerate(x, i)[2]
colcenterindices(x::SquareLatticePartition, j::Int) = colcenterenumerate(x, j)[2]

rowcenterenumerate(x::SquareLatticePartition, i::Int) = _get_center(rowindices(x, i), x.center[1])
colcenterenumerate(x::SquareLatticePartition, j::Int) = _get_center(colindices(x, j), x.center[2])


function get_center_positions(L::Int, a::Int)
	if L < a
		return Int[]
	else
		s1 = div(L - a, 2)
		return s1+1:s1+a
	end
end
function _get_center(rows, a)
	pos = get_center_positions(length(rows), a)
	return pos, rows[pos]
end

rowcenterindices(x::SquareLatticePartition) = vbonds_rowcenterindices(x)
colcenterindices(x::SquareLatticePartition) = hbonds_colcenterindices(x)
function vbonds_rowcenterindices(x::SquareLatticePartition) 
	r = Int[]
	for i in 1:nrows(x)
		append!(r, rowcenterindices(x, i)[1:end-1])
	end
	return mod1.(r, size(x, 1))
end
function vbonds_colcenterindices(x::SquareLatticePartition)
	r = Int[]
	for i in 1:ncols(x)
		append!(r, colcenterindices(x, i))
	end
	return mod1.(r, size(x, 2))
end
function hbonds_rowcenterindices(x::SquareLatticePartition)
	r = Int[]
	for i in 1:nrows(x)
		append!(r, rowcenterindices(x, i))
	end
	return mod1.(r, size(x, 1))
end
function hbonds_colcenterindices(x::SquareLatticePartition)
	r = Int[]
	for i in 1:ncols(x)
		append!(r, colcenterindices(x, i)[1:end-1])
	end
	return mod1.(r, size(x, 2))
end


function subblock(tn::AbstractMatrix, blk::SquareLatticePartition, i::Int, j::Int) 
	# @assert size(tn) == blk.shape
	# @assert (size(tn,1) <= size(blk, 1)) && (size(tn,2) <= size(blk, 2))
	@assert (i <= nrows(blk)) && (j <= ncols(blk))
	return tn[rowindices(blk, i), colindices(blk, j)]
end

Base.:(==)(x::SquareLatticePartition, y::SquareLatticePartition) = (x.shape == y.shape) && (x.rows == y.rows) && (x.cols == y.cols) && (x.center == y.center)
Base.size(x::SquareLatticePartition) = x.shape
Base.size(x::SquareLatticePartition, i::Int) = x.shape[i]

# partition_nrows(x::SquareLatticePartition) = x.rows[end] - x.rows[1]
# partition_ncols(x::SquareLatticePartition) = x.cols[end] - x.cols[1]
# partition_size(x::SquareLatticePartition) = (partition_nrows(x), partition_ncols(x))

Vshift(x::SquareLatticePartition, i::Int) = SquareLatticePartition(size(x), x.rows, x.cols .+ i, x.center)
Hshift(x::SquareLatticePartition, i::Int) = SquareLatticePartition(size(x), x.rows .+ i, x.cols, x.center)
VHshift(x::SquareLatticePartition, i::Int, j::Int) = SquareLatticePartition(size(x), x.rows .+ i, x.cols .+ j, x.center)


function block_partition(shape::Tuple{Int, Int}, block_size::Tuple{Int, Int}, start::Tuple{Int, Int}, center::Tuple{Int, Int})
	shape = max_shape(shape, block_size)
	m, n = shape
	a, b = block_size
	rows = compute_block_labels(m, a) 
	cols = compute_block_labels(n, b) 
	return SquareLatticePartition(shape, rows .+ (start[1]-1), cols .+ (start[2]-1), center)	
end
block_partition(shape::Tuple{Int, Int}, block_size::Tuple{Int, Int}; start::Tuple{Int, Int}=(1,1), center::Tuple{Int, Int}=(2,2)) = block_partition(
	shape, block_size, start, center)

function compute_block_labels(m::Int, a::Int)
	@assert (2 <= a <= m)
	# if a >= m
	# 	return [0, a]
	# end
	nrow, row_rest = div(m, a), m % a
	(row_rest != 1) || error("incompatible partition.")
	v = Int[]
	for i in 1:nrow+1
		push!(v, (i-1) * a)
	end
	if row_rest > 0
		push!(v, m)
	end
	return v
end

max_shape(x::Tuple{Int, Int}, y::Tuple{Int, Int}) = (max(x[1], y[1]), max(x[2], y[2]))
min_shape(x::Tuple{Int, Int}, y::Tuple{Int, Int}) = (min(x[1], y[1]), min(x[2], y[2]))

