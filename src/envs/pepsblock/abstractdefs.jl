# nonperiodic PEPS block
abstract type AbstractBlock{T} end
abstract type AbstractPEPSBlock{T} <: AbstractBlock{T} end
abstract type AbstractSquareTNBlock{T} <: AbstractBlock{T} end

Base.size(x::AbstractBlock) = size(x.peps)
Base.size(x::AbstractBlock, i::Int) = size(x.peps, i)
scalartype(::Type{<:AbstractBlock{T}}) where {T<:Number} = T
scalartype(x::AbstractBlock) = scalartype(typeof(x))


left_boundary(x::AbstractBlock) = _boundary_mps(x.left)
right_boundary(x::AbstractBlock) = _boundary_mps(x.right)
up_boundary(x::AbstractBlock) = _boundary_mps(x.up)
down_boundary(x::AbstractBlock) = _boundary_mps(x.down)


# these functions are also used for cyclic peps block so they are not typed
function sl_mpoleft_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,4,8,5,9,2,6] := conj(t[1,2,3,4,5]) * t[1,6,7,8,9]
		return tie(tmp, (2,2,2,2))
	end
	return get_tn.(x.peps[:, i])
end 
function sl_mporight_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,2,6,5,9,4,8] := conj(t[1,2,3,4,5]) * t[1,6,7,8,9]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[:, i])
end 
function sl_mpoup_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[2,6,5,9,4,8,3,7] := conj(t[1,2,3,4,5]) * t[1,6,7,8,9]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[i, :])
end
function sl_mpodown_util(x, i::Int)
	sandwich_single.(x.peps[i, :])
end

dl_mpoleft_util(x, i::Int) = [permute(item, (2,3,4,1)) for item in x.peps[:, i]]
dl_mporight_util(x, i::Int) = [permute(item, (2,1,4,3)) for item in x.peps[:, i]]
dl_mpoup_util(x, i::Int) = [permute(item, (1,4,3,2)) for item in x.peps[i, :]]
dl_mpodown_util(x, i::Int) = x.peps[i, :]

mpoleft_util(x::AbstractPEPSBlock, i::Int) = sl_mpoleft_util(x, i)
mporight_util(x::AbstractPEPSBlock, i::Int) = sl_mporight_util(x, i)
mpoup_util(x::AbstractPEPSBlock, i::Int) = sl_mpoup_util(x, i)
mpodown_util(x::AbstractPEPSBlock, i::Int) = sl_mpodown_util(x, i)

mpoleft_util(x::AbstractSquareTNBlock, i::Int) = dl_mpoleft_util(x, i)
mporight_util(x::AbstractSquareTNBlock, i::Int) = dl_mporight_util(x, i)
mpoup_util(x::AbstractSquareTNBlock, i::Int) = dl_mpoup_util(x, i)
mpodown_util(x::AbstractSquareTNBlock, i::Int) = dl_mpodown_util(x, i)

# apply mpo on the left
function mpoleft(x::AbstractBlock, i::Int)
	L = size(x.peps, 1)
	r = Vector{Array{scalartype(x), 4}}(undef, L+2)
	r[1] = permute(_insert_dim(x.up[i]), (1,4,3,2))
	r[L+2] = permute(_insert_dim(x.down[i]), (3,4,1,2))
	r[2:L+1] = mpoleft_util(x, i)
	return MPO(r)
end 
function mporight(x::AbstractBlock, i::Int)
	L = size(x.peps, 1)
	r = Vector{Array{scalartype(x), 4}}(undef, L+2)
	r[1] = _insert_dim(x.up[i])
	r[L+2] = permute(_insert_dim(x.down[i]), (3,2,1,4))
	r[2:L+1] = mporight_util(x, i)
	return MPO(r)
end 
function mpoup(x::AbstractBlock, i::Int)
	L = size(x.peps, 2)
	r = Vector{Array{scalartype(x), 4}}(undef, L+2)
	r[1] = permute(_insert_dim(x.left[i]), (1,4,3,2))
	r[L+2] = permute(_insert_dim(x.right[i]), (3,4,1,2))
	r[2:L+1] = mpoup_util(x, i)
	return MPO(r)
end
function mpodown(x::AbstractBlock, i::Int)
	L = size(x.peps, 2)
	r = Vector{Array{scalartype(x), 4}}(undef, L+2)
	r[1] = _insert_dim(x.left[i])
	r[L+2] = permute(_insert_dim(x.right[i]), (3,2,1,4))
	r[2:L+1] = mpodown_util(x, i)
	return MPO(r)
end

row_peps(x::AbstractBlock, i::Int) = x.peps[i, :]

sl_col_peps_as_row(x, i::Int) = [permute(item, (1,3,4,5,2)) for item in x.peps[:, i]]
dl_col_peps_as_row(x, i::Int) = [permute(item, (3,4,5,2)) for item in x.peps[:, i]]

col_peps_as_row(x::AbstractPEPSBlock, i::Int) = sl_col_peps_as_row(x, i)
col_peps_as_row(x::AbstractSquareTNBlock, i::Int) = dl_col_peps_as_row(x, i) 

function row(x::AbstractBlock, pos::Int, mult_alg::MPSCompression=SVDCompression())
	L = size(x, 1)

	up = up_boundary(x)
	down = down_boundary(x)

	for i in L:-1:pos+1
		mpo = mpodown(x, i) 
		down, err = mpompsmult(mpo, down, mult_alg)
		normalize!(down)
	end
	for i in 1:pos-1
		mpo = mpoup(x, i) 
		up, err = mpompsmult(mpo, up, mult_alg)
		normalize!(up)
	end
	return row_environments(up, row_peps(x, pos), down, x.left[pos], x.right[pos])
end


function col(x::AbstractBlock, pos::Int, mult_alg::MPSCompression=SVDCompression())
	L = size(x, 2)

	left = left_boundary(x)
	right = right_boundary(x)

	for i in 1:pos-1
		mpo = mpoleft(x, i)
		left, err = mpompsmult(mpo, left, mult_alg)
		normalize!(left)
	end
	for i in L:-1:pos+1
		mpo = mporight(x, i)
		right, err = mpompsmult(mpo, right, mult_alg)
		normalize!(right)
	end

	return row_environments(right, col_peps_as_row(x, i), left, permute(x.up[pos], (3,2,1)), permute(x.down[pos], (3,2,1)))
end

# insert a trivial dimension before the i-th dimension
function _insert_dim(t::AbstractArray, i::Int)
	@assert i <= ndims(t) + 1
	old_dims = size(t)
	new_dims = (old_dims[1:i-1]..., 1, old_dims[i:end]...)
	return reshape(t, new_dims)
end
_insert_dim(t::AbstractArray) = _insert_dim(t, 1)
_boundary_mps(x::MPS) = _boundary_mps(x.data)
_boundary_mps(x::Vector) = MPS([ones(1,1,1), x..., ones(1,1,1)])
