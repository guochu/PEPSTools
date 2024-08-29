"""
	struct BorderedPEPS{T, _MPS}

PEPS bordered by four MPSs
"""
struct BorderedPEPS{T, _MPS} <: AbstractBorderedPEPS{T}
	peps::Matrix{Array{T, 5}}
	left::_MPS
	right::_MPS
	up::_MPS
	down::_MPS

function BorderedPEPS(peps::AbstractMatrix{Array{T, 5}}, left::MPS, right::MPS, up::MPS, down::MPS) where {T<:Number}
	m, n = size(peps)
	@assert length(left) == length(right) == m
	@assert length(up) == length(down) == n
	for i in 1:m
		@assert size(left[i], 2) == size(peps[i, 1], 1)^2
		@assert size(right[i], 2) == size(peps[i, n], 3)^2
	end
	for i in 1:n
		@assert size(up[i], 2) == size(peps[1, i], 2)^2
		@assert size(down[i], 2) == size(peps[m, i], 4)^2
	end
	new{T, MPS{T, real(T)}}(peps, left, right, up, down)
end

end


function borderedpeps(peps::AbstractMatrix{Array{T, 5}}; left::MPS=trivial_mps(T, size(peps, 1)), right::MPS=trivial_mps(T, size(peps, 1)), 
	up::MPS=trivial_mps(T, size(peps, 2)), down::MPS=trivial_mps(T, size(peps, 2))) where {T <: Number}
	# @assert check(peps)
	return BorderedPEPS(peps, left, right, up, down)
end
borderedpeps(peps::PEPS, args...; kwargs...) = borderedpeps(raw_data(peps), args...; kwargs...)


function random_boundary_mps(::Type{T}, ds::AbstractVector{Int}; D::Int) where {T <: Number}
	L = length(ds)
	r = Vector{Array{T, 3}}(undef, L)
	for i in 1:L
		Dl = (i==1) ? 1 : D
		Dr = (i==L) ? 1 : D
		tmp = randn(T, Dl, ds[i], Dr, D) / D
		@tensor tmp2[1,5,2,6,3,7] := conj(tmp[1,2,3,4]) * tmp[5,6,7,4]
		r[i] = tie(tmp2, (2,2,2))
	end
	return MPS(r)
end
random_boundary_mps(::Type{T}, L::Int; D::Int) where {T <: Number} = random_boundary_mps(T, [D for i in 1:L]; D=D)

trivial_mps(::Type{T}, L::Int) where {T <: Number} = MPS([ones(T, 1, 1, 1) for i in 1:L])


mpoleft_util(x::BorderedPEPS, i::Int) = dl_mpoleft_util(x, i)
mporight_util(x::BorderedPEPS, i::Int) = dl_mporight_util(x, i)
mpoup_util(x::BorderedPEPS, i::Int) = dl_mpoup_util(x, i)
mpodown_util(x::BorderedPEPS, i::Int) = dl_mpodown_util(x, i)


function dl_mpoleft_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,4,8,5,9,2,6] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end
	return get_tn.(x.peps[:, i])
end 
function dl_mporight_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,2,6,5,9,4,8] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[:, i])
end 
function dl_mpoup_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[2,6,5,9,4,8,3,7] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[i, :])
end
function dl_mpodown_util(x, i::Int)
	sandwich_single.(x.peps[i, :])
end

col_peps_as_row(x::BorderedPEPS, i::Int) = dl_col_peps_as_row(x, i)
dl_col_peps_as_row(x, i::Int) = [permute(item, (2,3,4,1,5)) for item in x.peps[:, i]]

