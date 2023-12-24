
# single row environment
# see LubaschBanuls2014b PHYSICAL REVIEW B 90, 064425 (2014)

"""
	struct PEPSRowEnv{T}
	left to right (top down also converted to left right)
"""
struct PEPSRowEnv{T, _MPS} <: AbstractPEPSRowEnv
	up::_MPS
	middle::Vector{Array{T, 5}}
	down::_MPS
	hstorage::Vector{Array{T, 3}}
end

function PEPSRowEnv(up::MPS, middle::Vector{<:AbstractArray}, down::MPS, left::AbstractArray{<:Number, 3}, right::AbstractArray{<:Number, 3})
	T = eltype(up)
	middle = convert(Vector{Array{T, 5}}, middle)
	return PEPSRowEnv(up, middle, down, compute_hstorage_right(up, middle, down, left, right))
end 


row_environments(up::MPS, middle::Vector{<:AbstractArray{<:Number, 5}}, down::MPS, left::AbstractArray{<:Number, 3}, 
	right::AbstractArray{<:Number, 3}) = PEPSRowEnv(up, middle, down, left, right)


function update_left!(x::PEPSRowEnv, pos::Int)
    x.hstorage[pos+1] = normalize!(bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], x.middle[pos]))
end

# the last one is assumed to be nothing since this function is only use dto compute the energy
function update!(x::PEPSRowEnv, Us::Vector; kwargs...)
	@assert length(Us) == length(x)
	for pos in 1:length(Us)-1
		# println("updating the $pos pair...")
		update_single!(x, pos, Us[pos]; kwargs...)
	end
end


function row_expectation(x::PEPSRowEnv, Us::Vector)
	@assert length(Us) == length(x)
	return [row_expectation_single(x, pos, Us[pos]) for pos in 1:length(Us)-1]
end



function local_expectations(x::PEPSRowEnv, Us::Vector{M}) where {M <: Union{AbstractMatrix, Nothing}}
	@assert length(x) == length(Us)
	# n = compute_norm(x)
	return [local_expectation_single(x, pos, Us[pos]) for pos in 1:length(Us)]
end

function local_expectation_single(x::PEPSRowEnv, pos::Int, U::AbstractMatrix)
	tnj = sandwich_single(x.middle[pos], U, x.middle[pos]) 
	hleft = bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], tnj)
	nleft = bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], x.middle[pos])
	@tensor r = hleft[1,2,3] * x.hstorage[pos+1][1,2,3]
	@tensor n = nleft[1,2,3] * x.hstorage[pos+1][1,2,3]
	update_left!(x, pos)
	return r / n
end
function local_expectation_single(x::PEPSRowEnv, pos::Int, U::Nothing)
	update_left!(x, pos)
	return 0.
end

# function compute_norm(x::PEPSRowEnv)
# 	hright = bm_update_right(x.hstorage[2], x.up[2], x.down[2], x.middle[1])
# 	@tensor n = x.hstorage[1][1,2,3] * hright[1,2,3]
# 	return n
# end

function compute_center(x::PEPSRowEnv, pos::Int)
	# hleft = x.hstorage[pos]
	# hright = x.hstorage[pos+2]
	AL = x.middle[pos]
	AR = x.middle[pos+1]
	left_q, aL = tqr!(copy(AL), (2,3,5), (1,4))
	aR, right_q = tlq!(copy(AR), (2,1), (3,4,5))

	@tensor ml[1,5,2,6,4,8,3,7] := conj(left_q[1,2,3,4]) * left_q[5,6,7,8]
	ml4 = tie(ml, (2,2,2,2))
	@tensor mr[1,5,2,6,3,7,4,8] := conj(right_q[1,2,3,4]) * right_q[5,6,7,8]
	mr4 = tie(mr, (2,2,2,2))

	hleft = bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], ml4)
	hright = bm_update_right(x.hstorage[pos+2], x.up[pos+2], x.down[pos+2], mr4)

	@tensor N[2,4] := hleft[1,2,3] * hright[1,4,3]
	X = make_N_positive(N)

	# aL, aR, X = _regularize(aL, aR, X)
	return left_q, aL, right_q, aR, X
end

rdm1s(x::PEPSRowEnv) = [row_rdm1_single(x, pos) for pos in 1:length(x)]
function row_rdm1_single(x::PEPSRowEnv, pos::Int)
	AL = x.middle[pos]
	s1, s2 = size(AL, 1), size(AL, 2)
	@tensor tmp[6,1,2,7,3,8,4,9,5,10] := conj(AL[1,2,3,4,5]) * AL[6,7,8,9,10]
	tmp4 = tie(tmp, (4,2,2,2))
	hright = bm_update_right(x.hstorage[pos+1], x.up[pos+1], x.down[pos+1], tmp4)
	hright4 = reshape(hright, (size(hright, 1), s1, s1, s2*s2, size(hright, 3) ) )
	@tensor rho[4, 5] := x.hstorage[pos][1,2,3] * hright4[1,4,5, 2,3]

	update_left!(x, pos)
	return normalize_rho!(rho)
end
rdm2s(x::PEPSRowEnv) = [row_rdm2_single(x, pos) for pos in 1:length(x)-1]


function normalize_rho!(rho::AbstractArray{<:Number, 4})
	s1, s2 = size(rho, 1), size(rho, 2)
	rho2 = reshape(rho, (s1*s2, s1*s2))
	rho ./= tr(rho2)
	return rho
end
function normalize_rho!(rho::AbstractArray{<:Number, 2})
	rho ./= tr(rho)
	return rho
end

# boundary mps update from left
bm_update_left(hleft::AbstractArray{<:Number, 3}, mpsu::AbstractArray{<:Number, 3}, mpsd::AbstractArray{<:Number, 3}, m::AbstractArray{<:Number, 5}) = bm_update_left(
	hleft, mpsu, mpsd, sandwich_single(m, m))
function bm_update_left(hleft::AbstractArray{<:Number, 3}, mpsu::AbstractArray{<:Number, 3}, mpsd::AbstractArray{<:Number, 3}, mm::AbstractArray{<:Number, 4})
	@tensor a[8,7,5] := ((hleft[1,2,3] * mpsd[3,4,5]) * mm[2,6,7,4]) * mpsu[1,6,8]
	return a
end

bm_update_right(hright::AbstractArray{<:Number, 3}, mpsu::AbstractArray{<:Number, 3}, mpsd::AbstractArray{<:Number, 3}, m::AbstractArray{<:Number, 5}) = bm_update_right(
	hright, mpsu, mpsd, sandwich_single(m, m))
function bm_update_right(hright::AbstractArray{<:Number, 3}, mpsu::AbstractArray{<:Number, 3}, mpsd::AbstractArray{<:Number, 3}, mm::AbstractArray{<:Number, 4})
	@tensor a[1,6,8] := ((mpsu[1,2,3] * hright[3,4,5]) * mm[6,2,4,7]) * mpsd[8,7,5]
	return a
end


function compute_hstorage_right(mpsu::MPS, middle::Vector, mpsd::MPS, left::AbstractArray{<:Number, 3}, right::AbstractArray{<:Number, 3})
	@assert length(mpsu) == length(mpsd) == length(middle) + 2
	L = length(mpsu) - 2
	hstorage = Vector{Array{eltype(mpsu), 3}}(undef, L+1)
	mpsu_left = dropdims(mpsu[1], dims=1)
	mpsd_left = dropdims(mpsd[1], dims=1)
	@tensor tmpl[2,3,5] := mpsu_left[1,2] * left[1,3,4] * mpsd_left[4,5]
	hstorage[1] = normalize!(tmpl)
	mpsu_right = dropdims(mpsu[L+2], dims=3)
	mpsd_right = dropdims(mpsd[L+2], dims=3)
	@tensor tmpr[1,3,5] := mpsu_right[1,2] * right[2,3,4] * mpsd_right[5,4]
	hstorage[L+1] = normalize!(tmpr)
	for i in L:-1:2
		hstorage[i] = normalize!(bm_update_right(hstorage[i+1], mpsu[i+1], mpsd[i+1], middle[i]))
	end
	return hstorage
end
