
"""
	MPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct MPO{T<:Number} <: Abstract1DTN{T}
	data::Vector{Array{T, 4}}
	scaling::Ref{Float64}

"""
	MPO{A}(mpotensors::Vector)
Constructor entrance for MPO, which only supports strictly quantum number conserving operators

site tensor convention:
i mean in arrow, o means out arrow
    o 
    |
    2
o-1   3-i
	4
	|
	i
The left and right boundaries are always vacuum.
The case that the right boundary is not vacuum corresponds to operators which do not conserve quantum number, 
such as a†, this case is implemented with another MPO object.
"""
function MPO{T}(mpotensors::AbstractVector, scaling::Ref{Float64}) where {T<:Number}
	_check_mpo_space(mpotensors)
	return new{T}(convert(Vector{Array{T, 4}}, mpotensors), scaling)
end

end
MPO(data::AbstractVector{<:MPOTensor{T}}; scaling::Real=1) where {T <: Number} = MPO{T}(data, Ref(convert(Float64, scaling)))


# attributes
ophysical_dimensions(psi::MPO) = [size(item, 2) for item in psi.data]
iphysical_dimensions(psi::MPO) = [size(item, 4) for item in psi.data]


function _check_mpo_space(mpotensors::Vector)
	L = length(mpotensors)
	for i in 1:L-1
		(space_r(mpotensors[i]) == space_l(mpotensors[i+1])) || throw(DimensionMismatch())
	end
	# boundaries should be dimension 
	(space_l(mpotensors[1]) == 1) || throw(DimensionMismatch())
	(space_r(mpotensors[L]) == 1) || throw(DimensionMismatch())
	return true
end


# initializers

"""
	randommpo(::Type{T}, dy::Vector{Int}, dx::Vector{Int}; D::Int) where {T<:Number}
	dy are the input dimensions, dx are the output dimensions
"""
function randommpo(::Type{T}, dx::Vector{Int}, dy::Vector{Int}; D::Int) where {T<:Number}
	(length(dx) == length(dy)) || throw(DimensionMismatch())
	L = length(dx)
	r = Vector{Array{T, 4}}(undef, L)
	r[1] = randn(T, 1, dx[1], D, dy[1])
	r[L] = randn(T, D, dx[L], 1, dy[L])
	for i in 2:L-1
		r[i] = randn(T, D, dx[i], D, dy[i])
	end
	return MPO(r)
end 
randommpo(::Type{T}, physpaces::Vector{Int}; kwargs...) where {T<:Number} = randommpo(T, physpaces, physpaces; kwargs...)
randommpo(::Type{T}, L::Int; d::Int, D::Int) where {T<:Number} = randommpo(T, [d for i in 1:L], D=D)
randommpo(L::Int; kwargs...) = randommpo(Float64, L; kwargs...)
