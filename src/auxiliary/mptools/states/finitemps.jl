struct MPS{T<:Number, R<:Real} <: Abstract1DTN{T}
	data::Vector{Array{T, 3}}
	s::Vector{Union{Missing, Vector{R}}}
	scaling::Ref{Float64}

function MPS{T, R}(data::AbstractVector, svectors::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	(length(data)+1 == length(svectors)) || throw(DimensionMismatch("length of singular vectors must be length of site tensors+1"))
	_check_mps_space(data)
	new{T, R}(convert(Vector{Array{T, 3}}, data), convert(Vector{Union{Missing, Vector{R}}}, svectors), scaling)
end
end

function MPS{T, R}(data::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	_check_mps_space(data)
	svectors = Vector{Union{Missing, Vector{R}}}(undef, length(data)+1)
	svectors[1] = ones(space_l(data[1]))
	svectors[end] = ones(space_r(data[end]))
	return MPS{T, R}(convert(Vector{Array{T, 3}}, data), svectors, scaling)
end

function MPS(data::AbstractVector{<:MPSTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T <: Number}
	R = real(T)
	return MPS{T, R}(data, svectors, Ref(convert(R, scaling)))
end 
function MPS(data::AbstractVector{<:MPSTensor{T}}; scaling::Real=1) where {T <: Number}
	R = real(T)
	return MPS{T, R}(data, Ref(convert(R, scaling)))
end

Base.copy(psi::MPS) = MPS(copy(psi.data), copy(psi.s), scaling=scaling(psi))

svectors_uninitialized(psi::MPS) = any(ismissing, psi.s)
function unset_svectors!(psi::MPS)
	psi.s[2:end-1] .= missing
	return psi
end

function _check_mps_space(mpstensors::Vector)
	L = length(mpstensors)
	for i in 1:L-1
		(space_r(mpstensors[i]) == space_l(mpstensors[i+1])) || throw(DimensionMismatch())
	end
	(space_l(mpstensors[1]) == 1) || throw(DimensionMismatch("left boundary should be size 1."))
	(space_r(mpstensors[L]) == 1) || throw(DimensionMismatch("right boundary should be size 1."))
	return true
end

# attributes
physical_dimensions(psi::MPS) = [size(item, 2) for item in psi.data]


# initializers
function MPS(f, ::Type{T}, physpaces::Vector{Int}, virtualpaces::Vector{Int}) where {T <: Number}
	L = length(physpaces)
	(length(virtualpaces) == L+1) || throw(DimensionMismatch())
	any(virtualpaces .== 0) &&  @warn "auxiliary space is empty"
	mpstensors = [f(T, virtualpaces[i], physpaces[i], virtualpaces[i+1]) for i in 1:L]
	return MPS(mpstensors)
end


function randommps(::Type{T}, physpaces::Vector{Int}; D::Int) where {T <: Number}
	virtualpaces = max_bond_dimensions(physpaces, D)
	return MPS(randn, T, physpaces, virtualpaces)
end
randommps(physpaces::Vector{Int}; D::Int) = randommps(Float64, physpaces; D=D)
randommps(::Type{T}, L::Int; d::Int, D::Int) where {T <: Number} = randommps(T, [d for i in 1:L]; D=D)
randommps(L::Int; d::Int, D::Int) = randommps(Float64, L; d=d, D=D)


function increase_bond!(psi::MPS; D::Int)
	if bond_dimension(psi) < D
		virtualpaces = max_bond_dimensions(physical_dimensions(psi), D)
		for i in 1:length(psi)
			sl = max(min(virtualpaces[i], D), size(psi[i], 1))
			sr = max(min(virtualpaces[i+1], D), size(psi[i], 3))
			m = zeros(scalartype(psi), sl, size(psi[i], 2), sr)
			m[1:size(psi[i], 1), :, 1:size(psi[i], 3)] .= psi[i]
			psi[i] = m
		end
	end
	return psi
end


function max_bond_dimensions(physpaces::Vector{Int}, D::Int) 
	L = length(physpaces)
	left = 1
	right = 1
	virtualpaces = Vector{Int}(undef, L+1)
	virtualpaces[1] = left
	for i in 2:L
		virtualpaces[i] = min(virtualpaces[i-1] * physpaces[i-1], D)
	end
	virtualpaces[L+1] = right
	for i in L:-1:2
		virtualpaces[i] = min(virtualpaces[i], physpaces[i] * virtualpaces[i+1])
	end
	return virtualpaces
end

# check is canonical
isleftcanonical(a::MPS; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
isrightcanonical(a::MPS; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::MPS; kwargs...) = is_right_canonical(psi; kwargs...)
check if the state is right-canonical, the singular vectors are also checked that whether there are the correct Schmidt numbers or not
This form is useful for time evolution for stability issue and also efficient for computing observers of unitary systems
"""
function iscanonical(psi::MPS; kwargs...)
	isrightcanonical(psi) || return false
	# we also check whether the singular vectors are the correct Schmidt numbers
	svectors_uninitialized(psi) && return false
	hold = l_LL(psi)
	for i in 1:length(psi)-1
		hold = updateleft(hold, psi[i], psi[i])
		tmp = psi.s[i+1]
		isapprox(hold, Diagonal(tmp.^2); kwargs...) || return false
	end
	return true
end
