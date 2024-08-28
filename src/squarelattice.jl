
"""
	struct SquareLattice{M}
"""
struct SquareLattice{M}
	V::PeriodicArray{M,2}
	H::PeriodicArray{M,2}
end

SquareLattice(V::AbstractMatrix{M}, H::AbstractMatrix{M}) where M = SquareLattice(PeriodicArray(V), PeriodicArray(H))
SquareLattice(; V::AbstractMatrix{M}, H::AbstractMatrix{M}) where M = SquareLattice(V, H)

Base.size(x::SquareLattice) = size(x. H)
Base.size(x::SquareLattice, i::Int) = size(x.H, i)
Base.repeat(x::SquareLattice, i::Int...) = SquareLattice(repeat(x.V.data, i...), repeat(x.H.data, i...))


const SquareLatticeHamiltonianBase{T<:Number} = SquareLattice{Vector{Tuple{Matrix{T}, Matrix{T}}}}
const SquareLatticeOperatorBase{T<:Number} = SquareLattice{Union{Array{T, 4}, Nothing}}

scalartype(::Type{SquareLatticeHamiltonianBase{T}}) where T = T
scalartype(x::SquareLatticeHamiltonianBase) = scalartype(typeof(x))
scalartype(::Type{SquareLatticeOperatorBase{T}}) where T = T
scalartype(x::SquareLatticeOperatorBase) = scalartype(typeof(x))


is_h_periodic(x::SquareLatticeOperatorBase) = (size(x, 2)==1) || (!any(isnothing.(x.H[:, end])))
is_h_nonperiodic(x::SquareLatticeOperatorBase) = (size(x, 2)==1) || all(isnothing.(x.H[:, end]))
is_v_periodic(x::SquareLatticeOperatorBase) = (size(x, 1)==1) || (!any(isnothing.(x.V[end, :])))
is_v_nonperiodic(x::SquareLatticeOperatorBase) = (size(x, 1)==1) || all(isnothing.(x.V[end, :]))

is_periodic(x::SquareLatticeOperatorBase) = is_h_periodic(x) && is_v_periodic(x)
is_nonperiodic(x::SquareLatticeOperatorBase) = is_h_nonperiodic(x) && is_v_nonperiodic(x)
nontrivial_terms(x::SquareLattice{Union{M, Nothing}}) where M = nontrivial_terms(x.H) + nontrivial_terms(x.V)

function SquareLatticeHamiltonian(::Type{T}, m::Int, n::Int) where {T <: Number}
	H = PeriodicArray{Vector{Tuple{Matrix{T}, Matrix{T}}}, 2}(undef, m, n)
	V = PeriodicArray{Vector{Tuple{Matrix{T}, Matrix{T}}}, 2}(undef, m, n)
	for i in 1:length(H)
		H[i] = Vector{Tuple{Matrix{T}, Matrix{T}}}()
	end
	for i in 1:length(V)
		V[i] = Vector{Tuple{Matrix{T}, Matrix{T}}}()
	end
	return SquareLatticeHamiltonianBase{T}(V, H)
end

function SquareLatticeOperator(::Type{T}, m::Int, n::Int) where {T <: Number}
	H = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, m, n)
	V = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, m, n)
	return SquareLatticeOperatorBase{T}(V, H)
end

function exponential(h::SquareLatticeHamiltonianBase{T}, dt::Number) where {T <: Number}
	m, n = size(h)
	_T = promote_type(T, typeof(dt))
	# r = SquareLatticeOperator(promote_type(T, typeof(dt)), m, n)
	H = PeriodicArray{Union{Array{_T, 4}, Nothing}, 2}(nothing, size(h.H))
	V = PeriodicArray{Union{Array{_T, 4}, Nothing}, 2}(nothing, size(h.V))
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			H[i, j] = _get_mat_exp(h.H[i, j], dt)
		end
	end
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			V[i, j] = _get_mat_exp(h.V[i, j], dt)
		end
	end
	return SquareLattice(V, H)
end
function exponential(h::SquareLatticeOperatorBase{T}, dt::Number) where {T <: Number}
	m, n = size(h)
	_T = promote_type(T, typeof(dt))
	H = PeriodicArray{Union{Array{_T, 4}, Nothing}, 2}(nothing, size(h.H))
	V = PeriodicArray{Union{Array{_T, 4}, Nothing}, 2}(nothing, size(h.V))
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			item = h.H[i, j]
			H[i, j] = isnothing(item) ? nothing : reshape(exp(tie(item, (2,2)) * dt), size(item))
		end
	end
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			item = h.V[i, j]
			V[i, j] = isnothing(item) ? nothing : reshape(exp(tie(item, (2,2)) * dt), size(item))
		end
	end
	return SquareLattice(V, H)

end

function squeeze(h::SquareLatticeHamiltonianBase{T}) where {T <: Number}
	m, n = size(h)
	H = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, size(h.H))
	V = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, size(h.V))
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			H[i, j] = _get_mat(h.H[i, j])
		end
	end
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			V[i, j] = _get_mat(h.V[i, j])
		end
	end
	return SquareLattice(V, H)
end


function _get_mat(m::Vector{Tuple{Matrix{T}, Matrix{T}}}) where T
	isempty(m) && return nothing
	@tensor a[1,3,2,4] := m[1][1][1,2] * m[1][2][3,4]
	for i in 2:length(m)
		@tensor a[1,3,2,4] += m[i][1][1,2] * m[i][2][3,4]
	end
	return a
end

function _get_mat_exp(m::Vector{Tuple{Matrix{T}, Matrix{T}}}, dt::Number) where T
	a = _get_mat(m)
	isnothing(a) && return a
	s = size(a)
	a2 = exp(tie(a, (2, 2)) * dt)
	return reshape(a2, s)
end

nontrivial_terms(x::AbstractArray{Union{M, Nothing}}) where M = sum((!).(isnothing.(x)))

