
"""
	struct SquareLatticeBonds{M}
"""
struct SquareLatticeBonds{M}
	V::PeriodicArray{M,2}
	H::PeriodicArray{M,2}

function SquareLatticeBonds{M}(V::PeriodicArray, H::PeriodicArray) where {M}
	(size(V) == size(H)) || throw(ArgumentError("V and H size mismatch"))
	new{M}(V, H)
end

end

SquareLatticeBonds(V::PeriodicArray{M, 2}, H::PeriodicArray{M, 2}) where {M} = SquareLatticeBonds{M}(V, H)
SquareLatticeBonds(V::AbstractMatrix{M}, H::AbstractMatrix{M}) where M = SquareLatticeBonds(PeriodicArray(V), PeriodicArray(H))
SquareLatticeBonds(; V::AbstractMatrix{M}, H::AbstractMatrix{M}) where M = SquareLatticeBonds(V, H)

Base.size(x::SquareLatticeBonds) = size(x. H)
Base.size(x::SquareLatticeBonds, i::Int) = size(x.H, i)
Base.repeat(x::SquareLatticeBonds, i::Int...) = SquareLatticeBonds(repeat(x.V.data, i...), repeat(x.H.data, i...))
Base.copy(x::SquareLatticeBonds) = SquareLatticeBonds(copy(x.V), copy(x.H))

function Base.getindex(m::SquareLatticeBonds, bond::Tuple{Int, Int})
	s1, s2 = size(m)
	index = CartesianIndices(size(m))
	n1, n2 = index[bond[1]], index[bond[2]]
	row1, col1 = n1[1], n1[2]
	row2, col2 = n2[1], n2[2]
	if row1 == row2
		if col1 == mod1(col2 - 1, s2)
			return m.H[row1, col1]
		elseif col1 == mod1(col2 + 1, s2)
			return m.H[row1, col2]
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	elseif col1 == col2
		if row1 == mod1(row2 - 1, s1)
			return m.V[row1, col1]
		elseif row1 == mod1(row2 + 1, s1)
			return m.V[row2, col1]
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	else
		throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
	end
end

function Base.setindex!(m::SquareLatticeBonds, v, bond::Tuple{Int, Int})
	s1, s2 = size(m)
	index = CartesianIndices(size(m))
	n1, n2 = index[bond[1]], index[bond[2]]
	row1, col1 = n1[1], n1[2]
	row2, col2 = n2[1], n2[2]
	if row1 == row2
		if col1 == mod1(col2 - 1, s2)
			m.H[row1, col1] = v
		elseif col1 == mod1(col2 + 1, s2)
			m.H[row1, col2] = v
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	elseif col1 == col2
		if row1 == mod1(row2 - 1, s1)
			m.V[row1, col1] = v
		elseif row1 == mod1(row2 + 1, s1)
			m.V[row2, col1] = v
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	else
		throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
	end
end

const SquareLatticeHamiltonian{T<:Number} = SquareLatticeBonds{Vector{Tuple{Matrix{T}, Matrix{T}}}}
const SquareLatticeOperator{T<:Number} = SquareLatticeBonds{Union{Array{T, 4}, Nothing}}

scalartype(::Type{SquareLatticeHamiltonian{T}}) where T = T
scalartype(x::SquareLatticeHamiltonian) = scalartype(typeof(x))
scalartype(::Type{SquareLatticeOperator{T}}) where T = T
scalartype(x::SquareLatticeOperator) = scalartype(typeof(x))


is_h_periodic(x::SquareLatticeOperator) = (size(x, 2)==1) || (!any(isnothing.(x.H[:, end])))
is_h_nonperiodic(x::SquareLatticeOperator) = (size(x, 2)==1) || all(isnothing.(x.H[:, end]))
is_v_periodic(x::SquareLatticeOperator) = (size(x, 1)==1) || (!any(isnothing.(x.V[end, :])))
is_v_nonperiodic(x::SquareLatticeOperator) = (size(x, 1)==1) || all(isnothing.(x.V[end, :]))

is_periodic(x::SquareLatticeOperator) = is_h_periodic(x) && is_v_periodic(x)
is_nonperiodic(x::SquareLatticeOperator) = is_h_nonperiodic(x) && is_v_nonperiodic(x)
n_nontrivial_terms(x::SquareLatticeBonds{Union{M, Nothing}}) where M = n_nontrivial_terms(x.H) + n_nontrivial_terms(x.V)

function SquareLatticeHamiltonian(::Type{T}, m::Int, n::Int) where {T <: Number}
	H = PeriodicArray{Vector{Tuple{Matrix{T}, Matrix{T}}}, 2}(undef, m, n)
	V = PeriodicArray{Vector{Tuple{Matrix{T}, Matrix{T}}}, 2}(undef, m, n)
	for i in 1:length(H)
		H[i] = Vector{Tuple{Matrix{T}, Matrix{T}}}()
	end
	for i in 1:length(V)
		V[i] = Vector{Tuple{Matrix{T}, Matrix{T}}}()
	end
	return SquareLatticeHamiltonian{T}(V, H)
end

function SquareLatticeOperator(::Type{T}, m::Int, n::Int) where {T <: Number}
	H = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, m, n)
	V = PeriodicArray{Union{Array{T, 4}, Nothing}, 2}(nothing, m, n)
	return SquareLatticeOperator{T}(V, H)
end

function exponential(h::SquareLatticeHamiltonian{T}, dt::Number) where {T <: Number}
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
	return SquareLatticeBonds(V, H)
end
function exponential(h::SquareLatticeOperator{T}, dt::Number) where {T <: Number}
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
	return SquareLatticeBonds(V, H)

end

function squeeze(h::SquareLatticeHamiltonian{T}) where {T <: Number}
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
	return SquareLatticeBonds(V, H)
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

n_nontrivial_terms(x::AbstractArray{Union{M, Nothing}}) where M = sum((!).(isnothing.(x)))

