const VectorMessage{T<:Number} = Message{Vector{T}, Vector{T}}
Base.similar(m::VectorMessage) = Message(similar(m.i), similar(m.o))
const SquareLatticeBondMessages{T} = SquareLatticeBonds{VectorMessage{T}}

scalartype(::Type{SquareLatticeBondMessages{T}}) where {T} = T
scalartype(x::SquareLatticeBondMessages) = scalartype(typeof(x))
Base.similar(m::SquareLatticeBondMessages) = SquareLatticeBonds(similar.(m.V), similar.(m.H))

function Base.getindex(m::SquareLatticeBondMessages, bond::Pair{Int, Int})
	s1, s2 = size(m)
	index = CartesianIndices(size(m))
	n1, n2 = index[bond[1]], index[bond[2]]
	row1, col1 = n1[1], n1[2]
	row2, col2 = n2[1], n2[2]
	if row1 == row2
		if col1 == mod1(col2 - 1, s2)
			return m.H[row1, col1].o
		elseif col1 == mod1(col2 + 1, s2)
			return m.H[row1, col2].i
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	elseif col1 == col2
		if row1 == mod1(row2 - 1, s1)
			return m.V[row1, col1].o
		elseif row1 == mod1(row2 + 1, s1)
			return m.V[row2, col1].i
		else
			throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
		end
	else
		throw(ArgumentError("node $(n1) and node $(n2) are not neighbours"))
	end
end
function Base.setindex!(m::SquareLatticeBondMessages, v::AbstractVector, bond::Pair{Int, Int})
	copy!(m[bond], v)
end
