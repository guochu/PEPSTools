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




neighbors(m::SquareLatticeBondMessages, node::Int) = square_neighbors(size(m), node)
edges(m::SquareLatticeBondMessages) = square_edges(size(m))

function square_neighbors(shape::Tuple{Int, Int}, node::Int)
	row, col = shape
	((row > 2) && (col > 2)) || throw(ArgumentError("PEPS graph with min size less than 2 is not allowed"))
	n2 = CartesianIndices(shape)[node]
	i, j = n2[1], n2[2]
	idx = LinearIndices(shape)
	return idx[i, mod1(j-1, col)], idx[mod1(i-1, row), j], idx[i, mod1(j+1, col)], idx[mod1(i+1, row), j]
end

function square_edges(shape::Tuple{Int, Int})
	s1, s2 = shape
	index = CartesianIndices(shape)
	idx = LinearIndices(shape)
	_edges = Pair{Int, Int}[]
	for (node, nj) in enumerate(index)
		row, col = nj[1], nj[2]
		node_down = idx[mod1(row+1, s1), col]
		node_right = idx[row, mod1(col+1, s2)]
		push!(_edges, node=>node_down)
		push!(_edges, node_down=>node)
		push!(_edges, node=>node_right)
		push!(_edges, node_right=>node)
	end
	return _edges
end