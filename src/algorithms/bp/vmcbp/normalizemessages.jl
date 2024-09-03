# out of place version for message normalization


LinearAlgebra.normalize(m::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm) = normalize!(copy(m), alg)

Zygote.@adjoint LinearAlgebra.normalize(m::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm) = begin
	out = similar(m)
	backs = Dict{Pair{Int, Int}, Any}()
	for (k, v) in m.data
		out[k], backs[k] = Zygote.pullback(normalize, v)
	end
	return GraphMessage(m.graph, out), z -> begin
		data_back = Dict(k=>backs[k](v)[1] for (k, v) in z)
		return data_back, nothing
	end
end


normalize_sum(m::GraphMessage) = GraphMessage(m.graph, Dict(k=>v/sum(v) for (k, v) in m.data))

LinearAlgebra.normalize(m::GraphMessage, alg::FixedSum) = GraphMessage(m.graph, Dict(k=>v/sum(v) for (k, v) in m.data))
Zygote.@adjoint LinearAlgebra.normalize(m::GraphMessage, alg::FixedSum) = begin
	out = typeof(m.data)()
	backs = Dict{Pair{Int, Int}, Any}()
	for (k, v) in m.data
		out[k], backs[k] = Zygote.pullback(x->x/sum(x), v)
	end
	return GraphMessage(m.graph, out), z -> begin
		data_back = Dict(k=>backs[k](v)[1] for (k, v) in z)
		return data_back, nothing
	end
end


canonicalize(m::SquareLatticeBondMessages) = canonicalize!(copy(m))

Zygote.@adjoint canonicalize(m::SquareLatticeBondMessages) = begin
	out = similar(m)
	backs = Dict{Pair{Int, Int}, Any}()
	for edge in edges(m.graph)
		_edge = edge.src=>edge.dst
		r, backs[_edge] = Zygote.pullback(_canonicalize, m[_edge], m[reverse(_edge)])
		out[_edge] = r[1]
		out[reverse(_edge)] = r[2]
	end
	return out, z -> begin
		data_back = similar(m)
		for edge in edges(m.graph)
			_edge = edge.src=>edge.dst
			data_back[_edge], data_back[reverse(_edge)] = backs[_edge]((z[_edge], z[reverse(_edge)]))
		end
		return (data_back,)
	end		
end

function _canonicalize(ab::Vector, ba::Vector)
	sqrt_z = sqrt(mapreduce(*, +, ab, ba))
	return ab / sqrt_z, ba / sqrt_z
end