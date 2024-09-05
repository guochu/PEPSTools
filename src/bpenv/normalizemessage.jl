function LinearAlgebra.normalize!(m::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm)
	for item in m.H
		normalize!(item, alg)
	end
	for item in m.V
		normalize!(item, alg)
	end
	return m
end
function LinearAlgebra.normalize(m::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm)
	V, H = get_data(m)
	return SquareLatticeBonds([_normalize(x, alg) for x in V], [_normalize(x, alg) for x in H])
end 

function canonicalize!(m::SquareLatticeBondMessages)
	_canonicalize!.(m.H)
	_canonicalize!.(m.V)
	return m
end

LinearAlgebra.normalize(m::VectorMessage, alg::MessageNormalizationAlgorithm) = normalize!(copy(m), alg)

function LinearAlgebra.normalize!(m::VectorMessage, alg::MessageNormalizationAlgorithm)
	_normalize!(m.i, alg)
	_normalize!(m.o, alg)
	return m
end

_normalize!(v::Vector, alg::FixedNorm) = normalize!(v)
function _normalize!(v::Vector, alg::FixedSum)
	v ./= sum(v)
	return v
end

function canonicalize(m::SquareLatticeBondMessages)
	V, H = get_data(m)
	return SquareLatticeBonds(_canonicalize.(V), _canonicalize.(H))
end 
# function canonicalize_threaded(m::SquareLatticeBondMessages)
# 	V, H = get_data(m)
# 	V′, H′ = similar(V), similar(H)
# 	fetch.([Threads.@spawn V′[i] = _canonicalize(V[i]) for i in 1:length(V)])
# 	fetch.([Threads.@spawn H′[i] = _canonicalize(H[i]) for i in 1:length(H)])
# 	return SquareLatticeBonds(V′, H′)
# end
function _canonicalize!(m::VectorMessage)
	sqrt_z = sqrt(mapreduce(*, +, m.i, m.o))
	m.i ./= sqrt_z
	m.o ./= sqrt_z
	return m
end

function message_distance2(x::SquareLatticeBondMessages, y::SquareLatticeBondMessages) 
	# f(a, b) = message_distance2(a, b)
	_zero = zero(real(scalartype(x)))
	return mapreduce(message_distance2, +, x.H, y.H, init=_zero) + mapreduce(message_distance2, +, x.V, y.V, init=_zero)
end
message_distance2(x::VectorMessage, y::VectorMessage) = distance2(x.i, y.i) + distance2(x.o, y.o)


# Automatic differentiation
Zygote.@adjoint Message(i, o) = Message(i, o), z -> (z.i, z.o)
Zygote.@adjoint SquareLatticeBonds(V::AbstractMatrix, H::AbstractMatrix) = SquareLatticeBonds(V, H), z -> (z.V, z.H)
get_data(x::SquareLatticeBonds) = get_data(x.V), get_data(x.H)
Zygote.@adjoint get_data(x::SquareLatticeBonds) = get_data(x), z -> (SquareLatticeBonds(z[1], z[2]),)
get_data(x::Message) = x.i, x.o
Zygote.@adjoint get_data(x::Message) = get_data(x), z -> (Message(z[1], z[2]),)

Zygote.@adjoint PeriodicArray(data::Array) = PeriodicArray(data), z -> (z,)
get_data(x::PeriodicArray) = x.data
Zygote.@adjoint get_data(x::PeriodicArray) = get_data(x), z -> (z,)



_normalize(v::Vector, alg::FixedSum) = _normalize_sum(v)
_normalize(v::Vector, alg::FixedNorm) = normalize(v)
_normalize_sum(v::Vector) = v / sum(v)


function _normalize(v::VectorMessage, alg::MessageNormalizationAlgorithm)
	i, o = get_data(v)
	return Message(_normalize(i, alg), _normalize(o, alg))
end 

function _canonicalize(m::VectorMessage)
	ab, ba = get_data(m)
	sqrt_z = sqrt(mapreduce(*, +, ab, ba))
	return Message(ab / sqrt_z, ba / sqrt_z)
end