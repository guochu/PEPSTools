

function LinearAlgebra.normalize!(m::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm)
	for item in m.H
		_normalize!(item, alg)
	end
	for item in m.V
		_normalize!(item, alg)
	end
	return m
end

function canonicalize!(m::SquareLatticeBondMessages)
	_canonicalize!.(m.H)
	_canonicalize!.(m.V)
	return m
end

function _normalize!(m::VectorMessage, alg::MessageNormalizationAlgorithm)
	_normalize!(m.i, alg)
	_normalize!(m.o, alg)
	return m
end

_normalize!(v::Vector, alg::FixedNorm) = normalize!(v)
function _normalize!(v::Vector, alg::FixedSum)
	v ./= sum(v)
	return v
end
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

