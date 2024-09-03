

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

# message initializers
random_c_bondmessages(peps::Abstract2DTN) = _init_bond_messages(rand, peps)
unit_c_bondmessages(peps::Abstract2DTN) = _init_bond_messages(ones, peps)
random_q_bondmessages(peps::Abstract2DTN) = _init_bond_messages(rand_qvector, peps)
unit_q_bondmessages(peps::Abstract2DTN) = _init_bond_messages(unit_qvector, peps)


function _init_bond_messages(f, peps::Abstract2DTN)
	T = scalartype(peps)
	Hbonds = [_random_vector_message(T, f, size(item, 3)) for item in peps.data]
	Vbonds = [_random_vector_message(T, f, size(item, 4)) for item in peps.data]
	return SquareLatticeBonds(H=Hbonds, V=Vbonds)
end

function rand_qvector(::Type{T}, D::Int) where {T<:Number}
	m = rand(T, D, D)
	m = m' * m
	return reshape(m, D*D)
end
function unit_qvector(::Type{T}, D::Int) where {T<:Number}
	m = zeros(T, D, D)
	for i in 1:D
		m[i, i] = 1
	end
	return reshape(m, D*D)
end

_random_vector_message(::Type{T}, f, D::Int) where {T <: Number} = Message(f(T, D), f(T, D))

