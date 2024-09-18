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
	# sqrt_z = sqrt(mapreduce(*, +, m.i, m.o))
	# m.i ./= sqrt_z
	# m.o ./= sqrt_z

	z = mapreduce(*, +, m.i, m.o)
	m.o ./= z
	return m
end

function mixing(m_out::SquareLatticeBondMessages, m_in::SquareLatticeBondMessages, damping::Real)
	m_out_V, m_out_H = get_data(m_out)
	m_in_V, m_in_H = get_data(m_in)
	V = [mixing(o, i, damping) for (o, i) in zip(m_out_V, m_in_V)]
	H = [mixing(o, i, damping) for (o, i) in zip(m_out_H, m_in_H)]
	return SquareLatticeBonds(V, H)
end
function mixing!(m_out::SquareLatticeBondMessages, m_in::SquareLatticeBondMessages, damping::Real)
	for (o, i) in zip(m_out.V, m_in.V)
		mixing!(o, i, damping)
	end
	for (o, i) in zip(m_out.H, m_in.H)
		mixing!(o, i, damping)
	end	
	return m_out
end

function mixing!(m_out::VectorMessage, m_in::VectorMessage, damping::Real)
	mixing!(m_out.i, m_in.i, damping)
	mixing!(m_out.o, m_in.o, damping)
	return m_out
end

function mixing(m_out::VectorMessage, m_in::VectorMessage, damping::Real)
	m_out_i, m_out_o = get_data(m_out)
	m_in_i, m_in_o = get_data(m_in)
	i = mixing(m_out_i, m_in_i, damping)
	o = mixing(m_out_o, m_in_o, damping)
	return Message(i, o)
end

function mixing!(m_out::AbstractVector, m_in::AbstractVector, damping::Real)
	if damping > 0
		axpby!(damping, m_in, 1 - damping, m_out)
	end
	return m_out
end
mixing(m_out::AbstractVector, m_in::AbstractVector, damping::Real) = (1 - damping) .* m_out .+ damping .* m_in

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


function Base.:+(x::VectorMessage, y::VectorMessage)
	xi, xo = get_data(x)
	yi, yo = get_data(y)
	return Message(xi + yi, xo + yo)
end

function Base.:+(x::SquareLatticeBondMessages, y::SquareLatticeBondMessages)
	xV, xH = get_data(x)
	yV, yH = get_data(y)
	return SquareLatticeBonds(xV .+ yV, xH .+ yH)
end

_normalize(v::Vector, alg::FixedSum) = _normalize_sum(v)
_normalize(v::Vector, alg::FixedNorm) = normalize(v)
_normalize_sum(v::Vector) = v / sum(v)


function _normalize(v::VectorMessage, alg::MessageNormalizationAlgorithm)
	i, o = get_data(v)
	return Message(_normalize(i, alg), _normalize(o, alg))
end 

function _canonicalize(m::VectorMessage)
	# ab, ba = get_data(m)
	# sqrt_z = sqrt(mapreduce(*, +, ab, ba))
	# return Message(ab / sqrt_z, ba / sqrt_z)

	ab, ba = get_data(m)
	z = mapreduce(*, +, ab, ba)
	return Message(ab, ba / z)
end
