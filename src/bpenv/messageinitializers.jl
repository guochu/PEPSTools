# message initializers
random_c_bondmessages(peps::Abstract2DTN) = _init_bond_messages(rand, peps)
unit_c_bondmessages(peps::Abstract2DTN) = _init_bond_messages(ones, peps)
random_q_bondmessages(peps::Abstract2DTN) = _init_bond_messages(rand_qvector, peps)
unit_q_bondmessages(peps::Abstract2DTN) = _init_bond_messages(unit_qvector, peps)

function _init_c_bondmessages(peps::Abstract2DTN, initguess::Symbol, seed::Int)
	if initguess == :rand
		Random.seed!(seed)
		return random_c_bondmessages(peps)
	elseif initguess == :unit
		return unit_c_bondmessages(peps)
	else
		throw(ArgumentError("initguess $(initguess) not implemented"))
	end
end
function _init_q_bondmessages(peps::PEPS, initguess::Symbol, seed::Int)
	if initguess == :rand
		Random.seed!(seed)
		return random_q_bondmessages(peps)
	elseif initguess == :unit
		return unit_q_bondmessages(peps)
	else
		throw(ArgumentError("initguess $(initguess) not implemented"))
	end
end

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

