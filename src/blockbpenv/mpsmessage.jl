const MPSMessage{M} = Message{M, M} where {M <: MPS}

function random_mps_message(::Type{T}, physpaces::Vector{Int}; Di::Int=maximum(physpaces), Do::Int=Di) where {T <: Number}
	i = random_boundary_mps(T, physpaces, D=Di)
	# canonicalize!(i, normalize=true)
	# rightorth!(i, alg=QRFact())
	# normalize!(i, iscanonical=true)
	rightorth!(i, alg=Orthogonalize(alg=QR(), normalize=true))

	o = random_boundary_mps(T, physpaces, D=Do)
	# canonicalize!(o, normalize=true)
	# rightorth!(o, alg=QRFact())
	# normalize!(o, iscanonical=true)
	rightorth!(o, alg=Orthogonalize(alg=QR(), normalize=true))
	return Message(i, o)
	
end
random_mps_message(::Type{T}, L::Int; d::Int, kwargs...) where {T <: Number} = random_mps_message(T, [d for i in 1:L]; kwargs...)
trivial_mps_message(::Type{T}, L::Int) where {T <: Number} = Message(trivial_mps(T, L), trivial_mps(T, L))

# message_distance2(x::Message, y::Message) = distance2(x.i, y.i) + distance2(x.o, y.o)

function message_distance2(x::MPSMessage, y::MPSMessage)
	@assert physical_dimensions(x.i) == physical_dimensions(y.i)
	@assert physical_dimensions(x.o) == physical_dimensions(y.o)
	return distance2(x.i, y.i) + distance2(x.o, y.o)
end

# PEPSBlock contains all the input messages
function compute_out_messages(blk::AbstractBorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	# compute left output message
	mpsl = right_boundary(blk)
	for j in n:-1:1
		mpoj = mporight(blk, j)
		mpsl, err = mpompsmult(mpoj, mpsl, alg)
		setscaling!(mpsl, 1)
		(alg.verbosity >= 3) && println("error for computing the left out message at the $j-th column is $err")
	end
	mpsl = _mps_rm_boundary(mpsl)

	# compute the right output message
	mpsr = left_boundary(blk)
	for j in 1:n
		mpoj = mpoleft(blk, j)
		mpsr, err = mpompsmult(mpoj, mpsr, alg)
		setscaling!(mpsr, 1)
		(alg.verbosity >= 3) && println("error for computing the right out message at the $j-th column is $err")
	end
	mpsr = _mps_rm_boundary(mpsr)

	# compute the up output message
	mpsu = down_boundary(blk)
	for j in m:-1:1
		mpoj = mpodown(blk, j)
		mpsu, err = mpompsmult(mpoj, mpsu, alg)
		setscaling!(mpsu, 1)
		(alg.verbosity >= 3) && println("error for computing the up out message at the $j-th row is $err")
	end
	mpsu = _mps_rm_boundary(mpsu)

	# compute the down output message
	mpsd = up_boundary(blk)
	for j in 1:m
		mpoj = mpoup(blk, j)
		mpsd, err = mpompsmult(mpoj, mpsd, alg)
		setscaling!(mpsd, 1)
		(alg.verbosity >= 3) && println("error for computing the down out message at the $j-th row is $err")
	end
	mpsd = _mps_rm_boundary(mpsd)

	# canonicalize!(mpsl, normalize=true)
	# canonicalize!(mpsr, normalize=true)
	# canonicalize!(mpsu, normalize=true)
	# canonicalize!(mpsd, normalize=true)
	return mpsl, mpsr, mpsu, mpsd
end	

function _mps_rm_boundary(mps)
	L = length(mps)
	v = mps[2:L-1]
	# @assert length(mps[1]) == length(mps[L]) == 1
	v[end] .*= only(mps[L])
	v[1] .*= only(mps[1])
	v[1] ./= norm(v[1])
	r = MPS(v, scaling=scaling(mps))
	# @assert norm(r, iscanonical=true) ≈ 1.
	return r
end

