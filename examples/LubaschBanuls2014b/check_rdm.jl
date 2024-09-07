push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON


gen_peps_path(D2::Int, B::Real) = "data/Table_vii_D2$(D2)_B$(B).peps"

fidelity(x::AbstractMatrix, y::AbstractMatrix) = real(tr(sqrt(x) * sqrt(y)))

trace_distance(x::AbstractMatrix, y::AbstractMatrix) = 0.5*tr(abs.(x - y))

function make_positive(N2::AbstractArray{<:Number, 2})
	N2H = (N2 + N2') / 2
	if real(tr(N2H)) < 0.
		N2H = -N2H
	end
	eigvalues, eigvectors = eigen(Hermitian(N2H))
	# println("eigenvalues $eigvalues")

	eigvalue_threshold = 1.0e-12
	pos = findfirst(x->x >= eigvalue_threshold, eigvalues / norm(eigvalues))
	isnothing(pos) && error("environment becomes trivial")
	# println(eigvalues[1:pos-1])
	# pos = 1

	eigvalues = eigvalues[pos:end]
	# println("eigenvalues $eigvalues")
	eigvectors = eigvectors[:, pos:end]
	X = eigvectors * Diagonal(sqrt.(eigvalues))
	r = X * X'
	r ./= tr(r)
	err = sum(eigvalues[1:pos-1])
	return r, err
end

function make_positive(m::AbstractArray{<:Number, 4})
	m2, err = make_positive(reshape(m, (4,4)))
	# println("trace is $(tr(m2))")
	return reshape(m2, (2,2,2,2)), err
end 

function trivial_operator(m::Int, n::Int)
	iden4 = reshape(one(zeros(4, 4)), (2,2,2,2))
	H = Matrix{Union{Array{Float64, 4}, Nothing}}(nothing, m, n)
	V = Matrix{Union{Array{Float64, 4}, Nothing}}(nothing, m, n)
	for i in 1:m
		for j in 1:n-1
			H[i, j] = iden4
		end
	end

	for i in 1:m-1
		for j in 1:n
			V[i, j] = iden4
		end
	end
	return PEPSTools.SquareLattice(V=V, H=H)
end

function get_su_peps(peps::PEPS, D2::Int, nsweeps::Int=1000)
	m, n = size(peps)
	U = trivial_operator(m, n)

	alg = SimpleUpdate(D2=D2)
	cpeps = CanonicalPEPS(peps)

	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end

	println("the middle singular vector $(cpeps.Hbonds[10, 10])")

	return PEPS(cpeps)
end

function check_local_rmd1(D2::Int)
	block_size = (7, 7)
	m = 21
	n = 21
	B = 2.5
	peps_path = gen_peps_path(D2, B)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=block_size, update_alg=alg1)
	alg3 = SimpleUpdate(D2=D2)

	rho1s = rdm1s(peps, alg1)
	rho2s = rdm1s(peps, alg2)

	su_nsweeps = 100
	println("trivial iteration for $su_nsweeps sweeps")
	@time su_peps = get_su_peps(peps, D2, su_nsweeps)

	h = ising2D(m, n, J = -1, hz = -B, periodic=false)

	energy1 = expectation(h, peps, alg1)
	energy2 = expectation(h, su_peps, alg1)
	println("energy1=$energy1, energy2=$energy2")



	rho3s = rdm1s(su_peps, alg3)

	# f1 = fidelity(rho1, rho2)
	# f2 = fidelity(rho1, rho3)

	# println("quantum fidelity on site ($i, $j) is $f1 and $f2")

	return rho1s, rho2s, rho3s
end


function check_local_rmd2(D2::Int)
	block_size = (7, 7)
	m = 21
	n = 21
	B = 2.5
	peps_path = gen_peps_path(D2, B)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=block_size, update_alg=alg1)
	# alg3 = SimpleUpdate(D2=D2)


	rho1s = rdm2s(peps, alg1)
	rho2s = rdm2s(peps, alg2)
	# rho3 = rdm2H(peps, i, j, alg3)

	# rho1 = reshape(rho1, (4,4))
	# rho2 = reshape(rho2, (4,4))
	# rho3 = reshape(rho3, (4,4))

	# f1 = fidelity(rho1, rho2)
	# f2 = fidelity(rho1, rho3)

	# println("quantum fidelity on site ($i, $j) is $f1 and $f2")

	rho1s = [reshape(item, (4,4)) for item in rho1s["H"]]
	rho2s = [reshape(item, (4,4)) for item in rho2s["H"]]


	return rho1s, rho2s
end
