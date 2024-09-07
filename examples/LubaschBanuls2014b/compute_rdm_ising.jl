push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using QuantumSpins, PEPSTools
using Serialization, JSON


gen_peps_path(D2::Int, B::Real) = "data/Table_vii_D2$(D2)_B$(B).peps"

fidelity(x::AbstractMatrix, y::AbstractMatrix) = real(tr(sqrt(x) * sqrt(y)))

trace_distance(x::AbstractMatrix, y::AbstractMatrix) = 0.5*tr(abs.(x - y))

# function make_positive(N2::AbstractArray{<:Number, 2})
# 	N2H = (N2 + N2') / 2
# 	if real(tr(N2H)) < 0.
# 		N2H = -N2H
# 	end
# 	eigvalues, eigvectors = eigen(Hermitian(N2H))
# 	# println("eigenvalues $eigvalues")

# 	eigvalue_threshold = 1.0e-12
# 	pos = findfirst(x->x >= eigvalue_threshold, eigvalues / norm(eigvalues))
# 	isnothing(pos) && error("environment becomes trivial")
# 	# println(eigvalues[1:pos-1])
# 	# pos = 1

# 	eigvalues = eigvalues[pos:end]
# 	# println("eigenvalues $eigvalues")
# 	eigvectors = eigvectors[:, pos:end]
# 	X = eigvectors * Diagonal(sqrt.(eigvalues))
# 	r = X * X'
# 	# r ./= tr(r)
# 	err = sum(eigvalues[1:pos-1])
# 	return r, err
# end

# function make_positive(m::AbstractArray{<:Number, 4})
# 	m2, err = make_positive(reshape(m, (4,4)))
# 	return reshape(m2, (2,2,2,2)), err
# end 

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


function compute_rmd2s(D2::Int, B::Real)
	peps_path = gen_peps_path(D2, B)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10, mult_alg=IterativeCompression(D=D1, maxiter=10, tol=1.0e-6))
	alg2 = BlockBP(block_size=(7,7), update_alg=alg1)
	alg3 = SimpleUpdate(D2=D2)


	rho1s = rdm2s(peps, alg1)
	rho2s = rdm2s(peps, alg2)

	su_nsweeps = 100
	println("trivial iteration for $su_nsweeps sweeps")
	@time su_peps = get_su_peps(peps, D2, su_nsweeps)
	rho3s = rdm2s(su_peps, alg3)

	# h = ising2D(21, 21, J = -1, hz = -B, periodic=false)

	# energy1 = expectation(h, peps, alg1) / prod(size(peps))
	# energy2 = expectation(h, su_peps, alg1) / prod(size(peps))
	# energy3 = expectation(h, peps, alg2) / prod(size(peps))
	# println("check energy: E1=$energy1, E2=$energy2, E3=$energy3")



	rho1s = rho1s["H"]
	rho2s = rho2s["H"]
	rho3s = rho3s["H"]


	rho1s = [reshape(item, (4, 4)) for item in rho1s]
	rho2s = [reshape(item, (4, 4)) for item in rho2s]
	rho3s = [reshape(item, (4, 4)) for item in rho3s]


	f1s = [fidelity(a, b) for (a, b) in zip(rho1s, rho2s)]
	f2s = [fidelity(a, b) for (a, b) in zip(rho1s, rho3s)]

	d1s = [trace_distance(a, b) for (a, b) in zip(rho1s, rho2s)]
	d2s = [trace_distance(a, b) for (a, b) in zip(rho1s, rho3s)]

	# println("quantum fidelity on site ($i, $j) is $f1 and $f2")

	results = Dict("f_bp"=>f1s, "f_su"=>f2s, "d_bp"=>d1s, "d_su"=>d2s)

	result_path = "result/ising_rdm2_D2$(D2)_D1$(D1)_B$(B).json"

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end

	return rho1s, rho2s, rho3s
end


