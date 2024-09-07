push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON


gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table_iv_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"


fidelity(x::AbstractMatrix, y::AbstractMatrix) = real(tr(sqrt(x) * sqrt(y)))
trace_distance(x::AbstractMatrix, y::AbstractMatrix) = 0.5*tr(abs.(x - y))


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

function check_local_rmd2(D2::Int)
	m = 10
	n = 10
	block_size = (5, 5)
	peps_path = gen_peps_path(m, n, D2, block_size)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=block_size, update_alg=alg1)
	alg3 = SimpleUpdate(D2=D2)


	rho1s = rdm2s(peps, alg1)
	rho2s = rdm2s(peps, alg2)

	su_nsweeps = 100
	println("trivial iteration for $su_nsweeps sweeps")
	@time su_peps = get_su_peps(peps, D2, su_nsweeps)
	rho3s = rdm2s(su_peps, alg3)

	h = heisenberg2D(m, n, periodic=false)

	energy1 = expectation(h, peps, alg1) / prod(size(peps))
	energy2 = expectation(h, su_peps, alg1) / prod(size(peps))
	energy3 = expectation(h, peps, alg2) / prod(size(peps))
	println("check energy: E1=$energy1, E2=$energy2, E3=$energy3")

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


	results = Dict("f_bp"=>f1s, "f_su"=>f2s, "d_bp"=>d1s, "d_su"=>d2s)

	result_path = "result/heisenberg_bp_rdm2_D2$(D2)_D1$(D1).json"

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end

	return f1s, f2s
end
