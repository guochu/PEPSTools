push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON

gen_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)" * ".peps"



function compute_energy_bp(m::Int, D2::Int)
	peps_path = gen_peps_path(m, m, D2)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	mh = div(m, 2)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	alg2 = BlockBP(block_size=(mh,mh), update_alg=alg1)

	m, n = size(peps)
	h = heisenberg2D(m, n, periodic=false)

	energy2 = expectation(h, peps, alg2)

	return energy2 / (m*n)
end

