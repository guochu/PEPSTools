push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON


gen_peps_path(D2::Int, B::Real) = "data/Table_vii_D2$(D2)_B$(B).peps"


function compute_energies(D2::Int, B::Real)
	peps_path = gen_peps_path(D2, B)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=(7,7), update_alg=alg1)

	m, n = size(peps)
	h = ising2D(m, n, J = -1, hz = -B, periodic=false)

	energy1 = expectation(h, peps, alg1)
	energy2 = expectation(h, peps, alg2)

	return energy1 / (m*n), energy2 / (m*n)
end


function compute_energy_bp(D2::Int, B::Real)
	peps_path = gen_peps_path(D2, B)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	alg2 = BlockBP(block_size=(7,7), update_alg=alg1)

	m, n = size(peps)
	h = ising2D(m, n, J = -1, hz = -B, periodic=false)

	energy2 = expectation(h, peps, alg2)

	return energy2 / (m*n)
end

