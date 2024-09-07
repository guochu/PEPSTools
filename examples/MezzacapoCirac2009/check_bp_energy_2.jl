push!(LOAD_PATH, "../../src")

using Random
using PEPSTools, PeriodicMPS
using Serialization, JSON


gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"



function main(D2::Int, k::Int)
	m = 4
	n = m
	D1 = 2*D2^2 + 10

	block_size = (m, n)

	peps_path = gen_peps_path(m, n, D2, block_size)
	peps4 = Serialization.deserialize(peps_path)

	h4 = heisenberg2D(m, n, periodic=true)


	peps = repeat(peps4, k, k)
	h = repeat(h4, k, k)

	alg = BlockBP(block_size=size(peps), msg_D=D1, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 

	energy = expectation(h, peps, alg) / prod(size(peps))


	return energy
end
