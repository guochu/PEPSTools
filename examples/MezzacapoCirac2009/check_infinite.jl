push!(LOAD_PATH, "../../src")

using Random
using QuantumSpins, PEPSTools
using Serialization, JSON


function main(block_size::Tuple{Int, Int})
	B = 3.
	D2 = 2
	D1 = 2*D2^2+10

	m = 2
	n = 2

	h = ising2D(m, n, J = -1, hz = -B, periodic=true)

	Random.seed!(3598)
	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	peps = randompeps(Float64, m, n, d=2, D=1)


	# # blks = default_splitting(size(peps), block_size)

	# blks = center_splitting(size(peps), block_size, (2,2))

	# println(PEPSTools.is_valid_hamiltonian_splitting(size(peps), blks))

	# return blks


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=D1, update_alg=exp_alg)
	measure_alg = bp_alg


	U = exponential(h, -0.001)
	# Us = center_splitting(U, bp_alg.block_size)

	# Us = default_splitting(U, bp_alg.block_size)
	# return Us


	# gs_alg = ImaginaryTimePEPS([(1000, -0.001, 1.0e-8)], stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)

	for i in 1:50
		for j in 1:100
			sweep!(peps, U, bp_alg)
		end
		energy = expectation(h, peps, measure_alg) / prod(size(peps))
		println("energy is $energy")
	end

end
