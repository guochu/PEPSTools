push!(LOAD_PATH, "../../src")

using Random
using QuantumSpins, PEPSTools
using Serialization, JSON


gen_peps_path(D2::Int, B::Real) = "data/Table_vii_D2$(D2)_B$(B).peps"



function main(D2::Int, B::Real)
	D1 = 2*D2^2 + 10
	m = 21
	n = 21
	L = m * n
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, B=$(B)")

	peps_path = gen_peps_path(D2, B) 

	peps = Serialization.deserialize(peps_path)

	mh = div(m, 2)
	obs = LocalObservers{Matrix{Float64}}((m, n))
	p = spin_half_matrices()
	sz = p["z"]
	obs[mh, mh] = sz

	alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)	
	mags = local_expectations(obs, peps, alg)

	return mags[mh, mh]
end