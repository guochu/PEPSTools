# @everywhere push!(LOAD_PATH, "../../src")
# @everywhere using PEPSTools
push!(LOAD_PATH, "../../src")

using Random
using Serialization, JSON
using PEPSTools, NNQS
using Flux.Optimise
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int, B::Real) = "data/Table_vi_m$(m)_n$(n)_D2$(D2)_B$(B).peps"

function main(m::Int, n::Int;D2::Int, D1::Int=2*D2^2, B::Real)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, B=$B")
	L = m * n

	h = squeeze(ising2D(m, n, J = -1, hz = -B, periodic=false))

	# Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	peps_path = gen_peps_path(m, n, D2, B) 
	println("read initial peps from path $peps_path")
	ispath(peps_path) || error("peps path not provided")
	state = Serialization.deserialize(peps_path)

	alg = BoundaryMPS(D1=D1, D2=D2)
	bmps_energy = real(energy(h, state, alg))
	println("bmps energy is $(bmps_energy/L)")

	bp_alg = BP(FixedNorm(), msg_maxiter=100, msg_tol=1.0e-8, verbosity=1, damping=0)

	bp_energy = energy(h, state, bp_alg)
	println("bp energy is $(bp_energy/L)")

	# if scalartype(state) <: Real
	# 	state = complex(state)
	# end
	# state = randompeps(ComplexF64, m, n, d=2, D=D2, periodic=false)

	sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=1000, n_discard=10)

	ham = Ising2D((m, n), h=-B, J=-1)

	vmc_bp_energy = energy(ham, state, sampler)
	println("vmc bp energy is $(vmc_bp_energy/L)")

end


function main_2(m::Int, n::Int;D2::Int, D1::Int=2*D2^2+10)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2")
	L = m * n

	# Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	peps_path = gen_peps_path(m, n, D2, B) 
	println("read initial peps from path $peps_path")
	ispath(peps_path) || error("peps path not provided")
	state = Serialization.deserialize(peps_path)


	basis = rand((-1, 1), length(state))

	normalize_alg = FixedSum()

	alg = BP(normalize_alg, msg_maxiter=100, msg_tol=1.0e-8, verbosity=1, damping=0.2)

	bp_amp = Ψ(state, basis, alg)

	# msg, converged = fixedpoint_messages(tn, msg, alg)
	# bp_amp = _Ψ_util(state, basis, msg, alg)


	exact_amp = exact_amplitude(state, _state_to_index(basis))

	println("exact ", exact_amp, ", approximate ", bp_amp, ", error ", abs((exact_amp - bp_amp) / exact_amp) )

end






