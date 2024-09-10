# @everywhere push!(LOAD_PATH, "../../src")
# @everywhere using PEPSTools
include("../../src/includes.jl")

using Random
using Serialization, JSON
using Flux.Optimise
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int) = "result/Table_iv_m$(m)_n$(n)_D1$(D1)_D2$(D2).json"

gen_bp_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)_vmcbp" * ".peps"
gen_bp_result_path(m::Int, n::Int, D2::Int, lr::Real, epoches::Int) = "result/Table_iv_m$(m)_n$(n)_D2$(D2)_lr$(lr)_epoches$(epoches)_vmcbp.json"


function main(m::Int, n::Int;D2::Int, D1::Int=2*D2^2+10)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2")
	L = m * n

	h = squeeze(heisenberg2D(m, n))

	# Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	peps_path = gen_peps_path(m, n, D2) 
	println("read initial peps from path $peps_path")
	ispath(peps_path) || error("peps path not provided")
	state = Serialization.deserialize(peps_path)

	alg = BoundaryMPS(D1=D1, D2=D2)
	bmps_energy = real(energy(h, state, alg))
	println("bmps energy is $(bmps_energy/L)")

	bp_energy = energy(h, state, bp_environments(state, 10, 1.0e-8, verbosity=1))
	println("bp energy is $(bp_energy/L)")

	if scalartype(state) <: Real
		state = complex(state)
	end
	# state = randompeps(ComplexF64, m, n, d=2, D=D2, periodic=false)

	sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=1000, n_discard=10)

	ham = Heisenberg2D((m, n), J=0.25)

	vmc_bp_energy = NNQS.energy(ham, state, sampler)
	println("vmc bp energy is $(vmc_bp_energy/L)")


end






