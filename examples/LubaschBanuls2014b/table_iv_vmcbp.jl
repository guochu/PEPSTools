include("../../../NNQS/parallel/parallel_nqs.jl")

# @everywhere push!(LOAD_PATH, "../../src")
# @everywhere using PEPSTools
@everywhere include("../../src/includes.jl")

using Random
using Serialization, JSON
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int) = "result/Table_iv_m$(m)_n$(n)_D1$(D1)_D2$(D2).json"

gen_bp_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)_vmcbp" * ".peps"
gen_bp_result_path(m::Int, n::Int, D2::Int, lr::Real, epoches::Int) = "result/Table_iv_m$(m)_n$(n)_D2$(D2)_lr$(lr)_epoches$(epoches)_vmcbp.json"


function main(m::Int, n::Int, Dnew::Int, epoches::Int=100;D2::Int, D1::Int=2*D2^2+10, lr::Real=0.001)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, epoches=$(epoches), learn rate=$(lr)")
	L = m * n

	h = squeeze(heisenberg2D(m, n, periodic=false))

	# Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	bp_peps_path = gen_bp_peps_path(m, n, Dnew)
	if ispath(bp_peps_path)
		println("read initial vmcbp peps from path $bp_peps_path")
		state = Serialization.deserialize(bp_peps_path)
		peps_load_path = bp_peps_path
	else
		peps_path = gen_peps_path(m, n, D2) 
		println("read initial peps from path $peps_path")
		ispath(peps_path) || error("peps path not provided")
		state = Serialization.deserialize(peps_path)

		alg = BoundaryMPS(D1=D1, D2=D2)
		bmps_energy = real(energy(h, state, alg))
		println("bmps energy is $(bmps_energy)")

		state = increase_bond!(state, D=Dnew)
		peps_load_path = peps_path
	end

	if scalartype(state) <: Real
		state = complex(state)
	end
	# state = randompeps(ComplexF64, m, n, d=2, D=D2, periodic=false)

	sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)
	opt = ADAM(lr)


	losses = Float64[]

	x0, re = Flux.destructure(state)

	ham = Heisenberg2D((m, n), J=0.25)

    for i in 1:epoches
        @time train_loss, grad = parallel_energy_and_grad(ham, state, sampler, n_chain_per_rank=1, λ=0)

        Optimise.update!(opt, x0, grad)
        state = re(x0)

        push!(losses, train_loss)
        println("energy at the $i-th step is $(train_loss)")
    end


	println("save peps to path $(bp_peps_path)")
	Serialization.serialize(bp_peps_path, state)


	file_name = gen_bp_result_path(m, n, Dnew, lr, epoches) 

	results = Dict("energies"=>losses, "peps_path"=>peps_load_path)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return losses
end

# all the simulations in table IV

function main_10_10()
	m = 10
	n = 10
	main(m, n, 3, D2=2)
end

function main_14_14()
	m = 14
	n = 14
	# main(m, n, D2=2)
	main(m, n, D2=3)
	# main(m, n, D2=4)
end




