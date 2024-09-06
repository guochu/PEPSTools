include("../../NNQS/parallel/parallel_nqs.jl")

@everywhere push!(LOAD_PATH, "../src")
@everywhere using PEPSTools

using Random
using Serialization, JSON
using Flux, Flux.Optimise


function main_ising(L; hz=1, J=1, D=2)

	Random.seed!(1234)

	ham = Ising2D((L, L), h=hz, J=J)

	T = ComplexF64

	state = randompeps(T, L, L, d=2, D=D)

	sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)

	learn_rate = 0.01
	epoches = 200
	opt = ADAM(learn_rate)


	# println("initial loss is ", NNQS.energy(ham, state, sampler))

	losses = Float64[]

	x0, re = Flux.destructure(state)

	# bmps 
	alg = BoundaryMPS(D1=2*D^2+10, D2=D, als_maxiter=10)
	h2 = squeeze(ising2D(L, L, hz=hz, J=J))

    for i in 1:epoches
        @time train_loss, grad = parallel_energy_and_grad(ham, state, sampler, n_chain_per_rank=1, λ=1.0e-5)

        # bmps_energy = real(energy(h2, state, alg))

        Optimise.update!(opt, x0, grad)
        state = re(x0)

        push!(losses, train_loss)
        # println("energy at the $i-th step is $(train_loss), bmps energy is $(bmps_energy)")
        println("energy at the $i-th step is $(train_loss)")
    end
    return losses
end


function main_xxx(L; D=2)

	Random.seed!(1234)

	# there is a difference in convention in NNQS and PEPSTools
	ham = Heisenberg2D((L, L), J=0.25)

	T = ComplexF64

	state = randompeps(T, L, L, d=2, D=D)

	sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)

	learn_rate = 0.01
	epoches = 200
	opt = ADAM(learn_rate)


	# println("initial loss is ", NNQS.energy(ham, state, sampler))

	losses = Float64[]

	x0, re = Flux.destructure(state)

	# bmps 
	alg = BoundaryMPS(D1=2*D^2+10, D2=D, als_maxiter=10)
	h2 = squeeze(heisenberg2D(L, L))

    for i in 1:epoches
        @time train_loss, grad = parallel_energy_and_grad(ham, state, sampler, n_chain_per_rank=1, λ=1.0e-5)

        # bmps_energy = real(energy(h2, state, alg))

        Optimise.update!(opt, x0, grad)
        state = re(x0)

        push!(losses, train_loss)
        # println("energy at the $i-th step is $(train_loss), bmps energy is $(bmps_energy)")
        println("energy at the $i-th step is $(train_loss)")
    end
    return losses
end