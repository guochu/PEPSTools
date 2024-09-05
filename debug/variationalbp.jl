include("../src/includes.jl")


using Random
using Serialization, JSON
using Flux, Flux.Optimise



function main_bp(L; B::Real, D=2)

	Random.seed!(2344512)
	
	h = squeeze(ising2D(L, L, hz=-B, J=-1))

	T = ComplexF64

	state = randompeps(T, L, L, d=2, D=D)

	env = bp_environments(state, 5, verbosity=1)

	bp_tol = 1.0e-8

	# bp update

	α = 1.0e-1
	println("start variational BP with α=$(α)")
	for i in 1:10
	
		E, grad = withgradient(energy, h, state, env)
		axpy!(-α, grad[2], state) 
		println("energy at $(i)-th iteration is $(E/length(state)) with α=$(α)")

		env = bp_environments(state, 5, bp_tol, verbosity=1)
	end
	println()

	α = 1.0e-2
	println("start variational BP with α=$(α)")
	for i in 1:100

		E, grad = withgradient(energy, h, state, env)
		axpy!(-α, grad[2], state) 
		println("energy at $(i)-th iteration is $(E/length(state)) with α=$(α)")

		env = bp_environments(state, 15, bp_tol, verbosity=1)
	end
	println()





	# learn_rate = 0.005
	# epoches = 1000
	# opt = ADAM(learn_rate)

	# # x0, re = Flux.destructure(rbm)

	# losses = Float64[]

	# x0, re = Flux.destructure(state)

	# sampler = MetropolisLocal(length(ham), n_thermal=100, n_sample_per_chain=500, n_discard=10)
	# println("start NNQS BP...")
	# for i in 1:epoches
	
	# 	E1, grad = withgradient(expectation, h, state, env)

	# 	@time E2 = energy(ham, state, sampler, n_chain=10)

	# 	# axpy!(-α, grad[2], state) 

	# 	grad, _ = Flux.destructure(grad[2])

	# 	Optimise.update!(opt, x0, grad)
	# 	state = re(x0)


	# 	println("$(i)-th iteration: BP energy $(E1/length(h)), VMC energy $(E2/length(h))")

	# 	env = environments(state, 5, bp_tol, verbosity=1)
	# end
	# println()



end