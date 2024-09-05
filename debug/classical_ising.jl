include("../src/includes.jl")

using Random

function finite_ising(β)

	m = 20
	n = 20

	i = div(m, 2)
	j = div(n, 2)

	block_size = (5,5)

	D2 = 2

	D1 = 2*D2^2 + 10

	Random.seed!(1234)

	model = ClassicalIsing2D(m, n, J=1, hz=1.0e-3, periodic=false)


	# # peps = randomdoublelayerpeps(m, n, D=D2, periodic=true)
	# peps = PEPSTools.periodic_classical_ising_peps(m, n, J=-1, hz=-1.0e-3, β=β)

	# observer = LocalObservers(m, n, sz)

	obs1 = magnetization(model, i, j, BoundaryMPS(D2=D2, D1=D1), β=β)


	# pos = (1,1)

	# v = exact_contract(peps, pos[1], pos[2], sz)


	update_alg = BoundaryMPS(D2=D2,D1=D1,mult_alg=IterativeCompression(D=D1))
	

	obs2 = magnetization(model, i, j, BlockBp(msg_maxiter=30, msg_tol=1.0e-8, block_size=block_size, update_alg=update_alg), β=β)

	return obs1, obs2
	# return obs2, obs3

	# return obs1[pos...], v
end

function Onsager_m(beta)
	
	beta_C =  log(1+sqrt(2.0))/2
	
	if beta<= beta_C
		exact=0.0
	else
		exact = ( 1-1/(sinh(2*beta))^4)^(1.0/8)
	end
	return exact

end

function infinite_ising()

	m = 4
	i = j = div(m+1, 2)

	J = 1.0

	D2 = 2
	# D1 = 2*D2^2 + 10
	D1 = 8
	update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=50,mult_alg=IterativeCompression(D=D1,maxiter=50))
	exp_alg = BlockBP(msg_maxiter=50, msg_tol=1.0e-11, msg_D=6, block_size=(m, m), update_alg=update_alg)

	bp_alg = BP(msg_maxiter=50, msg_tol=1.0e-11)

	h = 0.
	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=true)


	betas_list = collect(0.2:0.05:1)

	tol = 1.0e-2

	for beta in betas_list
		obs1 = abs(magnetization(model, i, j, exp_alg, β=beta))
		obs2 = abs(magnetization(model, i, j, bp_alg, β=beta))
		obs_exact = Onsager_m(beta)

		# @test abs(obs1 - obs_exact) < tol
		println("beta=$beta, BlockBP=$(obs1), BP=$(obs2) Onsager=$(obs_exact)")
	end


end