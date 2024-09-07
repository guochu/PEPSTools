push!(LOAD_PATH, "../../src")

using Random, JSON
using QuantumSpins, PEPSTools

"""
	
	Returns Onsager's exact solution for the magnetization at inverse
	temperature beta, for an infinite and homogeneous system with J=1
	and h=0
	
"""
function Onsager_m(beta)
	
	beta_C =  log(1+sqrt(2.0))/2
	
	if beta<= beta_C
		exact=0.0
	else
		exact = ( 1-1/(sinh(2*beta))^4)^(1.0/8)
	end
	return exact

end


function main_periodic(m::Int; block_size::Int = m)
	i = j = div(m+1, 2)

	Random.seed!(34291)

	J = 1.0

	D2 = 2
	# D1 = 2*D2^2 + 10
	D1 = 8
	update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=50,mult_alg=IterativeCompression(D=D1,maxiter=50))
	exp_alg = BlockBP(msg_maxiter=50, msg_tol=1.0e-11, msg_D=6, block_size=(m, m), update_alg=update_alg)

	h = 0.
	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=true)
	model_obc = ClassicalIsing2D(m, m, J=J, hz=1.0e-3, periodic=false)


	betas_list = vcat(collect(0.2:0.05:0.4), collect(0.405:0.005:0.495), collect(0.5:0.05:1))


	bp_results = Float64[]
	Onsager_results = Float64[]
	boundary_mps_results = Float64[]

	for beta in betas_list
		obs1 = magnetization(model_obc, i, j, update_alg, β=beta)
		obs2 = magnetization(model, i, j, exp_alg, β=beta)
		m_Onsager = Onsager_m(beta)

		push!(boundary_mps_results, obs1)
		push!(bp_results, obs2)
		push!(Onsager_results, m_Onsager)

		println("beta=$beta, BoundaryMPS=$obs1, BP=$obs2, Onsager=$(m_Onsager)")
	end


	data_path = "result/infinite_ising_tn_N$(m)_h$(h)_blocksize$(block_size)_D$(D1).json"
	results = Dict("bmps"=>boundary_mps_results, "bp"=>bp_results, "Onsager"=>Onsager_results, "betas"=>betas_list)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

end

function main_periodic_critical(m::Int; block_size::Int = m)
	i = j = div(m+1, 2)

	J = 1.0

	D2 = 2
	D1 = 2*D2^2 
	update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=20,mult_alg=IterativeCompression(D=D1,maxiter=20))
	# update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=20)
	exp_alg = BlockBP(msg_maxiter=50, msg_tol=1.0e-10, msg_D=D1, block_size=(m, m), update_alg=update_alg)

	h = 0.
	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=true)
	model_obc = ClassicalIsing2D(m, m, J=J, hz=1.0e-3, periodic=false)

	betas_list = collect(0.41:0.002:0.46)

	bp_results = Float64[]
	Onsager_results = Float64[]
	boundary_mps_results = Float64[]

	for beta in betas_list
		obs1 = magnetization(model_obc, i, j, update_alg, β=beta)
		obs2 = magnetization(model, i, j, exp_alg, β=beta)
		m_Onsager = Onsager_m(beta)

		push!(boundary_mps_results, obs1)
		push!(bp_results, obs2)
		push!(Onsager_results, m_Onsager)

		println("beta=$beta, BoundaryMPS=$obs1, BP=$obs2, Onsager=$(m_Onsager)")
	end


	data_path = "result/infinite_ising_critical_tn_N$(m)_h$(h)_blocksize$(block_size)_D$(D1).json"
	results = Dict("bmps"=>boundary_mps_results, "bp"=>bp_results, "Onsager"=>Onsager_results, "betas"=>betas_list)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

end
