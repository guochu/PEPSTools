include("../../src/includes.jl")

using Random, JSON, Statistics


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

function main_non_periodic(m::Int, h::Real=1.0e-3; block_size::Int = div(m,2))
	i = j = div(m, 2)
	

	J = 1.0

	D2 = 2
	D1 = 2*D2^2 + 5
	# D1 = 25
	update_alg = BoundaryMPS(D2=D2,D1=D1,mult_alg=IterativeCompression(D=D1))
	exp_alg = BlockBP(msg_maxiter=50, msg_tol=1.0e-10, block_size=(block_size, block_size), update_alg=update_alg)

	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=false)



	beta_min = 0.2
	beta_max = 1.5
	beta_n = 40
	del_beta = (beta_max-beta_min)/(beta_n-1)

	betas_list = [beta_min + i*del_beta for i in 0:beta_n-1]

	boundary_mps_results = Float64[]
	bp_results = []
	Onsager_results = Float64[]

	N_bp = 100
	seeds = rand(1000:10000000, N_bp) 
	for beta in betas_list
		obs1 = magnetization(model, i, j, BoundaryMPS(D2=D2, D1=D1), β=beta)
		# obs2 = [magnetization(model, i, j, exp_alg, β=beta) for l in  1:N_bp]
		obs2 = Float64[]
		for l in 1:N_bp
			Random.seed!(seeds[l])
			push!(obs2, magnetization(model, i, j, exp_alg, β=beta) )
		end
		m_Onsager = Onsager_m(beta)

		push!(boundary_mps_results, obs1)
		push!(bp_results, obs2)
		push!(Onsager_results, m_Onsager)

		println("bp results $obs2")

		println("beta=$beta, BoundaryMPS=$obs1, BP mean=$(mean(obs2)), Onsager=$(m_Onsager)")
	end


	# data_path = "result/classical_bp_full_non_periodic_tn_N$(m)_h$(h)_blocksize$(block_size)_D$(D1).json"
	# results = Dict("boundarymps"=>boundary_mps_results, "bp"=>bp_results, "Onsager"=>Onsager_results, "betas"=>betas_list)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

end

