push!(LOAD_PATH, "../../src")

using Random
using PEPSTools
using Serialization, JSON
# include("util.jl")

gen_peps_path(D2::Int, B::Real) = "data/Table_vi_D2$(D2)_B$(B).peps"
gen_result_path(D1::Int, D2::Int, B::Real) = "result/Table_vi_D1$(D1)_D2$(D2)_B$(B).json"


function main(;D2::Int, B::Real, D1::Int=2*D2^2)
	m = 11
	n = 11
	L = m * n

	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, B=$(B)")

	h = ising2D(m, n, J = -1, hz = -B, periodic=false)

	Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	if D2 == 2
		println("generate initial peps randomly")
		peps = randompeps(Float64, m, n, d=2, D=1)
	else
		peps_path = gen_peps_path(D2-1, B) 
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	end

	
	alg = BoundaryMPS(D1=D1, D2=D2)

	# expec_peps = expectation(h, peps, alg) / L
	# println("initial energy is $expec_peps")
	# energies = [real(expec_peps)]


	# append!(energies, block_iteration(peps, h, alg, 50, -0.1))

	if D2 == 2
		parameters = [(2000, -0.05, 1.0e-6), (2000, -0.005, 1.0e-8), (500, -0.0005, 1.0e-8), (500, -0.00005, 1.0e-8)]

	else
		parameters = [(1000, -0.01, 1.0e-6), (1000, -0.001, 1.0e-8), (200, -0.0001, 1.0e-8), (200, -0.00001, 1.0e-8)]
	end

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=alg, verbosity=3)

	t = @elapsed energies, res = ground_state!(peps, h, gs_alg)

	# for (maxiter, dt) in parameters
	# 	append!(energies, block_iteration(peps, h, alg, maxiter, dt, tol))
	# end


	peps_path = gen_peps_path(D2, B) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)

	file_name = gen_result_path(D1, D2, B) 

	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return energies
end

Bs = [2., 2.5, 2.8, 2.9, 3., 3.1, 3.2, 3.5, 4.]

for B in Bs
	for D2 in [2,3,4]
		main(D2=D2, D1 = 2*D2^2+10, B=B)
	end
end


