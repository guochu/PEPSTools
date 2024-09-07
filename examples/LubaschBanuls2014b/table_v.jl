push!(LOAD_PATH, "../../src")

using Random
using PEPSTools
using Serialization, JSON
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int) = "data/Table_v_m$(m)_n$(n)_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int) = "result/Table_v_m$(m)_n$(n)_D1$(D1)_D2$(D2).json"


function main(m::Int, n::Int;D2::Int, D1::Int=2*D2^2)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2")
	L = m * n

	h = heisenberg2D(m, n, periodic=false)

	Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	if D2 == 2
		println("generate initial peps randomly")
		peps = CanonicalPEPS(randompeps(Float64, m, n, d=2, D=1))
	else
		peps_path = gen_peps_path(m, n, D2-1) 
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	end

	alg = SimpleUpdate(D2=D2)
	measure_alg = BoundaryMPS(D1=D1, D2=D2)

	# expec_peps = expectation(h, peps, alg) / prod(size(peps))
	# println("initial energy is $expec_peps")
	# energies = [real(expec_peps)]

	# tol = 1.0e-6
	if D2 == 2
		parameters = [(2000, -0.05, 1.0e-6), (2000, -0.005, 1.0e-8), (500, -0.0005, 1.0e-8), (500, -0.00005, 1.0e-8)]

	else
		parameters = [(100, -0.01, 1.0e-6), (1000, -0.001, 1.0e-8), (200, -0.0001, 1.0e-8), (200, -0.00001, 1.0e-8)]
	end

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=alg, measure_alg=measure_alg, sweeps_per_measure=50, verbosity=3)

	# for (maxiter, dt) in parameters
	# 	append!(energies, block_iteration(peps, h, alg, maxiter, dt, tol))
	# end

	t = @elapsed energies, res = ground_state!(peps, h, gs_alg)


	peps_path = gen_peps_path(m, n, D2) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)


	file_name = gen_result_path(m, n, D1, D2) 

	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return energies
end

# all the simulations in table IV

function main_10_10()
	m = 10
	n = 10
	main(m, n, D1=20, D2=2)
	main(m, n, D1=50, D2=3)
	main(m, n, D2=4, D1=2*4^2+10)
	# main(m, n, D1=100, D2=5)
	# main(m, n, D1=100, D2=6)
end

function main_14_14()
	m = 14
	n = 14
	#main(m, n, D1=20, D2=2)
	#main(m, n, D1=50, D2=3)
	#main(m, n, D1=50, D2=4)
	main(m, n, D1=100, D2=5)
	main(m, n, D1=100, D2=6)
end




