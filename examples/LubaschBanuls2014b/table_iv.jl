# push!(LOAD_PATH, "../../src")
include("../../src/includes.jl")

using Random
using Serialization, JSON
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int) = "data/Table_iv_m$(m)_n$(n)_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int) = "result/Table_iv_m$(m)_n$(n)_D1$(D1)_D2$(D2).json"

function simple_update_nsweeps(h, peps::PEPS, dt::Real, nsweeps::Int, alg::SimpleUpdate)
	U = exponential(h, dt)
	cpeps = CanonicalPEPS(peps)
	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end
	return PEPS(cpeps)
end

function main(m::Int, n::Int, parameters = [(10000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9), (1000, -0.00001, 1.0e-9)];D2::Int, D1::Int=2*D2^2+10)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2")
	L = m * n

	h = heisenberg2D(m, n, periodic=false)

	Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")


	peps = nothing
	peps_path = gen_peps_path(m, n, D2) 
	if ispath(peps_path)
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	else
		if D2 == 2
			println("generate initial peps randomly")
			peps = randompeps(Float64, m, n, d=2, D=1)
		else
			tmp = gen_peps_path(m, n, D2-1) 
			println("read initial peps from path $tmp")
			peps = Serialization.deserialize(tmp)
		end
	end

	if !ispath(peps_path)
		# do simple update first
		su_alg = SimpleUpdate(D2=D2)
		if D2 == 2
			peps = simple_update_nsweeps(h, peps, -0.01, 1000, su_alg)
		end
		su_nsweeps = 10000
		su_dt = -0.001
		println("do $su_nsweeps SimpleUpdate with dt=$(su_dt)...")
		t = @elapsed peps = simple_update_nsweeps(h, peps, su_dt, su_nsweeps, su_alg)
		println("simple update takes $t seconds.")
		# expec_peps = expectation(h, peps, exp_alg) / prod(size(peps))
		# println("initial energy is $expec_peps")
		# energies = [real(expec_peps)]
		
	end

	alg = BoundaryMPS(D1=D1, D2=D2)

	# expec_peps = expectation(h, peps, alg) / prod(size(peps))
	# println("initial energy is $expec_peps")
	# energies = [real(expec_peps)]

	# tol = 1.0e-6
	if (!ispath(peps_path)) && (D2 == 2)
		parameters = [(20000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9), (500, -0.00001, 1.0e-9)]
	end

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=alg, sweeps_per_measure=50, verbosity=3)

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
	main(m, n, D2=2)
	main(m, n, D2=3)
	main(m, n, D2=4)
	# main(m, n, D1=100, D2=5)
	# main(m, n, D1=100, D2=6)
end

function main_14_14()
	m = 14
	n = 14
	# main(m, n, D2=2)
	main(m, n, D2=3)
	# main(m, n, D2=4)
end




