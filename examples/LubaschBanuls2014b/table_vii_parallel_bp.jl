# push!(LOAD_PATH, "../../src")

include("../../parallel/parallel_bp.jl")

using Random
# using PEPSTools
using Serialization, JSON
# include("util.jl")

gen_peps_path(D2::Int, B::Real, block_size::Tuple{Int, Int}) = "data/Table_vii_parallel_bp_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)_B$(B).peps"
gen_result_path(D1::Int, D2::Int, B::Real, block_size::Tuple{Int, Int}) = "result/Table_vii_parallel_bp_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2)_B$(B).json"


function main(;block_size::Tuple{Int, Int}, D2::Int, B::Real, D1::Int=2*D2^2)

	(nworkers() == 9) || error("9 processes expected")

	m = 21
	n = 21
	L = m * n

	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, B=$(B), block_size=$(block_size)")

	h = ising2D(m, n, J = -1, hz = -B, periodic=false)

	Random.seed!(3598)

	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	if D2 == 2
		println("generate initial peps randomly")
		peps = randompeps(Float64, m, n, d=2, D=1)
	else
		peps_path = gen_peps_path(D2-1, B, block_size) 
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	end

	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=BoundaryMPS(D1=D1, D2=D2))
	# exp_alg = BoundaryMPS(D1=D1, D2=D2)
	exp_alg = bp_alg

	# expec_peps = expectation(h, peps, exp_alg) / prod(size(peps))
	# println("initial energy is $expec_peps")
	# energies = [real(expec_peps)]


	# append!(energies, block_iteration(peps, h, alg, 50, -0.1))

	if D2 == 2
		parameters = [(5000, -0.05, 1.0e-7), (5000, -0.005, 1.0e-9), (500, -0.001, 1.0e-9), (500, -0.0001, 1.0e-9)]

	else
		parameters = [(100, -0.005, 1.0e-9), (1000, -0.001, 1.0e-9), (200, -0.0001, 1.0e-9)]
	end

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=exp_alg, verbosity=3)

	t = @elapsed energies, res = parallel_ground_state!(peps, h, gs_alg)


	# for (maxiter, dt) in parameters
	# 	append!(energies, block_iteration(peps, h, bp_alg, maxiter, dt, tol, exp_alg))
	# end

	# final_energy = expectation(h, peps, final_exp_alg) / prod(size(peps))


	peps_path = gen_peps_path(D2, B, block_size) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)

	file_name = gen_result_path(D1, D2, B, block_size) 

	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return energies
end


Bs = [2., 2.5, 2.8, 2.9, 3., 3.1, 3.2, 3.5, 4.]

for B in Bs
	for D2 in [2,3,4]
		main(D2=D2, D1 = 2*D2^2+10, B=B, block_size=(7,7))
	end
end


