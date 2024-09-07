push!(LOAD_PATH, "../../src")

using Random
using PEPSTools, PeriodicMPS
using Serialization, JSON

gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int, block_size::Tuple{Int, Int}) = "result/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2).json"



function main(parameters = [(1000, -0.01, 1.0e-8), (1000, -0.001, 1.0e-9)];D2::Int, D1::Int=2*D2^2+10)
	m = 2
	n = 2
	block_size = (8, 8)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size)")

	h = heisenberg2D(m, n, periodic=true)

	Random.seed!(1356)
	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	peps = nothing
	old_peps_path = gen_peps_path(m, n, D2, block_size) 
	if ispath(old_peps_path)
		println("read initial peps from path $old_peps_path")
		peps = Serialization.deserialize(old_peps_path)
	else
		peps4_path = gen_peps_path(m, n, D2, (4,4)) 
		println("read initial 4 by 4 peps from path $peps4_path")
		peps = Serialization.deserialize(peps4_path)
	end


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=exp_alg)
	measure_alg = bp_alg


	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)

	t = @elapsed energies, res = ground_state!(peps, h, gs_alg)

	# final_exp_alg = BlockBP(block_size=(m, n), msg_D=D2^2, update_alg=exp_alg)
	# final_energy = PEPSTools.center_expectation(h, peps, final_exp_alg) / prod(size(peps))

	# println("final energy using Full Block is $final_energy")

	peps_path = gen_peps_path(m, n, D2, block_size) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)


	file_name = gen_result_path(m, n, D1, D2, block_size) 

	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end


	return energies
end
