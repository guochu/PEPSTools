push!(LOAD_PATH, "../../src")

using Random
using PEPSTools
using Serialization, JSON
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int, block_size::Tuple{Int, Int}) = "result/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2).json"

gen_final_energy_path(m::Int, n::Int, D1::Int, D2::Int, block_size::Tuple{Int, Int}) = "result/FinalEnergy_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2).json"
gen_final_energy_bpmpo_path(m::Int, n::Int, D1::Int, D2::Int, block_size::Tuple{Int, Int}) = "result/FinalEnergyBPBMPO_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2).json"

function simple_update_nsweeps(h, peps::PEPS, dt::Real, nsweeps::Int, alg::SimpleUpdate)
	U = exponential(h, dt)
	cpeps = CanonicalPEPS(peps)
	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end
	return PEPS(cpeps)
end

function main(m::Int, n::Int, parameters = [(1000, -0.01, 1.0e-8), (1000, -0.001, 1.0e-9)]; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size)")
	L = m * n

	h = heisenberg2D(m, n, periodic=true)

	Random.seed!(3598)
	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	peps = nothing
	peps_path = gen_peps_path(m, n, D2, block_size) 
	if ispath(peps_path)
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	else
		if D2 == 2
			println("generate initial peps randomly")
			peps = randompeps(Float64, m, n, d=2, D=1)
		else
			tmp = gen_peps_path(m, n, D2-1, block_size) 
			println("read initial peps from path $tmp")
			peps = Serialization.deserialize(tmp)
		end
	end


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=exp_alg)
	measure_alg = bp_alg

	if (!ispath(peps_path)) && (D2==2)
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

	if (!ispath(peps_path)) && (D2 == 2)
		parameters = [(1000, -0.01, 1.0e-8), (10000, -0.001, 1.0e-9)]
	end

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

function measure_final_energy(m::Int, n::Int; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2)
	peps_path = gen_peps_path(m, n, D2, block_size) 
	println("read peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	measure_alg = BoundaryMPO(D1=D1, D2=D2, als_maxiter=20)

	h = heisenberg2D(m, n, periodic=true)

	final_energy = expectation(h, peps, measure_alg) / prod(size(peps))

	println("final energy using BoundaryMPO is $final_energy")


	file_name = gen_final_energy_path(m, n, D1, D2, block_size) 

	results = Dict("final_energy"=>final_energy)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end	
end

function measure_final_energy_bp(m::Int, n::Int; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10)
	peps_path = gen_peps_path(m, n, D2, block_size) 
	println("read peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	measure_alg = BlockBP(block_size=(m, n), msg_D=D2^2, update_alg=BoundaryMPS(D1=D1, D2=D2, als_maxiter=20))
	 
	h = squeeze(heisenberg2D(m, n, periodic=true))

	final_energy = expectation(h, peps, measure_alg) / prod(size(peps))

	println("final energy using BlockBP is $final_energy")

	file_name = "result/FinalEnergy_bp_m$(m)_n$(n)_D1$(D1)_D2$(D2).json"

	results = Dict("final_energy"=>final_energy)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end	

	return final_energy
end

function check_final_energy(m::Int, n::Int; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10)
	peps_path = gen_peps_path(m, n, D2, block_size) 
	println("read peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	measure_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=BoundaryMPS(D1=D1, D2=D2))
	 
	h = squeeze(heisenberg2D(m, n, periodic=true))

	Us = center_splitting(h, measure_alg.block_size)

	final_energies = [expectationfull(U, blockbp_environments(peps, U.partition), measure_alg) for U in Us]

	# println("final energy using BoundaryMPO is $final_energy")

	file_name = "result/FinalEnergies_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2).json"

	results = Dict("final_energies"=>final_energies)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end	

	return final_energies
end

# function measure_final_energy_bpmpo(m::Int, n::Int; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2)
# 	peps_path = gen_peps_path(m, n, D2, block_size) 
# 	println("read peps from path $peps_path")
# 	peps = Serialization.deserialize(peps_path)

# 	measure_alg = BPBoundaryMPO(D1=D1, D2=D2, mult_alg=BPIterativeCompress(msg_maxiter=200, msg_tol=1.0e-10, D=D1))

# 	h = heisenberg2D(m, n, periodic=true)

# 	final_energy = expectation(h, peps, measure_alg) / prod(size(peps))

# 	println("final energy using BPBoundaryMPO is $final_energy")


# 	file_name = gen_final_energy_bpmpo_path(m, n, D1, D2, block_size) 

# 	results = Dict("final_energy"=>final_energy)

# 	open(file_name, "w") do f
# 		write(f, JSON.json(results))
# 	end	
# end

function main_10_10()
	m = 10
	n = 10
	block_size = (5, 5)
	main(m, n, block_size=block_size, D2=2)
	main(m, n, block_size=block_size, D2=3)
	# main(m, n, block_size=block_size, D2=4)
	# main(m, n, block_size=block_size, D2=5)
	# main(m, n, block_size=block_size, D2=6)
end

function main_8_8()
	m = 8
	n = 8
	block_size = (m, n)
	main(m, n, block_size=block_size, D2=2)
	main(m, n, block_size=block_size, D2=3)
	main(m, n, block_size=block_size, D2=4)
	# main(m, n, block_size=block_size, D2=5)
	# main(m, n, block_size=block_size, D2=6)
end

function main_5_5()
	m = 5
	n = 5
	block_size = (m, n)
	main(m, n, block_size=block_size, D2=2)
	main(m, n, block_size=block_size, D2=3)
	main(m, n, block_size=block_size, D2=4)
	# main(m, n, block_size=block_size, D2=5)
	# main(m, n, block_size=block_size, D2=6)
end


function main_6_6()
	m = 6
	n = 6
	block_size = (6, 6)
	main(m, n, block_size=block_size, D2=2)
	main(m, n, block_size=block_size, D2=3)
	main(m, n, block_size=block_size, D2=4)
	# main(m, n, block_size=block_size, D2=5)
	# main(m, n, block_size=block_size, D2=6)
end


function main_16_16()
	m = 16
	n = 16
	block_size = (8, 8)
	# main(m, n, block_size=block_size, D2=2)
	# main(m, n, block_size=block_size, D2=3)
	main(m, n, block_size=block_size, D2=4)
	# main(m, n, block_size=block_size, D2=5)
	# main(m, n, block_size=block_size, D2=6)
end

