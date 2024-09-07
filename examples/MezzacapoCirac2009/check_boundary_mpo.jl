push!(LOAD_PATH, "../../src")

using Random, JSON
using PEPSTools

function main(m::Int, n::Int; block_size::Tuple{Int, Int})
	# println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size)")
	L = m * n

	h = heisenberg2D(m, n, periodic=true)

	Random.seed!(3598)

	println("generate initial peps randomly")
	peps = randompeps(Float64, m, n, d=2, D=1)

	D2 = 2
	D1 = 2*D2^2

	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=BoundaryMPS(D1=D1, D2=D2))

	bp_energies = Float64[]
	bp_energies_central = Float64[]
	bo_energies = Float64[]

	U = exponential(h, -0.05)

	println("bond dimension 2...")
	for i in 1:10
		println("the $i-th sweep..")
		for j in 1:10
			sweep!(peps, U, bp_alg)	
		end

		# bp_energy = expectation(h, peps, bp_alg) / prod(size(peps))
		# bp_energy_central = expectation(h, peps, BlockBP(bp_alg)) / prod(size(peps))
		# bo_energy = expectation(h, peps, BoundaryMPO(D1=D1, D2=D2)) / prod(size(peps))

		# push!(bp_energies, bp_energy)
		# push!(bp_energies_central, bp_energy_central)
		# push!(bo_energies, bo_energy)

		# println("bp energy=$bp_energy, bp central=$bp_energy_central, Boundary MPO energy=$bo_energy")
	end

	# data_path = "result/periodic_expectations_comparison_m$(m)_n$(n)_D2.json"
	# results = Dict("bp_energies"=>bp_energies, "bp_energies_central"=>bp_energies_central, "bo_energies"=>bo_energies)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

	println()

	D2 = 3
	D1 = 2*D2^2

	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=BoundaryMPS(D1=D1, D2=D2))

	println("bond dimension 3...")
	bp_energies = Float64[]
	bp_energies_central = Float64[]
	bo_energies = Float64[]

	U = exponential(h, -0.01)

	for i in 1:10
		println("the $i-th sweep..")
		for j in 1:10
			sweep!(peps, U, bp_alg)	
		end
		bp_energy = expectation(h, peps, bp_alg) / prod(size(peps))
		bp_energy_central = expectation(h, peps, BlockBP(bp_alg)) / prod(size(peps))
		bo_energy = expectation(h, peps, BoundaryMPO(D1=D1, D2=D2)) / prod(size(peps))

		push!(bp_energies, bp_energy)
		push!(bp_energies_central, bp_energy_central)
		push!(bo_energies, bo_energy)

		println("bp energy=$bp_energy, bp central=$bp_energy_central, Boundary MPO energy=$bo_energy")
	end


	data_path = "result/periodic_expectations_comparison_m$(m)_n$(n)_D3.json"
	results = Dict("bp_energies"=>bp_energies, "bp_energies_central"=>bp_energies_central, "bo_energies"=>bo_energies)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


end
