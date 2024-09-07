push!(LOAD_PATH, "../../src")

using Random
using QuantumSpins, PEPSTools
using Serialization, JSON
# include("util.jl")

gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}, B::Real) = "data/ising_pbc_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)_B$(B)" * ".peps"
gen_result_path(m::Int, n::Int, D1::Int, D2::Int, block_size::Tuple{Int, Int}, B::Real) = "result/ising_pbc_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D1$(D1)_D2$(D2)_B$(B).json"

function simple_update_nsweeps(h, peps::PEPS, dt::Real, nsweeps::Int, alg::SimpleUpdate)
	U = exponential(h, dt)
	cpeps = CanonicalPEPS(peps)
	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end
	return PEPS(cpeps)
end

function main(m::Int, n::Int, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9), (200, -0.00001, 1.0e-9)]; B::Real, block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10)
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size), B=$B")
	L = m * n

	h = ising2D(m, n, J = -1, hz = -B, periodic=true)

	Random.seed!(3598)
	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	peps = nothing
	peps_path = gen_peps_path(m, n, D2, block_size, B) 
	if ispath(peps_path)
		println("read initial peps from path $peps_path")
		peps = Serialization.deserialize(peps_path)
	else
		if D2 == 2
			println("generate initial peps randomly")
			peps = randompeps(Float64, m, n, d=2, D=1)
		else
			tmp = gen_peps_path(m, n, D2-1, block_size, B) 
			println("read initial peps from path $tmp")
			peps = Serialization.deserialize(tmp)
		end
	end


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=div(D1, 2), update_alg=exp_alg)
	measure_alg = bp_alg

	if !ispath(peps_path)
		# do simple update first
		su_alg = SimpleUpdate(D2=D2)
		if D2 == 2
			peps = simple_update_nsweeps(h, peps, -0.01, 1000, su_alg)
		end
		# su_nsweeps = 10000
		# su_dt = -0.001
		# println("do $su_nsweeps SimpleUpdate with dt=$(su_dt)...")
		# t = @elapsed peps = simple_update_nsweeps(h, peps, su_dt, su_nsweeps, su_alg)
		# println("simple update takes $t seconds.")
		
	end

	if (!ispath(peps_path)) && (D2 == 2)
		parameters = [(10000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)]
	end

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)

	t = @elapsed energies, res = ground_state!(peps, h, gs_alg)

	# final_exp_alg = BlockBP(block_size=(m, n), msg_D=D2^2, update_alg=exp_alg)
	# final_energy = PEPSTools.center_expectation(h, peps, final_exp_alg) / prod(size(peps))

	# println("final energy using Full Block is $final_energy")

	peps_path = gen_peps_path(m, n, D2, block_size, B) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)


	file_name = gen_result_path(m, n, D1, D2, block_size, B) 

	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end


	return energies
end

# function main_load(m::Int, n::Int, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9), (200, -0.00001, 1.0e-9)]; 
# 	B::Real, peps_path::String, block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10)
# 	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size), B=$B")
# 	L = m * n

# 	h = ising2D(m, n, J = -1, hz = -B, periodic=true)

# 	Random.seed!(3598)
# 	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

# 	peps = Serialization.deserialize(peps_path)


# 	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
# 	bp_alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg=exp_alg)
# 	measure_alg = bp_alg


# 	println("parameters list $parameters")

# 	gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)

# 	t = @elapsed energies, res = ground_state!(peps, h, gs_alg)

# 	# final_exp_alg = BlockBP(block_size=(m, n), msg_D=D2^2, update_alg=exp_alg)
# 	# final_energy = PEPSTools.center_expectation(h, peps, final_exp_alg) / prod(size(peps))

# 	# println("final energy using Full Block is $final_energy")

# 	peps_path = gen_peps_path(m, n, D2, block_size, B) 
# 	println("save peps to path $peps_path")
# 	Serialization.serialize(peps_path, peps)


# 	file_name = gen_result_path(m, n, D1, D2, block_size, B) 

# 	results = Dict("parameters"=>parameters, "energies"=>energies, "runtime"=>t)

# 	open(file_name, "w") do f
# 		write(f, JSON.json(results))
# 	end


# 	return energies
# end

function main_load_single(peps_path::String, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)]; 
	block_size::Tuple{Int, Int}, B::Real, D2::Int, D1::Int=2*D2^2+10, msg_D=D2^2)
	m = 2
	n = 2	
	h = ising2D(m, n, J = -1, hz = -B, periodic=true)

	old_peps_path = gen_peps_path(m, n, D2, block_size, B) 
	if ispath(old_peps_path)
		peps = Serialization.deserialize(old_peps_path)
	else
		peps = Serialization.deserialize(peps_path)
	end

	
	println("run simulations for m=$m, n=$n, D1=$D1, D2=$D2, block_size=$(block_size), B=$B")

	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=msg_D, update_alg=exp_alg)
	measure_alg = bp_alg

	println("parameters list $parameters")

	gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)
	energies, res = ground_state!(peps, h, gs_alg)

	peps_path = gen_peps_path(m, n, D2, block_size, B) 

	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)
	file_name = gen_result_path(m, n, D1, D2, block_size, B) 

	results = Dict("parameters"=>parameters, "energies"=>energies)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end	
end

function main_load(parameters = [(1000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)]; B::Real, D2::Int, D1::Int=2*D2^2+10, msg_D=D2^2)

	m = 2
	n = 2

	block_sizes = [4,8,16,32]

	Random.seed!(3598)
	(D2 >= 2) || error("bond dimension D2 should be larger than 1.")

	# the first one
	block_size = (block_sizes[1], block_sizes[1])

	main(m, n, parameters, B=B, D2=D2, D1=D1, block_size=block_size)


	for i in 2:length(block_sizes)
		block_size = (block_sizes[i], block_sizes[i])
		peps_path = gen_peps_path(m, n, D2, (block_sizes[i-1], block_sizes[i-1]), B)
		main_load_single(peps_path, parameters, block_size=block_size, B=B, D2=D2, D1=D1, msg_D=msg_D)
	end


end


function compute_energy_obc_single(m::Int, block_size::Int, D2::Int, B::Real, k::Int)
	n = m
	D1 = 2*D2^2 + 10
	p = spin_half_matrices()
	sz = p["z"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)

	peps_path = gen_peps_path(m, n, D2, (block_size, block_size), B)
	peps4 = Serialization.deserialize(peps_path)

	L = m * k + 2

	peps = randompeps(Float64, L, L, d=2, D=D2)
	for i in 0:k-1
		for j in 0:k-1
			peps.data[m*i+2:m*(i+1)+1, m*j+2:m*(j+1)+1] = copy(peps4.data)
		end
	end

	obs = LocalObservers{Matrix{Float64}}((L, L))
	kh = div(k, 2)
	obs.data[m*kh+2:m*(kh+1)+1, m*kh+2:m*(kh+1)+1] = obs4.data	

	alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)	

	tmp = local_expectations(obs, peps, alg)

	println(tmp[m*kh+2:m*(kh+1)+1, m*kh+2:m*(kh+1)+1])

	energy = sum(tmp) / prod(m*n)

	println("sz=$energy for m=$m, B=$B, k=$k")
	return energy
end

function compute_energy_obc(D2::Int)
	m = 2
	n = m
	D1 = 2*D2^2 + 10
	# block_size = (m, n)

	# p = spin_half_matrices()
	# sz = p["z"]

	# obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	# for i in 1:length(obs4)
	# 	obs4[i] = sz
	# end
	# obs4 = LocalObservers(obs4)


	ks = collect(3:4:23)

	# Bs = vcat(collect(2.8:0.05:3), collect(3.005:0.005:3.035), collect(3.04:0.001:3.05), collect(3.055:0.005:3.1), collect(3.15:0.05:3.35) )	

	# Bs = vcat(collect(2.8:0.05:3), collect(3.01:0.01:3.14), collect(3.15:0.05:3.5) )
	B1s = vcat(collect(2.5:0.05:3.0), collect(3.01:0.01:3.04))
	B2s = collect(3.05:0.01:3.15)
	B3s = collect(3.2:0.05:3.5)

	Bs = vcat(B1s, B2s, B3s)

	mags = []
	for B in B1s
		energies = Float64[]
		for k in ks
			energy = compute_energy_obc_single(m, 4, D2, B, k)
			push!(energies, energy)
		end
		push!(mags, energies)
	end
	for B in B2s
		energies = Float64[]
		for k in ks
			energy = compute_energy_obc_single(m, 4, D2, B, k)
			push!(energies, energy)
		end
		push!(mags, energies)
	end
	for B in B3s
		energies = Float64[]
		for k in ks
			energy = compute_energy_obc_single(m, 4, D2, B, k)
			push!(energies, energy)
		end
		push!(mags, energies)
	end

	data_path = "result/ising_pbc_mag_obc_m$(m)_D1_$(D1)_D2_$(D2).json"

	println("save results to path $(data_path)")

	results = Dict("Bs"=>Bs, "ks"=>ks, "ms"=>mags)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end	
end

function compute_energy_obc_critical(D2::Int)
	m = 2
	n = m
	D1 = 2*D2^2 + 10
	# block_size = (m, n)

	# p = spin_half_matrices()
	# sz = p["z"]

	# obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	# for i in 1:length(obs4)
	# 	obs4[i] = sz
	# end
	# obs4 = LocalObservers(obs4)


	k = 25

	# Bs = vcat(collect(2.8:0.05:3), collect(3.005:0.005:3.035), collect(3.04:0.001:3.05), collect(3.055:0.005:3.1), collect(3.15:0.05:3.35) )	

	# Bs = vcat(collect(2.8:0.05:3), collect(3.01:0.01:3.14), collect(3.15:0.05:3.5) )
	# Bs = collect(3.05:0.01:3.15)

	# Bs = collect(3.09:0.01:3.15) 
	# blocksizes = [4,8,16,32]

	Bs = collect(3.04:0.01:3.11)
	blocksizes = [4, 8, 16, 32]

	mags = []
	for B in Bs
		energies = Float64[]
		for blocksize in blocksizes
			energy = compute_energy_obc_single(m, blocksize, D2, B, k)
			push!(energies, energy)
		end
		push!(mags, energies)
	end

	data_path = "result/ising_pbc_critical_mag_obc_m$(m)_D1_$(D1)_D2_$(D2).json"

	println("save results to path $(data_path)")

	results = Dict("Bs"=>Bs, "k"=>k, "blocksizes"=>blocksizes, "ms"=>mags)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end	
end

function compute_energy_pbc_single(m::Int, blocksize::Int, D2::Int, B::Real, k::Int)
	n = m
	D1 = 2*D2^2 + 10
	p = spin_half_matrices()
	sz = p["z"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)

	peps_path = gen_peps_path(m, n, D2, (blocksize, blocksize), B)
	peps = Serialization.deserialize(peps_path)

	block_size = (m*k, n*k)

	alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 

	# return center_splitting(obs4, alg.block_size, center=(m, n))

	tmp = local_expectations(center_splitting(obs4, alg.block_size, center=(m, n)), peps, alg)

	println(tmp)
	
	energy = sum(tmp) / prod(size(peps))
	println("sz=$energy for m=$m, blocksize=$(blocksize), B=$B, k=$k")

	return energy
end


function compute_energy_pbc_critical(D2::Int)
	m = 2
	n = m
	D1 = 2*D2^2 + 10
	# block_size = (m, n)

	# p = spin_half_matrices()
	# sz = p["z"]

	# obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	# for i in 1:length(obs4)
	# 	obs4[i] = sz
	# end
	# obs4 = LocalObservers(obs4)


	k = 15

	# Bs = vcat(collect(2.8:0.05:3), collect(3.005:0.005:3.035), collect(3.04:0.001:3.05), collect(3.055:0.005:3.1), collect(3.15:0.05:3.35) )	

	# Bs = vcat(collect(2.8:0.05:3), collect(3.01:0.01:3.14), collect(3.15:0.05:3.5) )
	# Bs = collect(3.05:0.01:3.15)

	# Bs = collect(3.09:0.01:3.15) 
	# blocksizes = [4,8,16,32]

	Bs = collect(3.04:0.01:3.11)
	blocksizes = [4, 8, 16]

	mags = []
	for B in Bs
		energies = Float64[]
		for blocksize in blocksizes
			energy = compute_energy_pbc_single(m, blocksize, D2, B, k)
			push!(energies, energy)
		end
		push!(mags, energies)
	end

	data_path = "result/ising_pbc_critical_mag_pbc_m$(m)_D1_$(D1)_D2_$(D2).json"

	println("save results to path $(data_path)")

	results = Dict("Bs"=>Bs, "k"=>k, "blocksizes"=>blocksizes, "ms"=>mags)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end	
end

function compute_energy_pbc(D2::Int)
	m = 2
	n = m
	D1 = 2*D2^2 + 10
	block_size = (m, n)

	p = spin_half_matrices()
	sz = p["z"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)


	ks = collect(2:3)

	# Bs = vcat(collect(2.8:0.05:3), collect(3.005:0.005:3.035), collect(3.04:0.001:3.05), collect(3.055:0.005:3.065) )	

	Bs = vcat(collect(2.8:0.05:3), collect(3.01:0.01:3.09), collect(3.1:0.05:3.5) )

	mags = []
	for B in Bs
		peps_path = gen_peps_path(m, n, D2, block_size, B)
		peps4 = Serialization.deserialize(peps_path)
		energies = Float64[]
		for k in ks
			# embed peps
			peps = repeat(peps4, k, k)
			obs = repeat(obs4, k, k)

			alg = BlockBP(block_size=size(peps), msg_D=D1, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 
			
			tmp = local_expectations(obs, peps, alg)
			# println(tmp[4*kh+2:4*(kh+1)+1, 4*kh+2:4*(kh+1)+1])
			energy = sum(tmp) / prod(size(peps))
			println("sz=$energy for B=$B, k=$k")
			push!(energies, energy)
		end
		push!(mags, energies)
	end

	data_path = "result/ising_pbc_mag_obc_m$(m)_D1_$(D1)_D2_$(D2).json"

	println("save results to path $(data_path)")

	results = Dict("Bs"=>Bs, "ks"=>ks, "ms"=>mags)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end	
end

# function compute_energy_pbc(D2::Int)
# 	m = 4
# 	n = m
# 	D1 = 2*D2^2 + 10
# 	block_size = (m, n)

# 	p = spin_half_matrices()
# 	sz = p["z"]

# 	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
# 	for i in 1:length(obs4)
# 		obs4[i] = sz
# 	end
# 	obs4 = LocalObservers(obs4)


# 	ks = collect(2:3)

# 	# Bs = vcat(collect(2.8:0.05:3), collect(3.005:0.005:3.035), collect(3.04:0.001:3.05), collect(3.055:0.005:3.065) )	

# 	Bs = vcat(collect(2.8:0.05:3), collect(3.01:0.01:3.09), collect(3.1:0.05:3.5) )

# 	mags = []
# 	for B in Bs
# 		peps_path = gen_peps_path(m, n, D2, block_size, B)
# 		peps4 = Serialization.deserialize(peps_path)
# 		energies = Float64[]
# 		for k in ks
# 			# embed peps
# 			peps = repeat(peps4, k, k)
# 			obs = repeat(obs4, k, k)

# 			alg = BlockBP(block_size=size(peps), msg_D=D1, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 
			
# 			tmp = local_expectations(obs, peps, alg)
# 			# println(tmp[4*kh+2:4*(kh+1)+1, 4*kh+2:4*(kh+1)+1])
# 			energy = sum(tmp) / prod(size(peps))
# 			println("sz=$energy for B=$B, k=$k")
# 			push!(energies, energy)
# 		end
# 		push!(mags, energies)
# 	end

# 	data_path = "result/ising_pbc_mag_obc_m$(m)_D1_$(D1)_D2_$(D2).json"

# 	println("save results to path $(data_path)")

# 	results = Dict("Bs"=>Bs, "ks"=>ks, "ms"=>mags)

# 	open(data_path, "w") do f
# 		write(f, JSON.json(results))
# 	end	
# end

function main_critical_a(m::Int, blocksize::Int, parameters = [(5000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])
	n = m
	block_size = (blocksize, blocksize)

	for B in 3.:0.01:3.07
		main(m, n, parameters, block_size=block_size, D2=2, B=B)
		# main(m, n, parameters, block_size=block_size, D2=3, B=B)
		# main(m, n, parameters, block_size=block_size, D2=4, B=B)
	end

end

function main_critical_b(m::Int, blocksize::Int, parameters = [(5000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])
	n = m
	block_size = (blocksize, blocksize)

	for B in 3.08:0.01:3.14
		main(m, n, parameters, block_size=block_size, D2=2, B=B)
		# main(m, n, parameters, block_size=block_size, D2=3, B=B)
		# main(m, n, parameters, block_size=block_size, D2=4, B=B)
	end

end

function main_critical_2(m::Int, parameters = [(10000, -0.0001, 1.0e-9)])
	n = m
	block_size = (m, n)

	for B in 3.08:0.003:3.11
		main(m, n, parameters, block_size=block_size, D2=2, B=B)
	end
end

function main_critical_3(m::Int, parameters = [(5000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])
	n = m
	block_size = (m, n)

	for B in 3.04:0.002:3.06
		main(m, n, parameters, block_size=block_size, D2=2, B=B)
		main(m, n, parameters, block_size=block_size, D2=3, B=B)
	end
end

function main_begin(m::Int, blocksize::Int, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)])
	n = m
	block_size = (blocksize, blocksize)

	for B in 2.5:0.05:2.95
		# main(m, n, parameters, block_size=block_size, D2=2, B=B)
		# main(m, n, parameters, block_size=block_size, D2=3, B=B)
		main(m, n, parameters, block_size=block_size, D2=4, B=B)
	end

end

function main_begin_begin(m::Int, blocksize::Int, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)])
	n = m
	block_size = (blocksize, blocksize)

	for B in 2.2:0.05:2.45
		main(m, n, parameters, block_size=block_size, D2=2, B=B)
		main(m, n, parameters, block_size=block_size, D2=3, B=B)
		# main(m, n, parameters, block_size=block_size, D2=4, B=B)
	end

end

function main_tail(m::Int, blocksize::Int, parameters = [(2000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)])
	n = m
	block_size = (blocksize, blocksize)

	for B in 3.15:0.05:3.5
		# main(m, n, parameters, block_size=block_size, D2=2, B=B)
		# main(m, n, parameters, block_size=block_size, D2=3, B=B)
		main(m, n, parameters, block_size=block_size, D2=4, B=B)
	end

end

# 

function main_load_critical(parameters = [(2000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])

	for B in 3.:0.01:3.14
		main_load(parameters, D2=2, B=B)
		main_load(parameters, D2=3, B=B)
	end

end

function main_load_critical_2(parameters = [(2000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])

	for B in 3.08:0.003:3.11
		main_load(parameters, D2=2, B=B)
	end
end

function main_load_critical_3(parameters = [(2000, -0.001, 1.0e-8), (2000, -0.0001, 1.0e-9)])

	for B in 3.04:0.002:3.06
		main_load(parameters, D2=2, B=B)
		main_load(parameters, D2=3, B=B)
	end
end

function main_load_begin(parameters = [(1000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)])

	for B in 2.5:0.05:2.95
		main_load(parameters, D2=2, B=B)
		main_load(parameters, D2=3, B=B)
	end

end

function main_load_tail(parameters = [(1000, -0.001, 1.0e-8), (1000, -0.0001, 1.0e-9)])

	for B in 3.15:0.05:3.5
		main_load(parameters, D2=2, B=B)
		main_load(parameters, D2=3, B=B)
	end

end
