push!(LOAD_PATH, "../../src")

using JSON
using QuantumSpins, PEPSTools


function read_entries(entries)
	A0 = zeros(2,3,3,3,3)
	for line in entries
		r = split(line, " ")
		@assert length(r) == 6
		phy, up, left, down, right, v = r
		phy = parse(Int, phy) + 1
		up = parse(Int, up) + 1
		left = parse(Int, left) + 1
		down = parse(Int, down) + 1
		right = parse(Int, right) + 1
		v = parse(Float64, v)
		A0[phy, left, up, right, down] = v
	end
	return A0
end

function read_critical_peps_site_tensor()
	peps_path = "data/BFGS100LS_D3-chi192-hx3.04438-ctm12-run0_c2ZsgdCTMe8_state.json"
	data = JSON.parsefile(peps_path)
	A0 = read_entries(data["sites"][1]["entries"])
	return A0
end


function gen_ipeps(m::Int, n::Int, A0=read_critical_peps_site_tensor())
	data = Matrix{typeof(A0)}(undef, m, n)
	for i in 1:length(data)
		data[i] = copy(A0)
	end
	return PEPS(data)
end


function measure_energy(m::Int, n::Int; block_size::Tuple{Int, Int})
	B = 3.04438
	h = ising2D(m, n, J = -1, hz = -B, periodic=true)

	D2 = 3
	D1 = 2*D2^2 + 10
	peps = gen_ipeps(m, n)
	alg = BlockBP(block_size=block_size, msg_D=D1, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 

	energy = expectation(h, peps, alg) / prod(size(peps))

	println("energy=$(energy) for block size $(block_size)")
end

function measure_sz(m::Int, n::Int; block_size::Tuple{Int, Int}, msg_D::Int)
	D2 = 3
	# D1 = 2*D2^2 + 10
	D1 = 50
	peps = gen_ipeps(m, n)

	p = spin_half_matrices()
	sz = p["z"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)

	alg = BlockBP(block_size=block_size, msg_D=msg_D, msg_tol=1.0e-6, msg_maxiter=20, update_alg = BoundaryMPS(D2=D2, D1=D1, mult_alg=IterativeCompression(D=D1, maxiter=20))) 

	tmp = local_expectations(center_splitting(obs4, alg.block_size, center=(m, n)), peps, alg)

	println(tmp)
	
	energy = sum(tmp) / prod(size(peps))
	println("sz=$energy for m=$m, blocksize=$(block_size)")

	# return energy
end

function measure_sz_obc(m::Int, k::Int)
	n = m
	D2 = 3
	D1 = 2*D2^2 + 10
	peps4 = gen_ipeps(m, n)

	p = spin_half_matrices()
	sz = p["z"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)

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

	println("sz=$energy for m=$m, k=$k")
end

function measure_sx(m::Int, n::Int; block_size::Tuple{Int, Int})
	D2 = 3
	D1 = 2*D2^2 + 10
	peps = gen_ipeps(m, n)

	p = spin_half_matrices()
	sz = p["x"]

	obs4 = Matrix{Matrix{Float64}}(undef, m, n)
	for i in 1:length(obs4)
		obs4[i] = sz
	end
	obs4 = LocalObservers(obs4)

	alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 

	tmp = local_expectations(center_splitting(obs4, alg.block_size, center=(m, n)), peps, alg)

	println(tmp)
	
	energy = sum(tmp) / prod(size(peps))
	println("sx=$energy for m=$m, blocksize=$(block_size)")

	# return energy
end
