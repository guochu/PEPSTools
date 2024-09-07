push!(LOAD_PATH, "../../src")

using Random
using PEPSTools
using Serialization, JSON


gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"

gen_peps_path(m::Int, n::Int, D2::Int, block_size::Tuple{Int, Int}) = "data/Table2_bp_m$(m)_n$(n)_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)" * ".peps"


function energy_finite_single(m::Int, blocksize::Int, D2::Int, k::Int, p::Int)
	n = m
	D1 = 2*D2^2 + p

	block_size = (blocksize, blocksize)

	peps_path = gen_peps_path(m, n, D2, block_size)
	peps4 = Serialization.deserialize(peps_path)

	h4 = heisenberg2D(m, n, periodic=true)

	Random.seed!(1234)
	# enlarge peps
	# k = 2
	L = m * k + 2
	peps = randompeps(Float64, L, L, d=2, D=D2)
	for i in 0:k-1
		for j in 0:k-1
			peps.data[m*i+2:m*(i+1)+1, m*j+2:m*(j+1)+1] = copy(peps4.data)
		end
	end

	# enlarge h
	h = SquareLatticeHamiltonian(Float64, L, L)
	kh = div(k, 2)
	h.H[m*kh+2:m*(kh+1)+1, m*kh+2:m*(kh+1)+1] = h4.H
	h.V[m*kh+2:m*(kh+1)+1, m*kh+2:m*(kh+1)+1] = h4.V

	alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)
	energy = expectation(h, peps, alg) / prod(m*n)

	return energy
end

function main_finite(m::Int, blocksize::Int, p::Int=10)
	Ds = [2,3,4]
	ks = collect(3:2:25)

	energies = []
	for D2 in Ds
		energies_tmp = Float64[]
		for k in ks
			energy = energy_finite_single(m, blocksize, D2, k, p)
			println("D2=$(D2), k=$(k), energy = $(energy)")
			push!(energies_tmp, energy)
		end
		push!(energies, energies_tmp)
	end

	file_name = "result/Table2_final_energy_m_$(m)_blocksize_$(blocksize)_p_$(p).json"

	results = Dict("Ds"=>Ds, "ks"=>ks, "energies"=>energies)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end


function energy_infinite_single(m::Int, blocksize::Int, D2::Int, k::Int, p::Int)
	n = m
	D1 = 2*D2^2 + p

	block_size = (blocksize, blocksize)

	peps_path = gen_peps_path(m, n, D2, block_size)
	peps = Serialization.deserialize(peps_path)

	h = heisenberg2D(m, n, periodic=true)

	block_size = (blocksize*k, blocksize*k)

	alg = BlockBP(block_size=block_size, msg_D=D2^2, update_alg = BoundaryMPS(D2=D2, D1=D1, als_maxiter=20)) 

	energy = expectation(h, peps, alg) / prod(size(peps))


	return energy
end

function main_infinite(m::Int, blocksize::Int, p::Int=10)
	Ds = [2,3,4]
	ks = collect(3:2:9)

	energies = []
	for D2 in Ds
		energies_tmp = Float64[]
		for k in ks
			energy = energy_infinite_single(m, blocksize, D2, k, p)
			push!(energies_tmp, energy)
			println("D2=$(D2), k=$(k), energy = $(energy)")
		end
		push!(energies, energies_tmp)
	end

	file_name = "result/Table2_bp_final_energy_m_$(m)_blocksize_$(blocksize)_p_$(p).json"

	results = Dict("Ds"=>Ds, "ks"=>ks, "energies"=>energies)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end
