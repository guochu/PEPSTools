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


function main(m::Int, n::Int; B::Real, block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10, msg_D=D2^2)
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
		error("path $peps_path not exist")
	end


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=msg_D, update_alg=exp_alg)
	measure_alg = bp_alg


	# gs_alg = ImaginaryTimePEPS(parameters, stepper=bp_alg, measure_alg=measure_alg, sweeps_per_measure=100, verbosity=3)

	# t = @elapsed energies, res = ground_state!(peps, h, gs_alg)
	δt = -0.01
	U = exponential(h, δt)

	sweep!(peps, U, bp_alg)

	t = @elapsed begin
		for i in 1:5
			sweep!(peps, U, bp_alg)
		end
	end

	t = t / 5

	println("it takes $t seconds per sweep.")

	return t
end

function main_vs_D()
	B = 3.1
	m = 2
	n = 2
	D2 = 3
	msg_D = D2^2

	Ds = [9, 18, 27, 36]

	ts = [main(m, n, B=B, block_size=(4, 4), D2=D2, D1=D1, msg_D=msg_D) for D1 in Ds]

	file_name = "result/runtime_vs_D.json"

	results = Dict("Ds"=>Ds, "ts"=>ts)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end


function main_vs_block_size()
	B = 3.1
	m = 2
	n = 2
	D2 = 3
	msg_D = D2^2
	D1 = 36
	blocksize = [4,8,16,32]


	ts = [main(m, n, B=B, block_size=(s, s), D2=D2, D1=D1, msg_D=msg_D) for s in blocksize]

	file_name = "result/runtime_vs_blocksize.json"

	results = Dict("blocksize"=>blocksize, "ts"=>ts)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end


function main_vs_msg_D()
	B = 3.1
	m = 2
	n = 2
	D2 = 3
	# msg_D = D2^2
	D1 = 36
	msgDs = [9, 18, 27, 36]

	ts = [main(m, n, B=B, block_size=(4, 4), D2=D2, D1=D1, msg_D=msg_D) for msg_D in msgDs]

	file_name = "result/runtime_vs_msgD.json"

	results = Dict("msgDs"=>msgDs, "ts"=>ts)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end


