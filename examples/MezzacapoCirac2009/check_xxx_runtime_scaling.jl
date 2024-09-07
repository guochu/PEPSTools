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

function main(m::Int, n::Int; block_size::Tuple{Int, Int}, D2::Int, D1::Int=2*D2^2+10, msg_D=D2^2)
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
		error("path $peps_path not exist")
	end


	exp_alg = BoundaryMPS(D1=D1, D2=D2, als_maxiter=10)
	bp_alg = BlockBP(block_size=block_size, msg_D=msg_D, update_alg=exp_alg)
	measure_alg = bp_alg


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
	m = 2
	n = 2

	Ds = [2,3,4,5,6]

	ts = [main(m, n, block_size=(4, 4), D2=D2) for D2 in Ds]

	file_name = "result/xxx_runtime_vs_D2.json"

	results = Dict("Ds"=>Ds, "ts"=>ts)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end
end

