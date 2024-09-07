const current_path = Base.@__DIR__
const PEPSTOOLS_mpath = dirname(dirname(current_path))

include(PEPSTOOLS_mpath * "/parallel/parallel_bp.jl")

using Random
using Serialization, JSON

gen_peps_path(m::Int, D2::Int) = current_path * "/data/xxx_su_state_m$(m)_D2$(D2)" * ".peps"

function donsweeps!(peps, U, alg, nsweeps)
	for i in 1:nsweeps
		@time parallel_sweep!(peps, U, alg)
	end
	return peps
end

function main(m::Int; block_size::Int, D2::Int, D1::Int=2*D2^2+10) 

	Random.seed!(1234)

	peps_path = gen_peps_path(m, D2)
	peps = Serialization.deserialize(peps_path)

	h = heisenberg2D(m, m, periodic=false)
	U = exponential(h, -0.01)

	measure_alg = BoundaryMPS(D1=D1, D2=D2)
	alg = BlockBP(block_size=(block_size, block_size), update_alg=measure_alg) 
	# gs_alg = ImaginaryTimePEPS(stepper=alg, verbosity=3)

	nsteps = 10
	nsweeps = 5

	times = Float64[]
	energies = [parallel_expectation(h, peps, measure_alg) / prod(size(peps))]

	for i in 1:nsteps
		t = @elapsed donsweeps!(peps, U, alg, nsweeps)
		energy = parallel_expectation(h, peps, measure_alg) / prod(size(peps))
		push!(times, t)
		push!(energies, energy)
	end

	np = nworkers()

	result_path = current_path * "/result/parallel_bp_np_$(np)_m$(m)_blocksize_$(block_size)_D2_$(D2)_D1_$(D1).json"

	results = Dict("times"=>times, "energies"=>energies)

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end
end

for D in 2:4
	main(40, block_size=5, D2=D)
end


