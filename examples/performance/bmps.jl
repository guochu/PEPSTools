
push!(LOAD_PATH, "../../src")

using PEPSTools, Random
using Serialization, JSON

gen_peps_path(m::Int, D2::Int) = "data/xxx_su_state_m$(m)_D2$(D2)" * ".peps"

function donsweeps!(peps, U, alg, nsweeps)
	for i in 1:nsweeps
		@time sweep!(peps, U, alg)
	end
	return peps
end

function main(m::Int; D2::Int, D1::Int=2*D2^2+10) 

	Random.seed!(1234)

	peps_path = gen_peps_path(m, D2)
	peps = Serialization.deserialize(peps_path)

	h = heisenberg2D(m, m, periodic=false)
	U = exponential(h, -0.01)

	alg = BoundaryMPS(D1=D1, D2=D2)
	# gs_alg = ImaginaryTimePEPS(stepper=alg, verbosity=3)

	nsteps = 10
	nsweeps = 5

	times = Float64[]
	energies = [expectation(h, peps, alg) / prod(size(peps))]

	for i in 1:nsteps
		t = @elapsed donsweeps!(peps, U, alg, nsweeps)
		energy = expectation(h, peps, alg) / prod(size(peps))
		push!(times, t)
		push!(energies, energy)
	end

	result_path = "result/bmps_m$(m)_D2_$(D2)_D1_$(D1).json"

	results = Dict("times"=>times, "energies"=>energies)

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end
end



for D in 2:4
	main(40, D2=D)
end