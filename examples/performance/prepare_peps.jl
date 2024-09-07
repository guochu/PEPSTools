push!(LOAD_PATH, "../../src")

using Random
using PEPSTools
using Serialization, JSON


gen_peps_path(m::Int, D2::Int) = "data/xxx_su_state_m$(m)_D2$(D2)" * ".peps"

function simple_update_nsweeps(h, peps::PEPS, dt::Real, nsweeps::Int, alg::SimpleUpdate)
	U = exponential(h, dt)
	cpeps = CanonicalPEPS(peps)
	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end
	return PEPS(cpeps)
end

function main(m::Int;D2::Int)
	n = m
	println("run simulations for m=$m, n=$n, D2=$D2")
	L = m * n

	h = heisenberg2D(m, n, periodic=false)

	Random.seed!(3598)

	peps = randompeps(Float64, m, n, d=2, D=1)
	su_alg = SimpleUpdate(D2=D2)
	peps = simple_update_nsweeps(h, peps, -0.01, 1000, su_alg)

	peps_path = gen_peps_path(m, D2) 
	println("save peps to path $peps_path")
	Serialization.serialize(peps_path, peps)

end


for D in 2:4
	main(40, D2=D)
end