println("------------------------------------")
println("|         PEPS expectation         |")
println("------------------------------------")


@testset "Expectation value and rdms" begin

_distance(x::Nothing, y::Nothing) = 0.
_distance(x::AbstractArray, y::AbstractArray) = norm(x - y)

	sz = [1. 0.; 0 -1]
	m = 8

	obs_data = Matrix{typeof(sz)}(undef, (m, m))
	for i in 1:m, j in 1:m
		obs_data[i, j] = sz
	end
	obs = SquareLatticeSites(obs_data)

	D2 = 2
	D1 = 20
	alg1 = BoundaryMPS(D2=D2,D1=D1,als_maxiter=50,mult_alg=IterativeCompression(D=D1,maxiter=10))
	m2 = div(m, 2)
	alg2 = BlockBP(msg_maxiter=10, msg_tol=1.0e-11, msg_D=10, block_size=(m2, m2), update_alg=alg1)

	tol = 1.0e-2
	for T in (Float64, ComplexF64)
		peps = randompeps(T, m, m, d=2, D=2, periodic=false)

		r1 = expectation(obs, peps, alg1)
		r2 = expectation(obs, peps, alg2)

		@test norm(r1 - r2) / m < tol

		o1 = rdm1s(peps, alg1)
		o2 = rdm1s(peps, alg2) 

		diss = [norm(x-y) for (x, y) in zip(o1, o2)]
		@test sum(diss) / length(diss) < tol


		p1 = rdm2s(peps, alg1)
		p2 = rdm2s(peps, alg2) 

		diss = [_distance(x, y) for (x, y) in zip(p1.H, p2.H)]
		@test sum(diss) / length(diss) < tol

		diss = [_distance(x, y) for (x, y) in zip(p1.V, p2.V)]
		@test sum(diss) / length(diss) < tol
	end

end