println("------------------------------------")
println("|              Amplitudes          |")
println("------------------------------------")


@testset "Parallel amplitudes" begin
	D = 2
	tol = 1.0e-8

	for T in (Float64, ComplexF64)
		state = randompeps(T, 3, 4, d=2, D=D)
		basis = rand((-1, 1), length(state), 4)

		amps1, back1 = Zygote.pullback(NNQS._Ψ, state, basis)
		state_back1, _ = back1(amps1)
		amps2, back2 = Zygote.pullback(Ψ_threaded, state, basis)
		state_back2, _ = back2(amps2)

		@test norm(amps1 - amps2) / norm(amps1) < tol

		for (a, b) in zip(state_back1, state_back2)
			@test norm(a-b) / norm(a) < tol
		end		

	end
end