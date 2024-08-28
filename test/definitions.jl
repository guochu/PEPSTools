println("------------------------------------")
println("|      basic peps definitions      |")
println("------------------------------------")


@testset "PEPS definitions" begin
	peps = randompeps(randn, Float64, 5, 6, d=2, D=3, periodic=false)
	@test scalartype(peps) == Float64
	@test size(peps, 1) == 5
	@test size(peps, 2) == 6
	@test size(peps) == (5, 6)
	@test is_periodic(peps) == false

	peps = randompeps(randn, ComplexF64, 5, 6, d=2, D=2, periodic=true)
	@test scalartype(peps) == ComplexF64
	@test is_periodic(peps) == true
end