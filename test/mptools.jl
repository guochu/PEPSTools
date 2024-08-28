println("------------------------------------")
println("|         mps operations           |")
println("------------------------------------")


@testset "mps definitions" begin
	L = 7
	for T in (Float64, ComplexF64)
		mps = randommps(T, L, d=2, D=3)
		@test bond_dimension(mps) == 3
		@test physical_dimensions(mps) == [2 for i in 1:L]
		@test scalartype(mps) == T
		@test length(mps) == 7
		@test scaling(mps) == 1
		setscaling!(mps, 2)
		@test scaling(mps) == 2
	end
end

@testset "mps canonicalize" begin
	L = 5
	for T in (Float64, ComplexF64)
		mps0 = randommps(T, L, d=3, D=5)
		mps = copy(mps0)
		n = norm(mps, iscanonical=false)
		leftorth!(mps, alg=Orthogonalize(alg=QR()))
		@test norm(mps, iscanonical=false) ≈ n
		@test isleftcanonical(mps) == true
		rightorth!(mps, alg=Orthogonalize(alg=SVD()))
		@test norm(mps, iscanonical=true) ≈ n
		@test norm(mps, iscanonical=false) ≈ n
		@test isrightcanonical(mps) == true
		@test iscanonical(mps) == true
		@test distance2(mps, mps0) ≈ 0 atol = 1.0e-8

		canonicalize!(mps, alg=Orthogonalize(alg=QR(), normalize=true))
		@test iscanonical(mps)
		@test scaling(mps) ≈ 1
	end
end


@testset "mpo mps multiplication" begin
	L = 5
	D1 = 3
	D2 = 2
	d = 2
	D = D1*D2
	tol = 1.0e-6
	for T in (Float64, ComplexF64)
		mpo = randommpo(T, L, d=d, D=D1)
		setscaling!(mpo, 1.23)
		mps = randommps(T, L, d=d, D=D2)

		r = mpo * mps
		r1, err = mult(mpo, mps, SVDCompression(trunc=truncdimcutoff(D=D, ϵ=1.0e-8)))
		@test distance2(r, r1) ≈ 0 atol = tol
		
		r2, err = mult(mpo, mps, IterativeCompression(D=D))
		@test distance2(r, r2) ≈ 0 atol = tol
	end
end
