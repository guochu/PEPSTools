push!(LOAD_PATH, "../src")

using Test
using QuantumCircuits, VQC, VQC.Utilities
using QuantumSpins
using PEPSTools, Random


function convert_to_qubitsop(h::PEPSTools.SquareLatticeHamiltonian)
	m, n = size(h)
	index = LinearIndices((m, n))
	r = QubitsOperator()
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			if (j == size(h,2)) && !isnothing(h.H[i, j])
				for (a, b) in h.H[i, j]
					# push!(terms, QubitsTerm(index[i, j]=>a, index[i, j+1]=>b))
					r += QubitsTerm(index[i, j]=>a, index[i, 1]=>b)
				end	
			else
				for (a, b) in h.H[i, j]
					r += QubitsTerm(index[i, j]=>a, index[i, j+1]=>b)
				end						
			end


		end
	end
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			if i == size(h,1) && !isnothing(h.V[i, j])
				for (a, b) in h.V[i, j]
					r += QubitsTerm(index[i, j]=>a, index[1, j]=>b)
				end					
			else
				for (a, b) in h.V[i, j]
					r += QubitsTerm(index[i, j]=>a, index[i+1, j]=>b)
				end					
			end
		end
	end
	return r
end


function check_expectation_value(::Type{T}) where T
	m = 4
	n = 4
	D = 20

	Random.seed!(3568)
	peps = randompeps(T, m, n, d=2, D=2)

	h = heisenberg2D(m, n)

	h2 = convert_to_qubitsop(h)
	sv = StateVector(normalize(densevector(peps)), m*n)
	expec_exact = VQC.expectation(h2, sv) 

	# println("exact expectation value is $expec_exact")

	expec_peps = PEPSTools.expectation(h, peps, IterativeArith(D=D))

	return abs((expec_exact - expec_peps) / expec_exact) < 1.0e-6
end


function check_expectation_value_2_util(::Type{T}, D::Int) where T

	Random.seed!(1234)

	m = 4
	n = 4
	D2 = 2

	peps = randompeps(T, m+2, n+2, d=2, D=D2)
	blk = PEPSBlock(PEPS(peps[2:m+1, 2:n+1]), left=PEPSTools.random_boundary_mps(T, m, D=D2), 
		right=PEPSTools.random_boundary_mps(T, m, D=D2), up = PEPSTools.random_boundary_mps(T, n, D=D2), 
		down=PEPSTools.random_boundary_mps(T, n, D=D2))
	h = heisenberg2D(m, n)
	return PEPSTools.expectation(h, blk, IterativeArith(D=D))
end

function check_expectation_value_2(::Type{T}) where T
	r1 = check_expectation_value_2_util(T, 40)
	r2 = check_expectation_value_2_util(T, 50)
	err = abs((r1 - r2) / r2)
	# println("error is $err")
	return err < 1.0e-4
end


function check_ising_gs_util(::Type{T}, D::Int) where T

	Random.seed!(1978)

	m = 4
	n = 4
	D2 = 3

	h = ising2D(m, n, hz=3.)

	h2 = convert_to_qubitsop(h)
	E_exact, gs_exact = Utilities.ground_state(h2) 
	E_exact /= (m*n)
	# println("exact gs energy is $E_exact")

	peps = randompeps(T, m, n, d=2, D=1)
	alg = BoundaryMPS(D2=D2, D1=D)
	gs_alg = ImaginaryTimePEPS([(50, -0.05), (100, -0.01), (100, -0.002)], stepper=alg)

	energies, res = ground_state!(peps, h, gs_alg)

	return E_exact, energies
end

function check_ising_gs()
	E_exact, energies = check_ising_gs_util(Float64, 10)
	err = abs((energies[end] - E_exact) / E_exact)
	# println("error is $err")
	return err < 1.0e-5
end


@testset "test expectation value using BoundaryMPS" begin
    @test check_expectation_value(Float64)
    @test check_expectation_value(ComplexF64)
    @test check_expectation_value_2(Float64)
    @test check_expectation_value_2(ComplexF64)
end

@testset "ground state search using BoundaryMPS" begin
	@test check_ising_gs()
end


function check_bp_default_split(shape, block_size)
	m, n = shape
	peps = randompeps(Float64, m, n, d=2, D=2)

	h = heisenberg2D(m, n)
	U = exponential(h, -0.1)
	blks = default_splitting(U, block_size)

	r = sum([sum(n_nontrivial_terms.(subblocks(blk))) for blk in blks] )

	return n_nontrivial_terms(U) == sum(n_nontrivial_terms.(blks)) == r
end


function check_bp_center_split(shape, block_size)
	m, n = shape
	peps = randompeps(Float64, m, n, d=2, D=2)

	h = heisenberg2D(m, n, periodic=true)
	U = exponential(h, -0.1)
	blks = center_splitting(U, block_size)

	r = sum([sum(n_nontrivial_terms.(subblocks(blk))) for blk in blks] )

	return n_nontrivial_terms(U) == sum(n_nontrivial_terms.(blks)) == r
end
@testset "BlockBP algorithm" begin
	@test check_bp_default_split((4,4),(2,2))
	# @test check_bp_default_split((4,5),(2,3))
	@test check_bp_default_split((6,3),(2,3))

	@test check_bp_center_split((8,8),(4,4))
	@test check_bp_center_split((10,12),(4,5))
	@test check_bp_center_split((20,18),(4,4))
end



