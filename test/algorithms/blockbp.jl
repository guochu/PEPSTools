

function check_bp_default_split(shape, block_size)
	m, n = shape
	peps = randompeps(Float64, m, n, d=2, D=2)

	h = heisenberg2D(m, n)
	U = exponential(h, -0.1)
	blks = default_splitting(U, block_size)

	r = sum([sum(nontrivial_terms.(subblocks(blk))) for blk in blks] )

	return nontrivial_terms(U) == sum(nontrivial_terms.(blks)) == r
end


function check_bp_center_split(shape, block_size)
	m, n = shape
	peps = randompeps(Float64, m, n, d=2, D=2)

	h = heisenberg2D(m, n, periodic=true)
	U = exponential(h, -0.1)
	blks = center_splitting(U, block_size)

	r = sum([sum(nontrivial_terms.(subblocks(blk))) for blk in blks] )

	return nontrivial_terms(U) == sum(nontrivial_terms.(blks)) == r
end

@testset "BlockBP splitting" begin
	@test check_bp_default_split((4,4),(2,2))
	# @test check_bp_default_split((4,5),(2,3))
	@test check_bp_default_split((6,3),(2,3))

	@test check_bp_center_split((8,8),(4,4))
	@test check_bp_center_split((10,12),(4,5))
	@test check_bp_center_split((20,18),(4,4))
end