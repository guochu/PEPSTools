
magnetizations(x::Classical2DModel, alg::AbstractBlockBPPEPSUpdateAlgorithm; β::Real) = local_expectations(
	MagnetizationTensors(magnetization_tensors(x, β=β)), SquareTN(site_tensors(x, β=β)), alg)

local_expectations(U::LocalClassicalObservers, peps::SquareTN, alg::BlockBP) = local_expectations(center_splitting(U, alg.block_size), peps, alg)


function local_expectations(Us::Vector{BlockLocalOperator{<:AbstractArray{T, 4}}}, peps::SquareTN, alg::AbstractBlockBPPEPSUpdateAlgorithm) where T
	r = zeros(eltype(peps), size(peps))
	for U in Us
		blk = BeliefSquareTNBlock(peps, U.partition)
		r += local_expectations(U, blk, alg)
	end
	return r
end

function local_expectations(U::BlockLocalOperator{<:AbstractArray{T, 4}}, blk::BeliefSquareTNBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm) where T
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	mult_alg = get_msg_mult_alg(alg)
	r = PeriodicArray(zeros(eltype(blk), size(blk.peps)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = PEPSBlock(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r[rowindices(blk, i), colindices(blk, j)] = local_expectations(sU, x, mult_alg)
		end
	end	
	return r.data
end


function local_expectation(mT::AbstractArray{<:Number, 4}, i::Int, j::Int, peps::SquareTN, alg::BlockBP)
	@assert size(mT) == size(peps[i, j])
	block_size = alg.block_size
	a = div(block_size[1]-1, 2)
	b = div(block_size[2]-1, 2)
	odd_partition = block_partition(size(peps), block_size, (i-a, j-b), (2,2))
	blk = peps_partition(peps, odd_partition)

	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	_peps, msgl, msgr, msgu, msgd = subblock(blk, 1, 1)
	subx = PEPSBlock(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
	row_i = row(subx, a+1, mult_alg)
	return row_magnetization(row_i, b+1, mT)

end

function magnetization(x::Classical2DModel, i::Int, j::Int, alg::BlockBP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	mT = magnetization_tensor(x, i, j, β=β)
	return local_expectation(mT, i, j, peps, alg)
end

function interactionH(x::Classical2DModel, i::Int, j::Int, alg::BlockBP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	mT1 = magnetization_tensor(x, i, j, β=β)
	mT2 = magnetization_tensor(x, i, j+1, β=β)
	block_size = alg.block_size
	a = div(block_size[1]-1, 2)
	b = div(block_size[2]-2, 2)
	odd_partition = block_partition(size(x), block_size, (i-a, j-b), (2,2))
	blk = BeliefSquareTNBlock(peps, odd_partition)

	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	_peps, msgl, msgr, msgu, msgd = subblock(blk, 1, 1)
	subx = PEPSBlock(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)

	row_i = row(subx, a+1, mult_alg)
	return row_interaction(row_i, b+1, mT1, mT2)
end


