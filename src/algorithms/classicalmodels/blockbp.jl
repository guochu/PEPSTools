magnetizations(x::Classical2DModel, alg::BlockBP; β::Real) = _local_expectations(
	SquareLatticeSites(magnetization_tensors(x, β=β)), SquareTN(site_tensors(x, β=β)), alg)

_local_expectations(U::LocalCObservers, peps::SquareTN, alg::BlockBP) = _local_expectations(center_splitting(U, alg.block_size), peps, alg)
function _local_expectations(Us::Vector{BlockLocalOperator{<:AbstractArray{T, 4}}}, peps::SquareTN, alg::BlockBP) where T
	r = zeros(scalartype(peps), size(peps))
	for U in Us
		blk = peps_partition(peps, U.partition)
		r += _local_expectations(U, blk, alg)
	end
	return r
end
function _local_expectations(U::BlockLocalOperator{<:AbstractArray{T, 4}}, blk::BlockBPPartitionSquareTN, alg::BlockBP) where T
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	mult_alg = get_msg_mult_alg(alg)
	r = PeriodicArray(zeros(scalartype(blk), size(blk.peps)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r[rowindices(blk, i), colindices(blk, j)] = _local_expectations(sU, x, mult_alg)
		end
	end	
	return r.data
end


function magnetization(x::Classical2DModel, i::Int, j::Int, alg::BlockBP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	mT = magnetization_tensor(x, i, j, β=β)
	return _local_expectation(mT, i, j, peps, alg)
end
function _local_expectation(mT::AbstractArray{<:Number, 4}, i::Int, j::Int, peps::SquareTN, alg::BlockBP)
	@assert size(mT) == size(peps[i, j])
	block_size = alg.block_size
	a = div(block_size[1]-1, 2)
	b = div(block_size[2]-1, 2)
	odd_partition = lattice_partition(size(peps), block_size, (i-a, j-b), (2,2))
	blk = peps_partition(peps, odd_partition)

	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	_peps, msgl, msgr, msgu, msgd = subblock(blk, 1, 1)
	subx = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
	row_i = row(subx, a+1, mult_alg)
	return expectation_site(row_i, b+1, mT)

end

function bondenergy(x::Classical2DModel, i::Int, j::Int, alg::BlockBP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	mT1 = magnetization_tensor(x, i, j, β=β)
	mT2 = magnetization_tensor(x, i, j+1, β=β)
	block_size = alg.block_size
	a = div(block_size[1]-1, 2)
	b = div(block_size[2]-2, 2)
	odd_partition = lattice_partition(size(x), block_size, (i-a, j-b), (2,2))
	blk = peps_partition(peps, odd_partition)

	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	_peps, msgl, msgr, msgu, msgd = subblock(blk, 1, 1)
	subx = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)

	row_i = row(subx, a+1, mult_alg)
	return expectation_bond(row_i, b+1, mT1, mT2)
end


