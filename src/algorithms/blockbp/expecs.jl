function energy(U::SquareLatticeOperator, peps::PEPS, alg::BlockBP) 
	r = expectation(U, peps, alg)
	return sum(r.H) + sum(r.V)
end
expectation(U::SquareLatticeOperator, peps::PEPS, alg::BlockBP) = expectation(center_splitting(U, alg.block_size), peps, alg)
function expectation(Us::Vector{<:BlockOperator}, peps::PEPS, alg::BlockBP)
	rH = PeriodicArray(zeros(scalartype(peps), size(peps)))
	rV = PeriodicArray(zeros(scalartype(peps), size(peps)))
	r = SquareLatticeBonds(H=rH, V=rV)

	for U in Us
		blk = peps_partition(peps, U.partition)
		expectation!(r, U, blk, alg)
	end
	return r
end

function expectation(U::BlockOperator, blk::BlockBPPartitionPEPS, alg::BlockBP)
	rH = PeriodicArray(zeros(scalartype(blk), size(blk)))
	rV = PeriodicArray(zeros(scalartype(blk), size(blk)))
	r = SquareLatticeBonds(H=rH, V=rV)
	return expectation!(r, U, blk, alg)
end

function expectation!(r::SquareLatticeBonds, U::BlockOperator, blk::BlockBPPartitionPEPS, alg::BlockBP)
	@assert blk.partition == U.partition
	@assert size(r) == size(blk)
	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)

	# rH = PeriodicArray(zeros(scalartype(blk), size(blk)))
	# rV = PeriodicArray(zeros(scalartype(blk), size(blk)))
	rH, rV = r.H, r.V

	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r = expectation(sU, x, mult_alg)
			rH[rowindices(blk, i), colindices(blk, j)] .+= r.H
			rV[rowindices(blk, i), colindices(blk, j)] .+= r.V
		end
	end	
	# return SquareLatticeBonds(H=rH.data, V=rV.data)
	return r
end

expectation(U::LocalQObservers, peps::PEPS, alg::BlockBP) = expectation(center_splitting(U, alg.block_size), peps, alg)


function expectation(Us::Vector{BlockLocalOperator{M}}, peps::PEPS, alg::BlockBP) where {M<:AbstractMatrix}
	r = zeros(scalartype(peps), size(peps))
	for U in Us
		blk = peps_partition(peps, U.partition)
		r .+= expectation(U, blk, alg)
	end
	return r
end

function expectation(U::BlockLocalOperator{M}, blk::BlockBPPartitionPEPS, alg::BlockBP) where {M<:AbstractMatrix}
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	mult_alg = get_msg_mult_alg(alg)
	r = PeriodicArray(zeros(scalartype(blk), size(blk.peps)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)

			tmp = expectation(sU, x, mult_alg)
			row_pos, row_idx = rowcenterenumerate(blk.partition, i)
			col_pos, col_idx = colcenterenumerate(blk.partition, j)

			r[row_idx, col_idx] = tmp[row_pos, col_pos]
		end
	end	
	return r.data
end

