expectation(U::SquareLatticeOperator, peps::PEPS, alg::BoundaryMPS) = expectation(U, borderedpeps(peps), get_mult_alg(alg))
expectation(U::SquareLatticeOperator, blk::BorderedPEPS, alg::BoundaryMPS) = expectation(U, blk, get_mult_alg(alg))

function expectation(U::SquareLatticeOperator, blk::BorderedPEPS, alg::MPSCompression)
	is_nonperiodic(U) || error("BoundaryMPS does not support periodic boundary, using BlockBP instead")
	@assert size(blk) == size(U)
	# m, n = size(blk)
	return SquareLatticeBonds(H=_expect_H(U.H, blk, alg), V=_expect_V(U.V, blk, alg))
end



function _expect_H(H, blk::BorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	rH = zeros(scalartype(blk), size(H, 1), size(H, 2))
	if n_nontrivial_terms(H) > 0
		mpsstorage = compute_H_mpsstorages(blk, alg)
		up = up_boundary(blk)
		for i in 1:m
			row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
			rH[i, 1:n-1] = expectation_bonds(row_i, H[i, :])
			if i != m
				mpo = mpoup(blk, i) 
				up, err = mpompsmult(mpo, up, alg)
				normalize!(up)
			end
		end		
	end

	return rH
end

function _expect_V(V, blk::BorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	rV= zeros(scalartype(blk), size(V, 1), size(V, 2))
	# update all the vertical terms
	if n_nontrivial_terms(V) > 0
		mpsstorage = compute_V_mpsstorages(blk, alg)
		left = left_boundary(blk)
		for i in 1:n
			# println("expectation of the $i-th column...")
			row_i = row_environments(mpsstorage[i+1], col_peps_as_row(blk, i), left, permute(blk.up[i], (3,2,1)), permute(blk.down[i], (3,2,1)))
			rV[1:m-1, i] = expectation_bonds(row_i, V[:, i])
			if i != n
				mpo = mpoleft(blk, i)
				left, err = mpompsmult(mpo, left, alg)
				normalize!(left)
			end
		end
	end

	return rV
end


expectation(U::LocalQObservers, peps::PEPS, alg::BoundaryMPS) = expectation(U, borderedpeps(peps), get_mult_alg(alg))
expectation(U::LocalQObservers, blk::BorderedPEPS, alg::MPSCompression) = _local_expectations(U.data, blk, alg)

function _local_expectations(U::AbstractMatrix{M}, blk::BorderedPEPS, alg::MPSCompression) where {M <: Union{AbstractMatrix, Nothing}}
	@assert size(U) == size(blk)
	m, n = size(blk)
	mpsstorage = compute_H_mpsstorages(blk, alg)
	up = up_boundary(blk)
	rH = zeros(scalartype(blk), size(blk))
	for i in 1:m
		row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
		rH[i, :] = expectation_sites(row_i, U[i, :])
		if i != m
			mpo = mpoup(blk, i) 
			up, err = mpompsmult(mpo, up, alg)
			normalize!(up)
		end
	end
	return rH
end



