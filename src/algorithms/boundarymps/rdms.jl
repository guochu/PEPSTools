# reduced density matrices
rdm1s(peps::PEPS, alg::BoundaryMPS) = rdm1s(peps, get_mult_alg(alg))
rdm1s(peps::PEPS, alg::MPSCompression) = rdm1s(borderedpeps(peps), alg)

"""
	rdm1s(blk::BorderedPEPS, alg::MPSCompression)

Return all the local reduced density matrices.
"""
function rdm1s(blk::BorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	rH = Matrix{Matrix{scalartype(blk)}}(undef, size(blk))

	mpsstorage = compute_H_mpsstorages(blk, alg)
	up = up_boundary(blk)
	for i in 1:m
		row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
		rH[i, :] = rdm1s(row_i)
		if i != m
			mpo = mpoup(blk, i) 
			up, err = mpompsmult(mpo, up, alg)
			normalize!(up)
		end
	end		

	return rH
end

rdm2s(peps::PEPS, alg::BoundaryMPS) = rdm2s(peps, get_mult_alg(alg))
rdm2s(peps::PEPS, alg::MPSCompression) = rdm2s(borderedpeps(peps), alg)
"""
	rdm2s(blk::BorderedPEPS, alg::MPSCompression)

Return all the nearest neighbour twobody reduced density matrices.
"""
rdm2s(blk::BorderedPEPS, alg::MPSCompression) = SquareLatticeBonds(H=rdm2sH(blk, alg), V=rdm2sV(blk, alg))


rdm2sH(peps::PEPS, alg::BoundaryMPS) = rdm2sH(peps, get_mult_alg(alg))
rdm2sH(peps::PEPS, alg::MPSCompression) = rdm2sH(borderedpeps(peps), alg)
function rdm2sH(blk::BorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	rH = Matrix{Array{scalartype(blk), 4}}(undef, size(blk,1), size(blk,2)-1)

	mpsstorage = compute_H_mpsstorages(blk, alg)
	up = up_boundary(blk)
	for i in 1:m
		row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
		rH[i, 1:n-1] = rdm2s(row_i)
		if i != m
			mpo = mpoup(blk, i) 
			up, err = mpompsmult(mpo, up, alg)
			normalize!(up)
		end
	end		

	return rH
end

rdm2sV(peps::PEPS, alg::BoundaryMPS) = rdm2sV(peps, get_mult_alg(alg))
rdm2sV(peps::PEPS, alg::MPSCompression) = rdm2sV(borderedpeps(peps), alg)
function rdm2sV(blk::BorderedPEPS, alg::MPSCompression)
	m, n = size(blk)

	rV= Matrix{Array{scalartype(blk), 4}}(undef, size(blk, 1)-1, size(blk, 2))
	# update all the vertical terms
	mpsstorage = compute_V_mpsstorages(blk, alg)
	left = left_boundary(blk)
	for i in 1:n
		# println("expectation of the $i-th column...")
		row_i = row_environments(mpsstorage[i+1], col_peps_as_row(blk, i), left, permute(blk.up[i], (3,2,1)), permute(blk.down[i], (3,2,1)))
		rV[1:m-1, i] = rdm2s(row_i)
		if i != n
			mpo = mpoleft(blk, i)
			left, err = mpompsmult(mpo, left, alg)
			normalize!(left)
		end
	end

	return rV
end

