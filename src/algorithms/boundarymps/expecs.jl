

expectation(U::SquareLatticeOperator, peps::PEPS, alg::BoundaryMPS) = expectation(U, borderedpeps(peps), get_mult_alg(alg))
expectation(U::SquareLatticeOperator, blk::BorderedPEPS, alg::BoundaryMPS) = expectation(U, blk, get_mult_alg(alg))
expectation(h::SquareLatticeHamiltonian, blk::BorderedPEPS, alg::MPSCompression) = expectation(squeeze(h), blk, alg)
expectation(h::SquareLatticeHamiltonian, blk::BorderedPEPS, alg::BoundaryMPS) = expectation(h, blk, get_mult_alg(alg))
expectation(h::SquareLatticeHamiltonian, peps::PEPS, alg::MPSCompression) = expectation(h, borderedpeps(peps), alg)
expectation(h::SquareLatticeHamiltonian, peps::PEPS, alg::BoundaryMPS) = expectation(h, peps, get_mult_alg(alg))


function expectation(U::SquareLatticeOperator, blk::BorderedPEPS, alg::MPSCompression)
	is_nonperiodic(U) || error("BoundaryMPS does not support periodic boundary, using BlockBP instead.")
	@assert size(blk) == size(U)
	# m, n = size(blk)
	return sum(_expect_H(U.H, blk, alg)) + sum(_expect_V(U.V, blk, alg))
end

function expectationfull(U::SquareLatticeOperator, blk::BorderedPEPS, alg::MPSCompression)
	is_nonperiodic(U) || error("BoundaryMPS does not support periodic boundary, using BlockBP instead.")
	@assert size(blk) == size(U)
	# m, n = size(blk)
	return Dict("H"=>_expect_H(U.H, blk, alg), "V"=>_expect_V(U.V, blk, alg))
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


local_expectations(U::LocalQObservers, peps::PEPS, alg::BoundaryMPS) = local_expectations(U, borderedpeps(peps), get_mult_alg(alg))
local_expectations(U::LocalQObservers, blk::BorderedPEPS, alg::MPSCompression) = local_expectations(U.data, blk, alg)
local_expectations(U::AbstractMatrix{M}, peps::PEPS, alg::BoundaryMPS) where {M <: Union{AbstractMatrix, Nothing}} = local_expectations(
	U, borderedpeps(peps), get_mult_alg(alg))

function local_expectations(U::AbstractMatrix{M}, blk::BorderedPEPS, alg::MPSCompression) where {M <: Union{AbstractMatrix, Nothing}}
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

function row_expectations(U::AbstractVector{M}, i::Int, peps::PEPS, alg::BoundaryMPS) where {M <: Union{AbstractMatrix, Nothing}}
	return row_expectations(U, i, borderedpeps(peps), get_mult_alg(alg))
end
function row_expectations(U::AbstractVector{M}, i::Int, blk::BorderedPEPS, alg::MPSCompression) where {M <: Union{AbstractMatrix, Nothing}}
	row_i = row(blk, i, alg)
	return local_expectations(row_i, U)
end



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
rdm2s(blk::BorderedPEPS, alg::MPSCompression) = Dict("H"=>rdm2sH(blk, alg), "V"=>rdm2sV(blk, alg))


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

