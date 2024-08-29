magnetizations(x::Classical2DModel, alg::BoundaryMPS; β::Real) = local_expectations(
	MagnetizationTensors(magnetization_tensors(x, β=β)), SquareTN(site_tensors(x, β=β)), alg)

function local_expectations(h::LocalCObservers, peps::SquareTN, alg::BoundaryMPS)
	is_nonperiodic(peps) || throw(ArgumentError("BoundaryMPS only supports OBC, use BoundaryMPO instead"))
	return local_expectations(h, SquareTNBlock(peps), alg)
end 
local_expectations(H::LocalCObservers, blk::SquareTNBlock, alg::BoundaryMPS) = local_expectations(H, blk, get_mult_alg(alg))
function local_expectations(ob::LocalCObservers, blk::SquareTNBlock, alg::MPSCompression)
	H = ob.data
	@assert size(H) == size(blk)
	m, n = size(blk)

	rH = zeros(scalartype(blk), m, n)
	
	if n_nontrivial_terms(ob) > 0
		mpsstorage = compute_H_mpsstorages(blk, alg)

		up = up_boundary(blk)
		for i in 1:m
			row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
			rH[i, :] = expectation_sites(row_i, H[i, :])
			if i != m
				mpo = mpoup(blk, i) 
				up, err = mpompsmult(mpo, up, alg)
				normalize!(up)
			end
		end		
		
	end

	return rH	
end

function row_expectations(U::AbstractVector{M}, i::Int, peps::SquareTN, alg::BoundaryMPS) where {M <: Union{AbstractArray{<:Number, 4}, Nothing}}
	is_nonperiodic(peps) || throw(ArgumentError("BoundaryMPS only supports OBC, use BoundaryMPO instead"))
	return row_expectations(U, i, SquareTNBlock(peps), get_mult_alg(alg))
end
function row_expectations(U::AbstractVector{M}, i::Int, blk::SquareTNBlock, alg::MPSCompression) where {M <: Union{AbstractArray{<:Number, 4}, Nothing}}
	row_i = row(blk, i, alg)
	return expectation_sites(row_i, U)
end

function local_expectation(mT::AbstractArray{<:Number, 4}, i::Int, j::Int, peps::SquareTN, alg::BoundaryMPS)
	is_nonperiodic(peps) || throw(ArgumentError("BoundaryMPS only supports OBC, use BoundaryMPO instead"))
	@assert size(mT) == size(peps[i, j])
	blk = SquareTNBlock(peps)
	return local_expectation(mT, i, j, blk, alg)
end

function local_expectation(mT::AbstractArray{<:Number, 4}, i::Int, j::Int, blk::SquareTNBlock, alg::BoundaryMPS)
	row_i = row(blk, i, get_mult_alg(alg))
	return expectation_site(row_i, j, mT)
end

# local magnetization
function magnetization(x::Classical2DModel, i::Int, j::Int, alg::BoundaryMPS; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	is_nonperiodic(peps) || throw(ArgumentError("BoundaryMPS only supports OBC, use BoundaryMPO instead"))
	blk = SquareTNBlock(peps)
	mT = magnetization_tensor(x, i, j, β=β)
	return local_expectation(mT, i, j, blk, alg)
end

function interactionH(x::Classical2DModel, i::Int, j::Int, alg::BoundaryMPS; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	is_nonperiodic(peps) || throw(ArgumentError("BoundaryMPS only supports OBC, use BoundaryMPO instead"))
	blk = SquareTNBlock(peps)
	mT1 = magnetization_tensor(x, i, j, β=β)
	mT2 = magnetization_tensor(x, i, j+1, β=β)

	row_i = row(blk, i, get_mult_alg(alg))
	return expectation_bond(row_i, j, mT1, mT2)
	
end
