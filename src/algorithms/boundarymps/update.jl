
# full update algorithm

struct BoundaryMPS{M<:AbstractMPSArith} <: AbstractPEPSUpdateAlgorithm
	D2::Int 
	D1::Int 
	ϵ::Float64 
	als_tol::Float64 
	als_maxiter::Int 
	verbosity::Int 
	mult_alg::M
end

BoundaryMPS(; D2::Int=3, D1::Int=2*D2^2, ϵ::Real=1.0e-8, als_tol::Real=1.0e-6, als_maxiter::Int=5, verbosity::Int=1, 
	mult_alg::AbstractMPSArith = IterativeArith(maxiter=als_maxiter, tol=als_tol, D=D1, verbosity=verbosity)) = BoundaryMPS(
	D2, D1, ϵ, als_tol, als_maxiter, verbosity, mult_alg)

get_mult_alg(x::BoundaryMPS) = x.mult_alg
get_trunc(x::BoundaryMPS) = MPSTruncation(D=x.D2, ϵ=x.ϵ)

function QuantumSpins.sweep!(peps::PEPS, U::SquareLatticeOperatorBase, alg::BoundaryMPS)
	is_nonperiodic(peps) || error("BoundaryMPS does not support periodic boundary, using SimpleUpdate or BlockBP instead.")
	sweep!(PEPSBlock(peps), U, alg)
end 
function QuantumSpins.sweep!(blk::PEPSBlock, U::SquareLatticeOperatorBase, alg::BoundaryMPS)
	is_nonperiodic(U) || error("BoundaryMPS does not support periodic boundary, using SimpleUpdate or BlockBP instead.")
	@assert size(blk) == size(U)
	m, n = size(blk)

	trunc = get_trunc(alg)
	mult_alg = get_mult_alg(alg)

	# update all the horizontal terms
	als_tol, als_maxiter, als_verbosity = alg.als_tol, alg.als_maxiter, alg.verbosity
	if nontrivial_terms(U.H) > 0
		mpsstorage = compute_H_mpsstorages(blk, mult_alg)
		up = up_boundary(blk)
		for i in 1:m
			(alg.verbosity >= 3) && println("updating the $i-th row...")
			if nontrivial_terms(U.H[i, :]) > 0
				row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
				update!(row_i, U.H[i, :], trunc=trunc, maxiter=als_maxiter, tol=als_tol, verbosity=als_verbosity, normalize=true)
				blk.peps[i, :] = row_i.middle
			end
			if i != m
				mpo = mpoup(blk, i) 
				up, err = mpompsmult(mpo, up, mult_alg)
				normalize!(up, iscanonical=true)
			end
		end		
	end


	# update all the vertical terms
	if nontrivial_terms(U.V) > 0
		mpsstorage = compute_V_mpsstorages(blk, mult_alg)
		left = left_boundary(blk)
		for i in 1:n
			(alg.verbosity >= 3) && println("updating the $i-th column...")
			if nontrivial_terms(U.V[:, i]) > 0
				row_i = row_environments(mpsstorage[i+1], col_peps_as_row(blk, i), left, permute(blk.up[i], (3,2,1)), permute(blk.down[i], (3,2,1)))
				update!(row_i, U.V[:, i], trunc=trunc, maxiter=als_maxiter, tol=als_tol, verbosity=als_verbosity, normalize=true)
				blk.peps[:, i] = [permute(item, (1,5,2,3,4)) for item in row_i.middle]
			end
			if i != n
				mpo = mpoleft(blk, i)
				left, err = mpompsmult(mpo, left, mult_alg)
				normalize!(left, iscanonical=true)
			end
		end
	end
end

function compute_H_mpsstorages(blk::AbstractBlock, mult_alg::AbstractMPSArith)
	m = size(blk, 1)
	mpsstorage = Vector{Any}(undef, m+1)
	mpsstorage[m+1] = down_boundary(blk) 
	for i in m:-1:2
		mpo = mpodown(blk, i) 
		mpsstorage[i], err = mpompsmult(mpo, mpsstorage[i+1], mult_alg)
		normalize!(mpsstorage[i], iscanonical=true)
	end
	return mpsstorage
end

function compute_V_mpsstorages(blk::AbstractBlock, mult_alg::AbstractMPSArith)
	n = size(blk, 2)
	mpsstorage = Vector{Any}(undef, n+1)
	mpsstorage[n+1] = right_boundary(blk)
	for i in n:-1:2
		mpo = mporight(blk, i)
		mpsstorage[i], err = mpompsmult(mpo, mpsstorage[i+1], mult_alg)
		normalize!(mpsstorage[i], iscanonical=true)
	end	
	return mpsstorage
end


