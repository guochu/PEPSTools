

# compute the new two site tensors in the center
function center_minimization(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, 
	X::AbstractArray{<:Number}, U::AbstractArray{<:Number, 4}; normalize::Bool=true, trunc::TruncationScheme, maxiter::Int, tol::Real, verbosity::Int)
	# regularization
	aL, aR, X, invL, invR = _regularize_aL_aR_X(aL, aR, X)

	# compute initial loss
	@tensor down[2,4,6] := aL[1,2,3] * X[1,4,5] * aR[3,6,5] 
	@tensor bx = conj(down[1,2,3]) * U[1,3,4,5] * down[4,2,5]
	init_loss = real(dot(down, down) - 2 * real(bx))
	# println("initial loss is $init_loss")
	

	# compute initial guess for aLn and aRn
	@tensor m[1,6,7,5] := aL[1,2,3] * aR[3,4,5] * U[6,7,2,4]
	s1, s2, s3, s4 = size(m)
	u, s, v, err = tsvd!(reshape(m, (s1*s2, s3*s4)); trunc=trunc)
	# normalize!(s)
	Dn = length(s)
	sqrt_s = Diagonal(sqrt.(s))
	u = u * sqrt_s
	v = sqrt_s * v
	aLn = reshape(u, (s1, s2, Dn))
	aRn = reshape(v, (Dn, s3, s4))

	# compute loss after SVD
	@tensor down2[2,4,6] := aLn[1,2,3] * X[1,4,5] * aRn[3,6,5] 
	@tensor bx2 = conj(down2[1,2,3]) * U[1,3,4,5] * down[4,2,5]
	loss_svd = real(dot(down2, down2) - 2 * real(bx2))

	# println("loss after svd is $loss_svd")

	iter = 1
	errs = [init_loss, loss_svd]
	# how to compute the loss here?
	while iter <= maxiter
		Nl = compute_Nl(aRn, X)
		bl = compute_bl(aL, aR, X, U, aRn)
		aLn, d1 = local_solve(Nl, bl)

		# push!(errs, d1)

		Nr = compute_Nr(aLn, X)
		br = compute_br(aL, aR, X, U, aLn)
		aRn, d2 = local_solve(Nr, br)

		err = abs((d2 - errs[end]) / init_loss)

		push!(errs, d2)

		# println("relative error in $iter-th iteration is $err")

		if err <= tol
			(verbosity > 2) && println("early converge in $iter with relative error $err")
			break
		end

		iter += 1
	end

	if (verbosity > 1) && (iter > maxiter)
		println("fail to converge to precision $tol in $maxiter iterations.")
	end

	# println("distances: $errs")

	@tensor aLt[2,3,4] := invL[1,2] * aLn[1,3,4]
	@tensor aRt[3,4,1] := invR[1,2] * aRn[3,4,2]

	return _regularize_aL_aR(aLt, aRt, X, normalize=normalize)
end

# compute the left and right norm matrices
function compute_Nl(aR::AbstractArray{<:Number, 3}, X::AbstractArray{<:Number, 3})
	@tensor tmp[1,2,4,5] := X[1,2,3] * aR[4,5,3]
	@tensor Nl[1,3,5,6] := conj(tmp[1,2,3,4]) * tmp[5,2,6,4]
	return Nl
end
function compute_Nr(aL::AbstractArray{<:Number, 3}, X::AbstractArray{<:Number, 3})
	@tensor tmp[2,3,4,5] := aL[1,2,3] * X[1,4,5] 
	@tensor Nr[2,4,5,6] := conj(tmp[1,2,3,4]) * tmp[1,5,3,6]
	return Nr
end

function compute_bl(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, 
	X::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}, aRn::AbstractArray{<:Number, 3})
	@tensor down[7,4,8] := aL[1,2,3] * X[1,4,5] * aR[3,6,5] * U[7,8,2,6]
	@tensor up[1,2,4,5] := X[1,2,3] * aRn[4,5,3]
	@tensor bl[1,3,5] := conj(up[1,2,3,4]) * down[5,2,4]
	return bl
end

function compute_br(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, 
	X::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}, aLn::AbstractArray{<:Number, 3})
	@tensor down[7,4,8] := aL[1,2,3] * X[1,4,5] * aR[3,6,5] * U[7,8,2,6]
	@tensor up[2,3,4,5] := aLn[1,2,3] * X[1,4,5]
	@tensor br[2,4,5] := conj(up[1,2,3,4]) * down[1,3,5]
	return br
end

function _compute_distance(N::AbstractMatrix, b::AbstractMatrix, x::AbstractMatrix)
	return real(-dot(b, x))
end

function local_solve(N::AbstractArray{<:Number, 4}, b::AbstractArray{<:Number, 3})
	s1, s2, s3 = size(b)
	L = s1 * s2
	N2 = reshape(N, (L, L))
	b2 = reshape(b, (L, s3))
	# println("here *********************")
	# println(size(N2))
	# println(size(b2))
	# println("conditional number is $(cond(N2))")
	r2 = N2 \ b2
	return permute(reshape(r2, size(b)), (1,3,2) ), _compute_distance(N2, b2, r2)
end


function make_positive_util(N2::AbstractArray{<:Number, 2})
	N2H = (N2 + N2') / 2
	if real(tr(N2H)) < 0.
		N2H = -N2H
	end
	eigvalues, eigvectors = eigen(Hermitian(N2H))
	# println("eigenvalues $eigvalues")

	eigvalue_threshold = 1.0e-12
	pos = findfirst(x->x >= eigvalue_threshold, eigvalues / norm(eigvalues))
	isnothing(pos) && error("environment becomes trivial")
	# pos = 1

	eigvalues = eigvalues[pos:end]
	# println("eigenvalues $eigvalues")
	eigvectors = eigvectors[:, pos:end]
	X = eigvectors * Diagonal(sqrt.(eigvalues))
	return X
end

function make_N_positive(N::AbstractArray{<:Number, 2})
	Dl, Dr = size(N)
	sqrt_Dl = round(Int, sqrt(Dl))
	sqrt_Dr = round(Int, sqrt(Dr))
	@assert (sqrt_Dl^2 == Dl) && (sqrt_Dr^2 == Dr)
	N2 = reshape(permute(reshape(N, (sqrt_Dl, sqrt_Dl, sqrt_Dr, sqrt_Dr)), (1,3,2,4)), (sqrt_Dl*sqrt_Dr, sqrt_Dl*sqrt_Dr))

	# # make N2 hermitian
	# N2H = (N2 + N2') / 2
	# # make N2 positive
	# eigvalues, eigvectors = eigen(Hermitian(N2H))
	# # println("eigenvalues $eigvalues")

	# eigvalue_threshold = 1.0e-12
	# pos = findfirst(x->x >= eigvalue_threshold, eigvalues / norm(eigvalues))
	# isnothing(pos) && error("environment becomes trivial")
	# # pos = 1

	# eigvalues = eigvalues[pos:end]
	# # println("eigenvalues $eigvalues")
	# eigvectors = eigvectors[:, pos:end]
	# X = eigvectors * Diagonal(sqrt.(eigvalues))
	X = make_positive_util(N2)
	return permute(reshape(X', (size(X, 2), sqrt_Dl, sqrt_Dr) ), (2,1,3))
end

# fig11
function _regularize_aL_aR_X(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, X::AbstractArray{<:Number, 3})
	s1, s2, s3 = size(X)
	Q1, R = LinearAlgebra.qr!(copy(reshape(X, (s1*s2, s3))))
	L, Q2 = LinearAlgebra.lq!(copy(reshape(X, (s1, s2*s3))))
	invL = inv(L)
	invR = inv(R)
	@tensor Xt[1,3,5] := invL[1,2] * X[2,3,4] * invR[4,5]
	@tensor aLt[2,3,4] := L[1,2] * aL[1,3,4]
	@tensor aRt[3,4,1] := R[1,2] * aR[3,4,2]
	return aLt, aRt, Xt, invL, invR
end

# fig12
function _regularize_aL_aR(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}; normalize::Bool)
	# @tensor down[2,4,6] := aL[1,2,3] * X[1,4,5] * aR[3,6,5]
	# total_norm = norm(down)
	# println("total norm is $total_norm")

	l1, l2, l3 = size(aL)
	r1, r2, r3 = size(aR)
	Q1, R = LinearAlgebra.qr!(copy(reshape(aL, (l1*l2, l3))))
	L, Q2 = LinearAlgebra.lq!(copy(reshape(aR, (r1, r2*r3))))
	m = R * L
	u, s, v, err = tsvd!(m)
	normalize && LinearAlgebra.normalize!(s)
	# s ./= total_norm 
	sqrt_s = Diagonal(sqrt.(s))
	u = u * sqrt_s
	v = sqrt_s * v

	return reshape(Q1 * u, (l1, l2, size(u, 2))), reshape(v * Q2, (size(v, 1), r2, r3))
end
_regularize_aL_aR(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, X::AbstractArray{<:Number, 3}; normalize::Bool) = _regularize_aL_aR(
	aL, aR; normalize=normalize)

