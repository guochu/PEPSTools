

@with_kw struct PEPSSimpleUpdate <: AbstractPEPSUpdateAlgorithm
	D2::Int = 3
	ϵ::Float64 = 1.0e-8
	als_tol::Float64 = 1.0e-6
	als_maxiter::Int = 5
	verbosity::Int = 1
end

QuantumSpins.expectation(U::Union{SquareLatticeOperatorBase, SquareLatticeHamiltonianBase}, peps::CanonicalPEPS, alg::AbstractPEPSUpdateAlgorithm) = expectation(
	U, PEPS(peps), alg)


function QuantumSpins.sweep!(peps::CanonicalPEPS, U::SquareLatticeOperatorBase, alg::PEPSSimpleUpdate)
	@assert size(peps) == size(U)
	m, n = size(peps)
	trunc = MPSTruncation(D=alg.D2, ϵ=alg.ϵ)
	Γs = peps.Γs
	Hbonds = peps.Hbonds
	Vbonds = peps.Vbonds
	# update all the horizontal terms
	for i in 1:size(U,1)
		for j in 1:size(U,2)
			Γs[i, j], Hbonds[i, j], Γs[i, j+1] = evolve_h_single(U.H[i, j], Hbonds[i,j-1],Vbonds[i-1,j],Vbonds[i,j], Γs[i, j],Hbonds[i,j], Γs[i, j+1],
				Vbonds[i-1,j+1],Hbonds[i,j+1],Vbonds[i,j+1],
				trunc=trunc, normalize=true, maxiter=alg.als_maxiter, tol=alg.als_tol, verbosity=alg.verbosity)
		end
	end

	# update all the vertical terms
	for i in 1:size(U,1)
		for j in 1:size(U,2)
			Γs[i, j],Vbonds[i,j], Γs[i+1, j] = evolve_v_single(U.V[i, j],Hbonds[i,j-1],Vbonds[i-1,j],Hbonds[i,j], Γs[i, j],Vbonds[i,j], Γs[i+1, j], 
				Hbonds[i+1,j-1],Hbonds[i+1,j],Vbonds[i+1,j],
				trunc=trunc, normalize=true, maxiter=alg.als_maxiter, tol=alg.als_tol, verbosity=alg.verbosity)
		end
	end
end



function evolve_h_single(U::AbstractArray{<:Number, 4}, sLl, sLu, sLd, AL, sM, AR, sRu, sRr, sRd; 
	normalize::Bool=true, trunc::TruncationScheme, maxiter::Int, tol::Real, verbosity::Int)
	invLl = QuantumSpins.diag(1 ./ sLl)
	invLu = QuantumSpins.diag(1 ./ sLu)
	invLd = QuantumSpins.diag(1 ./ sLd)
	invRu = QuantumSpins.diag(1 ./ sRu)
	invRr = QuantumSpins.diag(1 ./ sRr)
	invRd = QuantumSpins.diag(1 ./ sRd)
	sLl = QuantumSpins.diag(sLl)
	sLu = QuantumSpins.diag(sLu)
	sLd = QuantumSpins.diag(sLd)
	sRr = QuantumSpins.diag(sRr)
	sRu = QuantumSpins.diag(sRu)
	sRd = QuantumSpins.diag(sRd)
	# sM = QuantumSpins.diag(sM)

	# ALn = _absorb_bonds(AL, sLl, sM, sLu, sLd) 
	@tensor ALn[1,6,7,4,8] := AL[1,2,3,4,5] * sLl[6,2] * sLu[7,3] * sLd[8,5]
	@tensor ARn[1,2,6,7,8] := AR[1,2,3,4,5] * sRu[6,3] * sRr[7,4] * sRd[8,5]

	left_q, aL = tqr!(ALn, (2,3,5), (1,4))
	aR, right_q = tlq!(ARn, (2,1), (3,4,5))	
	left_q = @tensor tmp[5,6,7,4] := (left_q[1,2,3,4] * invLl[5,1]) * invLu[6,2] * invLd[7,3]
	right_q = @tensor tmp[1,5,6,7] := ((right_q[1,2,3,4] * invRu[5,2]) * invRr[6,3]) * invRd[7,4]


	aLn, s, aRn = bond_evolve(aL, sM, aR, U; trunc=trunc)

	normalize && LinearAlgebra.normalize!(s)

	@tensor r1[5,1,2,6,3] := left_q[1,2,3,4] * aLn[4,5,6]
	@tensor r2[2,1,4,5,6] := aRn[1,2,3] * right_q[3,4,5,6]

	return r1, s, r2
end
evolve_h_single(U::Nothing,sLl, sLu, sLd, AL, sM, AR, sRu, sRr, sRd;kwargs...) = AL, sM, AR

function evolve_v_single(U::AbstractArray{<:Number, 4},sLl, sLu, sLr, AL, sM, AR, sRl, sRr, sRd;
	normalize::Bool=true, trunc::TruncationScheme, maxiter::Int, tol::Real, verbosity::Int)
	invLl = QuantumSpins.diag(1 ./ sLl)
	invLu = QuantumSpins.diag(1 ./ sLu)
	invLr = QuantumSpins.diag(1 ./ sLr)
	invRl = QuantumSpins.diag(1 ./ sRl)
	invRr = QuantumSpins.diag(1 ./ sRr)
	invRd = QuantumSpins.diag(1 ./ sRd)
	sLl = QuantumSpins.diag(sLl)
	sLu = QuantumSpins.diag(sLu)
	sLr = QuantumSpins.diag(sLr)
	sRl = QuantumSpins.diag(sRl)
	sRr = QuantumSpins.diag(sRr)
	sRd = QuantumSpins.diag(sRd)
	# sM = QuantumSpins.diag(sM)

	# ALn = _absorb_bonds(AL, sLl, sLr, sLu, sM)
	@tensor ALn[1,6,7,8,5] := AL[1,2,3,4,5] * sLl[6,2] * sLu[7,3] * sLr[8,4]
	@tensor ARn[1,6,3,7,8] := AR[1,2,3,4,5] * sRl[6,2] * sRr[7,4] * sRd[8,5]

	left_q, aL = tqr!(ALn, (2,3,4), (1,5))
	aR, right_q = tlq!(ARn, (3,1), (2,4,5))
	left_q = @tensor tmp[5,6,7,4] := left_q[1,2,3,4] * invLl[5,1] * invLu[6,2] * invLr[7,3]
	right_q = @tensor tmp[1,5,6,7] := right_q[1,2,3,4] * invRl[5,2] * invRr[6,3] * invRd[7,4]


	aLn, s, aRn = bond_evolve(aL, sM, aR, U; trunc=trunc)

	normalize && LinearAlgebra.normalize!(s)

	@tensor r1[5,1,2,3,6] := left_q[1,2,3,4] * aLn[4,5,6]
	@tensor r2[2,4,1,5,6] := aRn[1,2,3] * right_q[3,4,5,6]

	return r1, s, r2
end
evolve_v_single(U::Nothing,sLl, sLu, sLr, AL, sM, AR, sRl, sRr, sRd;kwargs...) = AL, sM, AR


function bond_evolve(aL::AbstractArray{<:Number, 3}, sM::AbstractVector, aR::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}; trunc::TruncationScheme)
	sM = QuantumSpins.diag(sM)
	@tensor m[1,6,7,5] := ((aL[1,2,3] * sM[3,8]) * aR[8,4,5]) * U[6,7,2,4]
	s1, s2, s3, s4 = size(m)
	u, s, v, err = tsvd!(reshape(m, (s1*s2, s3*s4)); trunc=trunc)
	Dn = length(s)
	return reshape(u, (s1, s2, Dn)), s, reshape(v, (Dn, s3, s4))
end

# function simple_update(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}; 
# 	normalize::Bool=true, trunc::TruncationScheme, maxiter::Int, tol::Real, verbosity::Int)

# 	# compute initial loss
# 	init_loss = compute_loss_simple(aL, aR, U, aL, aR)
# 	# println("initial loss is $init_loss")
	
# 	# compute initial guess for aLn and aRn
# 	aLn, aRn = bond_evolve(aL, aR, U, trunc=trunc)
# 	loss_svd = compute_loss_simple(aL, aR, U, aLn, aRn)
# 	# println("simple update loss is $loss_svd")

# 	iter = 1
# 	errs = [init_loss, loss_svd]
# 	# how to compute the loss here?
# 	while iter <= maxiter
# 		bl = compute_bl_simple(aL, aR, U, aRn)
# 		aLn, d1 = local_solve_simple_left(updateright_simple(aRn), bl)

# 		# push!(errs, d1)
# 		br = compute_br_simple(aL, aR, U, aLn)
# 		aRn, d2 = local_solve_simple_right(updateleft_simple(aLn), br)

# 		err = abs((d2 - errs[end]) / init_loss)

# 		push!(errs, d2)

# 		# println("relative error in $iter-th iteration is $err")

# 		if err <= tol
# 			(verbosity > 2) && println("early converge in $iter with relative error $err")
# 			break
# 		end

# 		iter += 1
# 	end
# 	# println("final loss is $(errs[end])")

# 	if (verbosity > 1) && (iter > maxiter)
# 		println("fail to converge to precision $tol in $maxiter iterations.")
# 	end

# 	return _regularize_aL_aR(aLn, aRn, normalize=normalize)
# end

# function compute_loss_simple(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}, 
# 	aLn::AbstractArray{<:Number, 3}, aRn::AbstractArray{<:Number, 3})
# 	@tensor center[1,6,7,5] := aL[1,2,3] * aR[3,4,5] * U[6,7,2,4]
# 	@tensor centern[1,2,4,5] := aLn[1,2,3] * aRn[3,4,5]
# 	return real(dot(center, center)) - 2 * real(dot(centern, center))
# end

# function compute_bl_simple(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, 
# 	U::AbstractArray{<:Number, 4}, aRn::AbstractArray{<:Number, 3})
# 	@tensor tmp[7,1,8,5] := (conj(aRn[1,2,3]) * aR[5,6,3]) * U[7,2,8,6]
# 	@tensor r[2,5,6] := aL[2,3,4] * tmp[5,6,3,4]
# 	return r
# end

# function compute_br_simple(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, 
# 	U::AbstractArray{<:Number, 4}, aLn::AbstractArray{<:Number, 3})
# 	@tensor tmp[4,7,6,8] := conj(aLn[2,3,4]) * aL[2,5,6] * U[3,7,5,8]
# 	@tensor r[1,2,6] := tmp[1,2,3,4] * aR[3,4,6] 
# 	return r
# end

# function local_solve_simple_left(right::AbstractMatrix, bl::AbstractArray{<:Number, 3})
# 	invright = inv(right)
# 	@tensor tmp[1,3,5] := bl[1,3,4] * invright[5,4]
# 	return tmp, _compute_distance(bl, tmp)
# end

# function local_solve_simple_right(left::AbstractMatrix, bl::AbstractArray{<:Number, 3})
# 	# println("conditional numbers are $(cond(left)), $(cond(right))")
# 	invleft = inv(left)
# 	@tensor tmp[1,3,5] := invleft[1,2] * bl[2,3,5] 
# 	return tmp, _compute_distance(bl, tmp)
# end

# function updateright_simple(a::AbstractArray{<:Number, 3})
# 	@tensor tmp[1,4] := conj(a[1,2,3]) * a[4,2,3]
# 	return tmp
# end

# function updateleft_simple(a::AbstractArray{<:Number, 3})
# 	@tensor tmp[3,4] := conj(a[1,2,3]) * a[1,2,4]
# 	return tmp
# end
# _compute_distance(b::AbstractArray{<:Number, 3}, x::AbstractArray{<:Number, 3}) = real(-dot(b, x))

# function bond_evolve(aL::AbstractArray{<:Number, 3}, aR::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}; trunc::TruncationScheme)
# 	@tensor m[1,6,7,5] := aL[1,2,3] * aR[3,4,5] * U[6,7,2,4]
# 	s1, s2, s3, s4 = size(m)
# 	u, s, v, err = tsvd!(reshape(m, (s1*s2, s3*s4)); trunc=trunc)
# 	Dn = length(s)
# 	sqrt_s = Diagonal(sqrt.(s))
# 	u = u * sqrt_s
# 	v = sqrt_s * v
# 	aLn = reshape(u, (s1, s2, Dn))
# 	aRn = reshape(v, (Dn, s3, s4))
# 	return aLn, aRn	
# end
