

struct PEPSSimpleUpdate <: AbstractPEPSUpdateAlgorithm
	D2::Int 
	ϵ::Float64 
	als_tol::Float64 
	als_maxiter::Int 
	verbosity::Int 
end
PEPSSimpleUpdate(; D2::Int=3, ϵ::Real=1.0e-8, als_tol::Real=1.0e-6, als_maxiter::Int=5, verbosity::Int=1) = PEPSSimpleUpdate(
				D2, convert(Float64, ϵ), convert(Float64, als_tol), als_maxiter, verbosity)

expectation(U::Union{SquareLatticeOperator, SquareLatticeHamiltonian}, peps::CanonicalPEPS, alg::AbstractPEPSUpdateAlgorithm) = expectation(
	U, PEPS(peps), alg)


function sweep!(peps::CanonicalPEPS, U::SquareLatticeOperator, alg::PEPSSimpleUpdate)
	@assert size(peps) == size(U)
	m, n = size(peps)
	trunc = truncdimcutoff(D=alg.D2, ϵ=alg.ϵ)
	Γs = peps.data
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
	invLl = diagm(1 ./ sLl)
	invLu = diagm(1 ./ sLu)
	invLd = diagm(1 ./ sLd)
	invRu = diagm(1 ./ sRu)
	invRr = diagm(1 ./ sRr)
	invRd = diagm(1 ./ sRd)
	sLl = diagm(sLl)
	sLu = diagm(sLu)
	sLd = diagm(sLd)
	sRr = diagm(sRr)
	sRu = diagm(sRu)
	sRd = diagm(sRd)

	# ALn = _absorb_bonds(AL, sLl, sM, sLu, sLd) 
	@tensor ALn[6,7,4,8,1] := AL[2,3,4,5,1] * sLl[6,2] * sLu[7,3] * sLd[8,5]
	@tensor ARn[2,6,7,8,1] := AR[2,3,4,5,1] * sRu[6,3] * sRr[7,4] * sRd[8,5]

	left_q, aL = tqr!(ALn, (1,2,4), (5,3))
	aR, right_q = tlq!(ARn, (1,5), (2,3,4))	
	left_q = @tensor tmp[5,6,7,4] := (left_q[1,2,3,4] * invLl[5,1]) * invLu[6,2] * invLd[7,3]
	right_q = @tensor tmp[1,5,6,7] := ((right_q[1,2,3,4] * invRu[5,2]) * invRr[6,3]) * invRd[7,4]


	aLn, s, aRn = bond_evolve(aL, sM, aR, U; trunc=trunc)

	normalize && LinearAlgebra.normalize!(s)

	@tensor r1[1,2,6,3,5] := left_q[1,2,3,4] * aLn[4,5,6]
	@tensor r2[1,4,5,6,2] := aRn[1,2,3] * right_q[3,4,5,6]

	return r1, s, r2
end
evolve_h_single(U::Nothing,sLl, sLu, sLd, AL, sM, AR, sRu, sRr, sRd;kwargs...) = AL, sM, AR

function evolve_v_single(U::AbstractArray{<:Number, 4},sLl, sLu, sLr, AL, sM, AR, sRl, sRr, sRd;
	normalize::Bool=true, trunc::TruncationScheme, maxiter::Int, tol::Real, verbosity::Int)
	invLl = diagm(1 ./ sLl)
	invLu = diagm(1 ./ sLu)
	invLr = diagm(1 ./ sLr)
	invRl = diagm(1 ./ sRl)
	invRr = diagm(1 ./ sRr)
	invRd = diagm(1 ./ sRd)
	sLl = diagm(sLl)
	sLu = diagm(sLu)
	sLr = diagm(sLr)
	sRl = diagm(sRl)
	sRr = diagm(sRr)
	sRd = diagm(sRd)

	# ALn = _absorb_bonds(AL, sLl, sLr, sLu, sM)
	@tensor ALn[6,7,8,5,1] := AL[2,3,4,5,1] * sLl[6,2] * sLu[7,3] * sLr[8,4]
	@tensor ARn[6,3,7,8,1] := AR[2,3,4,5,1] * sRl[6,2] * sRr[7,4] * sRd[8,5]

	left_q, aL = tqr!(ALn, (1,2,3), (5,4))
	aR, right_q = tlq!(ARn, (2,5), (1,3,4))
	left_q = @tensor tmp[5,6,7,4] := left_q[1,2,3,4] * invLl[5,1] * invLu[6,2] * invLr[7,3]
	right_q = @tensor tmp[1,5,6,7] := right_q[1,2,3,4] * invRl[5,2] * invRr[6,3] * invRd[7,4]


	aLn, s, aRn = bond_evolve(aL, sM, aR, U; trunc=trunc)

	normalize && LinearAlgebra.normalize!(s)

	@tensor r1[1,2,3,6,5] := left_q[1,2,3,4] * aLn[4,5,6]
	@tensor r2[4,1,5,6,2] := aRn[1,2,3] * right_q[3,4,5,6]

	return r1, s, r2
end
evolve_v_single(U::Nothing,sLl, sLu, sLr, AL, sM, AR, sRl, sRr, sRd;kwargs...) = AL, sM, AR


function bond_evolve(aL::AbstractArray{<:Number, 3}, sM::AbstractVector, aR::AbstractArray{<:Number, 3}, U::AbstractArray{<:Number, 4}; trunc::TruncationScheme)
	sM = diagm(sM)
	@tensor m[1,6,7,5] := ((aL[1,2,3] * sM[3,8]) * aR[8,4,5]) * U[6,7,2,4]
	s1, s2, s3, s4 = size(m)
	u, s, v, err = tsvd!(reshape(m, (s1*s2, s3*s4)); trunc=trunc)
	Dn = length(s)
	return reshape(u, (s1, s2, Dn)), s, reshape(v, (Dn, s3, s4))
end
