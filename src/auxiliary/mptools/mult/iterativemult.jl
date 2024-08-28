struct MPOMPSIterativeMultCache{_MPO, _IMPS, _OMPS, _H} 
	mpo::_MPO
	imps::_IMPS
	omps::_OMPS
	hstorage::_H
end
scalartype(::Type{MPOMPSIterativeMultCache{_MPO, _IMPS, _OMPS, _H}}) where {_MPO, _IMPS, _OMPS, _H} = promote_type(scalartype(_MPO), scalartype(_IMPS), scalartype(_OMPS))
scalartype(x::MPOMPSIterativeMultCache) = scalartype(typeof(x))

sweep!(m::MPOMPSIterativeMultCache, alg::IterativeCompression, workspace = scalartype(m)[]) = vcat(_leftsweep!(m, alg, workspace), _rightsweep!(m, alg, workspace))

function compute!(m::MPOMPSIterativeMultCache, alg::IterativeCompression, workspace = scalartype(m)[])
    kvals = Float64[]
    iter = 0
    tol = 1.
    while (iter < alg.maxiter) && (tol >= alg.tol)
        _kvals = sweep!(m, alg, workspace)
        tol = iterative_error_2(_kvals)
        push!(kvals, tol)
        iter += 1
        (alg.verbosity > 1) && println("finish the $iter-th sweep with error $tol", "\n")
    end
    if (alg.verbosity >= 2) && (iter < alg.maxiter)
        println("early converge in $iter-th sweeps with error $tol")
    end
    if (alg.verbosity > 0) && (tol >= alg.tol)
        println("fail to converge, required precision: $(alg.tol), actual precision $tol in $iter sweeps")
    end
    return kvals
end

function iterativemult(mpo::MPO, mps::MPS, alg::IterativeCompression = IterativeCompression())
    if alg.initguess == :svd
        mpsout = _svd_guess(mpo, mps, alg.D)
    elseif alg.initguess == :rand
        mpsout = randommps(promote_type(scalartype(mpo), scalartype(mps)), ophysical_dimensions(mpo), D=alg.D)
    elseif alg.initguess == :pre
        mpsout = increase_bond!(copy(mps), alg.D)
    else
        error("unsupported initguess $(alg.initguess)")
    end
	# canonicalize!(mpsout, normalize=true)
    rightorth!(mpsout, alg=Orthogonalize(normalize=true))
	m = MPOMPSIterativeMultCache(mpo, mps, mpsout, init_hstorage_right(mpsout, mpo, mps))
	kvals = compute!(m, alg)
	# return m.omps, kvals[end]
    z = m.omps
    setscaling!(z, scaling(mpo) * scaling(mps))
    _rescaling!(z)
    return z, kvals[end]
end


function _leftsweep!(m::MPOMPSIterativeMultCache, alg::IterativeCompression, workspace = scalartype(m)[])
    mpsA = m.imps
    mpo = m.mpo
    mpsB = m.omps
    Cstorage = m.hstorage
    L = length(mpo)
    kvals = Float64[]
    for site in 1:L-1
        (alg.verbosity > 3) && println("Sweeping from left to right at bond: $site.")
        mpsj = reduceD_single_site(mpsA[site], mpo[site], Cstorage[site], Cstorage[site+1])
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
		mpsB[site], r = tqr!(mpsj, (1,2), (3,), workspace)
        Cstorage[site+1] = updateleft(Cstorage[site], mpsB[site], mpo[site], mpsA[site])
    end
    return kvals	
end

function _rightsweep!(m::MPOMPSIterativeMultCache, alg::IterativeCompression, workspace = scalartype(m)[])
    mpsA = m.imps
    mpo = m.mpo
    mpsB = m.omps
    Cstorage = m.hstorage
    L = length(mpo)
    kvals = Float64[]
    r = zeros(scalartype(mpsB), 0, 0)

    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at bond: $site.")
        mpsj = reduceD_single_site(mpsA[site], mpo[site], Cstorage[site], Cstorage[site+1])
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        r, mpsB[site] = tlq!(mpsj, (1,), (2,3), workspace)
        Cstorage[site] = updateright(Cstorage[site+1], mpsB[site], mpo[site], mpsA[site])
    end
    # println("norm of r is $(norm(r))")
    @tensor tmp[1,2,4] := mpsB[1][1,2,3] * r[3,4]
    mpsB[1] = tmp
    return kvals	
end


function reduceD_single_site(A::AbstractArray{<:Number, 3}, X::AbstractArray{<:Number, 4}, Cleft::AbstractArray{<:Number, 3}, 
	Cright::AbstractArray{<:Number, 3})
    @tensor tmp[1,6,8] := ((Cleft[1,2,3] * A[3,4,5]) * X[2,6,7,4]) * Cright[8,7,5]
    return tmp
end

# function reduceD_single_bond(A::AbstractArray{<:Number, 4}, X1::AbstractArray{<:Number, 4}, X2::AbstractArray{<:Number, 4}, 
# 	Cleft::AbstractArray{<:Number, 3}, Cright::AbstractArray{<:Number, 3})
# 	@tensor temp[1,2,4,5,6] := Cleft[1,2,3] * A[3,4,5,6]
# 	@tensor temp2[1,6,7,4,5] := temp[1,2,3,4,5] * X1[2,6,7,3]
# 	@tensor temp3[1,2,6,7,5] := temp2[1,2,3,4,5] * X2[3,6,7,4]
# 	@tensor temp4[1,2,3,6] := temp3[1,2,3,4,5] * Cright[6,4,5]
# 	return temp4
# end

function init_hstorage_right(mpsB, mpo, mpsA)
	T = promote_type(scalartype(mpsA), scalartype(mpo), scalartype(mpsB))
	L = length(mpo)
	hstorage = Vector{Array{T, 3}}(undef, L+1)
	hstorage[1] = ones(T, 1, 1, 1)
	hstorage[L+1] = ones(T, 1, 1, 1)
	for i in L:-1:2
		hstorage[i] = updateright(hstorage[i+1], mpsB[i], mpo[i], mpsA[i])
	end
	return hstorage
end

function _svd_guess(mpo::MPO, mps::MPS, D::Int)
    (length(mpo) != length(mps)) && throw(DimensionMismatch("mpo and mps should have the same size"))
    isempty(mpo) && error("mpo is empty")
    L = length(mps)
    T = promote_type(scalartype(mpo), scalartype(mps))
    res = Vector{Array{T, 3}}(undef, L)

    @tensor tmp[1,5,2,3,6] := mpo[1][1,2,3,4] * mps[1][5,4,6]
    tmp3 = tie(tmp, (2,1,2))
    workspace = T[]
    trunc = truncdim(D)
    u, s, v, err = tsvd!(tmp3, (1,2), (3,), workspace, trunc=trunc)
    res[1] = u
    v = Diagonal(s) * v

    for i in 2:L-1
        @tensor tmp[1,5,2,3,6] := mpo[i][1,2,3,4] * mps[i][5,4,6]
        tmp3 = tie(tmp, (2,1,2))
        @tensor tmp3c[1,3,4] := v[1,2] * tmp3[2,3,4]
        u, s, v, err = tsvd!(tmp3c, (1,2), (3,), workspace, trunc=trunc)
        res[i] = u
        v = Diagonal(s) * v
    end
    i = L
    @tensor tmp[1,5,2,3,6] := mpo[i][1,2,3,4] * mps[i][5,4,6]
    tmp3 = tie(tmp, (2,1,2))
    @tensor tmp3c[1,3,4] := v[1,2] * tmp3[2,3,4]
    res[L] = tmp3c
    return MPS(res)
end

function iterative_error_2(kvals::AbstractVector{<:Real})
    aver = sum(kvals) / length(kvals)
    std_devi = sum(x->(x-aver)^2, kvals)
    return std_devi / abs(aver)
end
