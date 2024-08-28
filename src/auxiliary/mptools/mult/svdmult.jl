
function svdmult(mpo::MPO, mps::MPS, trunc::TruncationScheme)
	(length(mpo) != length(mps)) && throw(DimensionMismatch("mpo and mps should have the same size"))
	isempty(mpo) && error("mpo is empty")
	L = length(mps)
	T = promote_type(scalartype(mpo), scalartype(mps))
	res = Vector{Array{T, 3}}(undef, L)

	@tensor tmp[1,5,2,3,6] := mpo[1][1,2,3,4] * mps[1][5,4,6]
	tmp3 = tie(tmp, (2,1,2))
	workspace = T[]
	q, r = tqr!(tmp3, (1,2), (3,), workspace)
	res[1] = q

	for i in 2:L-1
	    @tensor tmp[1,5,2,3,6] := mpo[i][1,2,3,4] * mps[i][5,4,6]
	    tmp3 = tie(tmp, (2,1,2))
	    @tensor tmp3c[1,3,4] := r[1,2] * tmp3[2,3,4]
	    q, r = tqr!(tmp3c, (1,2), (3,), workspace)
	    res[i] = q
	end
	i = L
	@tensor tmp[1,5,2,3,6] := mpo[i][1,2,3,4] * mps[i][5,4,6]
	tmp3 = tie(tmp, (2,1,2))
	@tensor tmp3c[1,3,4] := r[1,2] * tmp3[2,3,4]
	res[L] = tmp3c
	z = MPS(res, scaling = scaling(mpo) * scaling(mps))
	return rightorth!(z, workspace, alg=Orthogonalize(SVD(), trunc, normalize=false))
end
