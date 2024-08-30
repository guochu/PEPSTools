

# function rdm1(peps::CanonicalPEPS, i::Int, j::Int, alg::SimpleUpdate)
# 	bond_left = QuantumSpins.diag(sqrt.(peps.Hbonds[i, j-1]))
# 	bond_right = QuantumSpins.diag(sqrt.(peps.Hbonds[i, j]))
# 	bond_up = QuantumSpins.diag(sqrt.(peps.Vbonds[i-1, j]))
# 	bond_down = QuantumSpins.diag(sqrt.(peps.Vbonds[i, j]))
# 	tmp = _absorb_bonds(peps.Γs[i, j], bond_left, bond_right, bond_up, bond_down)
# 	@tensor rho[6,1] := conj(tmp[1,2,3,4,5]) * tmp[6,2,3,4,5]
# 	return rho
# end


rdm1(peps::PEPS, i::Int, j::Int, alg::SimpleUpdate) = rdm1_util(peps[i, j])

function rdm1s(peps::PEPS, alg::SimpleUpdate)
	m, n = size(peps)
	r = Matrix{Matrix{scalartype(peps)}}(undef, size(peps))
	for i in 1:m
		for j in 1:n
			r[i, j] = rdm1_util(peps[i, j])
		end
	end
	return r
end

function rdm1_util(m::AbstractArray{<:Number, 5})
	@tensor rho[6,1] := conj(m[2,3,4,5,1]) * m[2,3,4,5,6]
	rho ./= tr(rho)
	return rho
end

rdm2s(peps::PEPS, alg::SimpleUpdate; kwargs...) = Dict("H"=>rdm2sH(peps, alg; kwargs...), "V"=>rdm2sV(peps, alg; kwargs...))


rdm2H(peps::PEPS, i::Int, j::Int, alg::SimpleUpdate) = rdm2H_util(peps[i, j], peps[i, j+1])

function rdm2sH(peps::PEPS, alg::SimpleUpdate; periodic::Bool=!is_nonperiodic(peps))
	m, n = size(peps)

	_ncols = periodic ? n : n-1
	rH = Matrix{Array{scalartype(peps), 4}}(undef, m, _ncols)
	for i in 1:m
		for j in 1:_ncols
			rH[i,j] = rdm2H_util(peps[i, j], peps[i, j+1])
		end
	end

	return rH
end


function rdm2H_util(m1::AbstractArray{<:Number, 5}, m2::AbstractArray{<:Number, 5})
	@tensor rho[6,12,1,8] := (conj(m1[2,3,4,5,1]) * m1[2,3,7,5,6]) * conj(m2[4,9,10,11,8]) * m2[7,9,10,11,12]
	return normalize_rho!(rho)
end

rdm2V(peps::PEPS, i::Int, j::Int, alg::SimpleUpdate) = rdm2V_util(peps[i, j], peps[i+1, j])

function rdm2sV(peps::PEPS, alg::SimpleUpdate; periodic::Bool=!is_nonperiodic(peps))
	m, n = size(peps)

	_nrows = periodic ? m : m-1
	rV = Matrix{Array{scalartype(peps), 4}}(undef, _nrows, n)
	for i in 1:_nrows
		for j in 1:n
			rV[i,j] = rdm2V_util(peps[i, j], peps[i+1, j])
		end
	end

	return rV
end

function rdm2V_util(m1::AbstractArray{<:Number, 5}, m2::AbstractArray{<:Number, 5})
	@tensor rho[6,12,1,8] := (conj(m1[2,3,4,5,1]) * m1[2,3,4,7,6]) * conj(m2[9,5,10,11,8]) * m2[9,7,10,11,12]
	return normalize_rho!(rho)
end



