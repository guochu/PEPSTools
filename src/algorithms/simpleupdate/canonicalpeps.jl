


struct CanonicalPEPS{T<:Number, R<:Real} <: AbstractPEPS{T}
	data::PeriodicArray{Array{T, 5},2}
	Hbonds::PeriodicArray{Vector{R},2}
	Vbonds::PeriodicArray{Vector{R},2}
end

function CanonicalPEPS(Γs::PeriodicArray{Array{T, 5},2}) where {T <: Number}
	m, n = size(Γs)
	R = real(T)
	Hbonds = PeriodicArray{Vector{R},2}(undef, m, n)
	Vbonds = PeriodicArray{Vector{R},2}(undef, m, n)
	for i in 1:m
		for j in 1:n
			D4, D5 = size(Γs[i, j], 4), size(Γs[i, j], 5)
			Hbonds[i, j] = ones(R, D4) ./ sqrt(D4)
			Vbonds[i, j] = ones(R, D5) ./ sqrt(D5)
		end
	end
	return CanonicalPEPS(Γs, Hbonds, Vbonds)
end

CanonicalPEPS(peps::PEPS) = CanonicalPEPS(copy(peps.data))

function PEPS(x::CanonicalPEPS{T, R}) where {T, R}
	Hbonds = PeriodicArray(diagm.([sqrt.(item) for item in x.Hbonds.data]))
	Vbonds = PeriodicArray(diagm.([sqrt.(item) for item in x.Vbonds.data]))
	Γs = x.data
	m, n = size(x)
	r = PEPS(T, m, n)
	for i in 1:m
		for j in 1:n
			r[i, j] = _absorb_bonds(Γs[i, j], Hbonds[i, j-1], Hbonds[i, j], Vbonds[i-1, j], Vbonds[i, j])
		end
	end
	return r
end

function _absorb_bonds(m::AbstractArray{<:Number, 5}, l::AbstractMatrix, r::AbstractMatrix, u::AbstractMatrix, d::AbstractMatrix)
	@tensor tmp[1,6,7,8,9] := (((m[1,2,3,4,5] * l[6,2]) * u[7,3]) * r[4,8]) * d[5,9]
	return tmp
end