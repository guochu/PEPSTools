abstract type TruncationScheme end


struct NoTruncation <: TruncationScheme end
notrunc() = NoTruncation()

struct TruncateDim <: TruncationScheme
	D::Int
end
truncdim(D::Int) = TruncateDim(D)


struct TruncateCutoff <: TruncationScheme
	ϵ::Float64
end

# reserve at least add_back singular values
struct TruncateDimCutoff <: TruncationScheme
	D::Int
	ϵ::Float64
	add_back::Int

function TruncateDimCutoff(D::Int, ϵ::Real, add_back::Int)
	new(D, convert(Float64, ϵ), min(add_back, D))
end
end
TruncateDimCutoff(;D::Int, ϵ::Real, add_back::Int=1) = TruncateDimCutoff(D, ϵ, add_back)
truncdimcutoff(D::Int, ϵ::Real; add_back::Int=1) = TruncateDimCutoff(D, ϵ, add_back)
truncdimcutoff(; D::Int, ϵ::Real, add_back::Int=1) = TruncateDimCutoff(D, ϵ, add_back)

_truncate!(v::AbstractVector{<:Real}, trunc::NoTruncation, p::Real=2) = v, 0.

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncateDim, p::Real=2)
	dtrunc = min(length(v), trunc.D)
	truncerr = norm(view(v, dtrunc+1:length(v)), p)
	resize!(v, dtrunc)
	return v, truncerr
end

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncateCutoff, p::Real=2)
	sca = norm(v, p)
	dtrunc = findlast(Base.Fix2(>, sca * trunc.ϵ), v)
	if isnothing(dtrunc)
		dtrunc = 0
	end
	return _truncate!(v, TruncateDim(dtrunc), p)
end

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncateDimCutoff, p::Real=2)
	sca = norm(v, p)
	dtrunc = findlast(Base.Fix2(>, sca * trunc.ϵ), v)
	if isnothing(dtrunc)
		dtrunc = trunc.add_back
	end
	v, err = _truncate!(v, TruncateDim(min(trunc.D, dtrunc)), p)
	return v, err / sca
end


const DefaultTruncation = TruncateDimCutoff(D=100, ϵ=1.0e-6)
