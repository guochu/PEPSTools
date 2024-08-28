abstract type Abstract1DTN{T<:Number} end

const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

scalartype(::Type{<:Abstract1DTN{T}}) where {T} = T
scalartype(t::Abstract1DTN) = scalartype(typeof(t))
Base.length(t::Abstract1DTN) = length(t.data)
Base.getindex(t::Abstract1DTN, i::ValidIndices) = getindex(t.data, i)
Base.setindex!(t::Abstract1DTN, v, i::ValidIndices) = setindex!(t.data, v, i)
Base.firstindex(x::Abstract1DTN) = firstindex(x.data)
Base.lastindex(t::Abstract1DTN) = lastindex(t.data)
Base.isempty(t::Abstract1DTN) = isempty(t.data)

scaling(x::Abstract1DTN) = x.scaling[]
setscaling!(x::Abstract1DTN, scaling::Real) = (x.scaling[] = scaling)

function LinearAlgebra.normalize!(x::Abstract1DTN)
	setscaling!(x, 1)
	return x
end


const MPSTensor{T} = AbstractArray{T, 3} where {T<:Number}
const MPOTensor{T} = AbstractArray{T, 4} where {T<:Number}

space_l(m::MPSTensor) = size(m, 1)
space_r(m::MPSTensor) = size(m, 3)
space_l(m::MPOTensor) = size(m, 1)
space_r(m::MPOTensor) = size(m, 3)
space_l(t::Abstract1DTN) = size(t[1], 1)
space_r(t::Abstract1DTN) = size(t[end], 3)

r_RR(psiA::Abstract1DTN, psiB::Abstract1DTN) = Matrix{promote_type(scalartype(psiA), scalartype(psiB))}(I, space_r(psiA), space_r(psiB))
r_RR(psi::Abstract1DTN) = r_RR(psi, psi)
l_LL(psiA::Abstract1DTN, psiB::Abstract1DTN) = Matrix{promote_type(scalartype(psiA), scalartype(psiB))}(I, space_l(psiA), space_l(psiB))
l_LL(psi::Abstract1DTN) = l_LL(psi, psi)


bond_dimension(psi::Abstract1DTN, bond::Int) = begin
	((bond >= 1) && (bond <= length(psi))) || throw(BoundsError())
	space_r(psi[bond])
end 
bond_dimensions(psi::Abstract1DTN) = [bond_dimension(psi, i) for i in 1:length(psi)]
bond_dimension(psi::Abstract1DTN) = maximum(bond_dimensions(psi))


abstract type MPSCompression end
struct SVDCompression{T <: TruncationScheme} <: MPSCompression 
	trunc::T
end
SVDCompression(; trunc::TruncationScheme=DefaultTruncation) = SVDCompression(trunc)
changeD(x::SVDCompression{TruncateDimCutoff}; D::Int) = SVDCompression(truncdimcutoff(D=D, ϵ=x.trunc.ϵ))

const AllowedInitGuesses = (:svd, :pre, :rand)

struct IterativeCompression <: MPSCompression
	D::Int 
	maxiter::Int 
	tol::Float64 
	initguess::Symbol 
	verbosity::Int 
end

function IterativeCompression(;D::Int=100, maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0)
    (initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
    return IterativeCompression(D, maxiter, tol, initguess, verbosity)
end
changeD(x::IterativeCompression; D::Int) = IterativeCompression(D=D, maxiter=x.maxiter, tol=x.tol, initguess=x.initguess, verbosity=x.verbosity) 



function isleftcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[1,2,-1]) * psij[1,2,-2]
	return isapprox(r, one(r); kwargs...) 
end
function isrightcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[-1,1,2]) * psij[-2,1,2]
	return isapprox(r, one(r); kwargs...) 
end
function isleftcanonical(psij::MPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1,3]) * psij[1,2,-2,3]
	return isapprox(r, one(r); kwargs...) 
end
function isrightcanonical(psij::MPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[-1,1,2,3]) * psij[-2,1,2,3]
	return isapprox(r, one(r); kwargs...) 
end