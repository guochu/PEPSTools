# orthogonalize mps to be left-canonical or right-canonical
abstract type MatrixProductOrthogonalAlgorithm  end

"""
	struct MatrixProductOrthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}
"""
struct Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
	orth::A
	trunc::T
	normalize::Bool
end
Orthogonalize(a::Union{QR, SVD}, trunc::TruncationScheme; normalize::Bool=false) = Orthogonalize(a, trunc, normalize)
Orthogonalize(a::Union{QR, SVD}; trunc::TruncationScheme=NoTruncation(), normalize::Bool=false) = Orthogonalize(a, trunc, normalize)
Orthogonalize(; alg::Union{QR, SVD} = SVD(), trunc::TruncationScheme=NoTruncation(), normalize::Bool=false) = Orthogonalize(alg, trunc, normalize)


leftorth!(psi::MPS, workspace::AbstractVector=Vector{scalartype(psi)}(); alg::Orthogonalize = Orthogonalize()) = _leftorth!(psi, alg.orth, alg.trunc, alg.normalize, workspace)
function _leftorth!(psi::MPS, alg::QR, trunc::TruncationScheme, normalize::Bool, workspace::AbstractVector=Vector{scalartype(psi)}())
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in 1:L-1
		q, r = tqr!(psi[i], (1, 2), (3,), workspace)
		_renormalize!(psi, r, normalize)
		psi[i] = q
		psi[i+1] = reshape(r * tie(psi[i+1], (1, 2)), size(r, 1), size(psi[i+1], 2), size(psi[i+1], 3))
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi, 0.
end

function _leftorth!(psi::MPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, workspace::AbstractVector=Vector{scalartype(psi)}())
	L = length(psi)
	errs = 0.
	for i in 1:L-1
		# u, s, v, err = stable_tsvd(psi[i], (1, 2), (3,), trunc=trunc)
		u, s, v, err = tsvd!(tie(psi[i], (2, 1)), workspace, trunc=trunc)
		_renormalize!(psi, s, normalize)
		d = length(s)
		psi[i] = reshape(u, size(psi[i],1), size(psi[i],2),d)
		v = Diagonal(s) * v
		psi[i+1] = reshape(v * tie(psi[i+1], (1, 2)), d, size(psi[i+1], 2), size(psi[i+1], 3))
		psi.s[i+1] = s
		errs = max(errs, err)
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi, errs
end

rightorth!(psi::MPS, workspace::AbstractVector=Vector{scalartype(psi)}(); alg::Orthogonalize = Orthogonalize()) = _rightorth!(psi, alg.orth, alg.trunc, alg.normalize, workspace)
function _rightorth!(psi::MPS, alg::QR, trunc::TruncationScheme, normalize::Bool, workspace::AbstractVector=Vector{scalartype(psi)}())
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in L:-1:2
		l, q = tlq!(psi[i], (1,), (2, 3), workspace)
		_renormalize!(psi, l, normalize)
		psi[i] = q
		psi[i-1] = reshape(tie(psi[i-1], (2, 1)) * l, size(psi[i-1], 1), size(psi[i-1], 2), size(l, 2))
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi, 0.
end

function _rightorth!(psi::MPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, workspace::AbstractVector=Vector{scalartype(psi)}())
	L = length(psi)
	errs = 0.
	for i in L:-1:2
		# u, s, v, err = stable_tsvd(psi[i], (1,), (2, 3), trunc=trunc)
		u, s, v, err = tsvd!(tie(psi[i], (1, 2)), workspace, trunc=trunc)
		_renormalize!(psi, s, normalize)
		d = length(s)
		psi[i] = reshape(v, d, size(psi[i], 2), size(psi[i], 3))
		u = u * Diagonal(s)
		psi[i-1] = reshape(tie(psi[i-1], (2, 1)) * u, size(psi[i-1], 1), size(psi[i-1], 2), d)
		psi.s[i] = s
		errs = max(errs, err)
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi, errs
end

function rightcanonicalize!(psi::MPS, workspace::AbstractVector=Vector{scalartype(psi)}(); alg::Orthogonalize = Orthogonalize(SVD(), DefaultTruncation, normalize=false))
	_leftorth!(psi, QR(), NoTruncation(), alg.normalize, workspace)
	return rightorth!(psi, workspace, alg=alg)
end

canonicalize!(psi::MPS, args...; kwargs...) = rightcanonicalize!(psi, args...; kwargs...)


function _rescaling!(psi::MPS, n::Real)
	L = length(psi)
	scale1 = n^(1/L)
	setscaling!(psi, scaling(psi) * scale1)
	return psi
end
function _rescaling!(psi::MPS)
	nrm1 = norm(psi[1])
	psi[1] = rmul!(psi[1], 1/nrm1)
	return _rescaling!(psi, nrm1)
end


function _renormalize!(psi::MPS, r, normalize::Bool)
	nr = norm(r)
	if nr != zero(nr)
		if !normalize
			_rescaling!(psi, nr)
		end
		lmul!(1/nr, r)  
  	end
end

function _renormalize_coeff!(psi::MPS, normalize::Bool)
	normalize && setscaling!(psi, 1)
end
