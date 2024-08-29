

# struct BlockBPSimpleExp{M<:BoundaryMPS} <: AbstractBlockBPPEPSUpdateAlgorithm
# 	block_size::Tuple{Int, Int} 
# 	msg_tol::Float64 
# 	msg_maxiter::Int
# 	msg_D::Int 
# 	verbosity::Int 
# 	update_alg::M 
# end

# BlockBPSimpleExp(; block_size::Tuple{Int, Int}=(3,3), msg_maxiter::Int=10, msg_tol::Real=1.0e-5, verbosity::Int=1, 
# 	update_alg::BoundaryMPS=BoundaryMPS(D2=3, D1=2*D2^2, verbosity=verbosity), msg_D=(update_alg.D2)^2) = BlockBPSimpleExp(
# 	block_size, msg_tol, msg_maxiter, msg_D, verbosity, update_alg)
# BlockBPSimpleExp(x::AbstractBlockBPPEPSUpdateAlgorithm) = BlockBPSimpleExp(x.block_size, x.msg_tol, x.msg_maxiter, x.msg_D, x.verbosity, x.update_alg)

# struct BlockBPCentralExp{M<:BoundaryMPS} <: AbstractBlockBPPEPSUpdateAlgorithm
# 	block_size::Tuple{Int, Int} 
# 	msg_tol::Float64 
# 	msg_maxiter::Int
# 	msg_D::Int 
# 	verbosity::Int 
# 	update_alg::M 
# end

# BlockBPCentralExp(; block_size::Tuple{Int, Int}=(3,3), msg_maxiter::Int=10, msg_tol::Real=1.0e-5, verbosity::Int=1, 
# 	update_alg::BoundaryMPS=BoundaryMPS(D2=3, D1=2*D2^2, verbosity=verbosity), msg_D=(update_alg.D2)^2) = BlockBPCentralExp(
# 	block_size, msg_tol, msg_maxiter, msg_D, verbosity, update_alg)
# BlockBPCentralExp(x::AbstractBlockBPPEPSUpdateAlgorithm) = BlockBPCentralExp(x.block_size, x.msg_tol, x.msg_maxiter, x.msg_D, x.verbosity, x.update_alg)



expectation(U::SquareLatticeHamiltonian, peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm) = expectation(squeeze(U), peps, alg)
# expectation(U::SquareLatticeOperator, peps::PEPS, alg::BlockBPSimpleExp) = expectation(default_splitting(U, alg.block_size), peps, alg)
# expectation(U::SquareLatticeOperator, peps::PEPS, alg::BlockBPCentralExp) = expectation(center_splitting(U, alg.block_size), peps, alg)
expectation(U::SquareLatticeOperator, peps::PEPS, alg::BlockBP) = expectation(center_splitting(U, alg.block_size), peps, alg)


function expectation(Us::Vector{<:BlockOperator}, peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm)
	r = 0.
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		r += expectation(U, blk, alg)
	end
	return r
end

function expectation(U::BlockOperator, blk::BeliefPEPSBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm)
	@assert blk.partition == U.partition
	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	r = 0.
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r += expectation(sU, x, mult_alg)
		end
	end	
	return r
end

function expectationfull(U::BlockOperator, blk::BeliefPEPSBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm)
	@assert blk.partition == U.partition
	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)

	rH = PeriodicArray(zeros(scalartype(blk), size(blk)))
	rV = PeriodicArray(zeros(scalartype(blk), size(blk)))

	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r = expectationfull(sU, x, mult_alg)
			rH[rowindices(blk, i), colindices(blk, j)] .+= r["H"]
			rV[rowindices(blk, i), colindices(blk, j)] .+= r["V"]
		end
	end	
	return Dict("H"=>rH.data, "V"=>rV.data)
end

local_expectations(U::LocalQObservers, peps::PEPS, alg::BlockBP) = local_expectations(center_splitting(U, alg.block_size), peps, alg)


function local_expectations(Us::Vector{BlockLocalOperator{M}}, peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm) where {M<:AbstractMatrix}
	r = zeros(scalartype(peps), size(peps))
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		r += local_expectations(U, blk, alg)
	end
	return r
end

function local_expectations(U::BlockLocalOperator{M}, blk::BeliefPEPSBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm) where {M<:AbstractMatrix}
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	mult_alg = get_msg_mult_alg(alg)
	r = PeriodicArray(zeros(scalartype(blk), size(blk.peps)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)

			# r[rowindices(blk, i), colindices(blk, j)] = local_expectations(sU, x, mult_alg)

			tmp = local_expectations(sU, x, mult_alg)
			row_pos, row_idx = rowcenterenumerate(blk.partition, i)
			col_pos, col_idx = colcenterenumerate(blk.partition, j)

			r[row_idx, col_idx] = tmp[row_pos, col_pos]
		end
	end	
	return r.data
end

function rdm1s(peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm)
	U = rdm1_trivial_operator(size(peps))
	return rdm1s_util(center_splitting(U, alg.block_size), peps, alg)
end

function rdm1s_util(Us::Vector{BlockLocalOperator{Int}}, peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm) 
	T = scalartype(peps)
	r = Matrix{Union{Matrix{T}, Nothing}}(nothing, size(peps))
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		# r += local_expectations(U, blk, alg)
		r = _merge_rdms!(r, rdm1s_single_block(U, blk, alg))
	end
	return r
end

function rdm1s_single_block(U::BlockLocalOperator{Int}, blk::BeliefPEPSBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm) 
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	mult_alg = get_msg_mult_alg(alg)
	T = scalartype(blk)
	r = PeriodicArray(Matrix{Union{Matrix{T}, Nothing}}(nothing, size(blk)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r[rowindices(blk, i), colindices(blk, j)] = _rdm1s(sU, x, mult_alg)
		end
	end	
	return r.data
end

function _rdm1s(U::SquareLatticeSites{Int}, blk::BoundaryPEPS, alg::MPSCompression)
	m, n = size(blk)

	rH = Matrix{Union{Matrix{scalartype(blk)}, Nothing}}(nothing, size(blk))

	if n_nontrivial_terms(U) > 0
		mpsstorage = compute_H_mpsstorages(blk, alg)
		up = up_boundary(blk)
		for i in 1:m
			row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
			rH[i, :] = _rdm1s(U[i, :], row_i)
			if i != m
				mpo = mpoup(blk, i) 
				up, err = mpompsmult(mpo, up, alg)
				normalize!(up)
			end
		end	
	end	

	return rH
end

function _rdm1s(Us::Vector, x::DoubleLayerSandwichEnv)
	@assert length(x) == length(Us)
	return [_row_rdm1_single(x, pos, Us[pos]) for pos in 1:length(x)]
end 
function _row_rdm1_single(x::DoubleLayerSandwichEnv, pos::Int, m::Nothing)
	update_left!(x, pos)
	return nothing
end 
_row_rdm1_single(x::DoubleLayerSandwichEnv, pos::Int, m::Int) = row_rdm1_single(x, pos)

function rdm1_trivial_operator(shape::Tuple{Int, Int})
	data = Matrix{Union{Int, Nothing}}(nothing, shape)
	m, n = shape
	for i in 1:m
		for j in 1:n
			data[i, j] = 1
		end
	end
	return SquareLatticeSites(PeriodicArray(data))
end

function rdm2s(peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm; periodic::Bool=!is_nonperiodic(peps))
	U = rdm2_trivial_operator(size(peps), periodic)
	T = scalartype(peps)
	r = rdm2s_util(center_splitting(U, alg.block_size), peps, alg)
	rH = r["H"]
	rV = r["V"]
	if !periodic
		rH = rH[:, 1:end-1]
		rV = rV[1:end-1, :]
	end
	rH = convert(Matrix{Array{T, 4}}, rH)
	rV = convert(Matrix{Array{T, 4}}, rV)
	return Dict("H"=>rH, "V"=>rV)
end

function rdm2s_util(Us::Vector{BlockOperator{Int}}, peps::PEPS, alg::AbstractBlockBPPEPSUpdateAlgorithm) 
	T = scalartype(peps)
	rH = Matrix{Union{Array{T, 4}, Nothing}}(nothing, size(peps))
	rV = Matrix{Union{Array{T, 4}, Nothing}}(nothing, size(peps))
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		# r += local_expectations(U, blk, alg)
		r = rdm2s_single_block(U, blk, alg)
		rH = _merge_rdms!(rH, r["H"])
		rV = _merge_rdms!(rV, r["V"])
	end
	return Dict("H"=>rH, "V"=>rV)
end

function rdm2s_single_block(U::BlockOperator, blk::BeliefPEPSBlock, alg::AbstractBlockBPPEPSUpdateAlgorithm)
	@assert blk.partition == U.partition
	mult_alg = get_msg_mult_alg(alg)
	compute_messages!(blk, alg)
	T = scalartype(blk)
	rH = PeriodicArray(Matrix{Union{Array{T, 4}, Nothing}}(nothing, size(blk)) )
	rV = PeriodicArray(Matrix{Union{Array{T, 4}, Nothing}}(nothing, size(blk)) )

	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock(U, i, j)
			r = _rdm2s(sU, x, mult_alg)
			rH[rowindices(blk, i), colindices(blk, j)] = r["H"]
			rV[rowindices(blk, i), colindices(blk, j)] = r["V"]
		end
	end	
	return Dict("H"=>rH.data, "V"=>rV.data)
end

_rdm2s(U::SquareLatticeBonds{Union{Int, Nothing}}, blk::BoundaryPEPS, alg::MPSCompression) = Dict(
	"H"=>_rdm2sH(U.H, blk, alg), "V"=>_rdm2sV(U.V, blk, alg))


function _rdm2sH(H, blk::BoundaryPEPS, alg::MPSCompression)
	m, n = size(blk)

	rH = Matrix{Union{Array{scalartype(blk), 4}, Nothing}}(nothing, size(blk))

	if n_nontrivial_terms(H) > 0
		mpsstorage = compute_H_mpsstorages(blk, alg)
		up = up_boundary(blk)
		for i in 1:m
			row_i = row_environments(up, row_peps(blk, i), mpsstorage[i+1], blk.left[i], blk.right[i]) 
			rH[i, 1:n-1] = _rdm2s(H[i, 1:n-1], row_i)
			if i != m
				mpo = mpoup(blk, i) 
				up, err = mpompsmult(mpo, up, alg)
				normalize!(up)
			end
		end		
	end

	return rH
end


function _rdm2sV(V, blk::BoundaryPEPS, alg::MPSCompression)
	m, n = size(blk)

	rV= Matrix{Union{Array{scalartype(blk), 4}, Nothing}}(nothing, size(blk))
	
	if n_nontrivial_terms(V) > 0
		mpsstorage = compute_V_mpsstorages(blk, alg)
		left = left_boundary(blk)
		for i in 1:n
			# println("expectation of the $i-th column...")
			row_i = row_environments(mpsstorage[i+1], col_peps_as_row(blk, i), left, permute(blk.up[i], (3,2,1)), permute(blk.down[i], (3,2,1)))
			rV[1:m-1, i] = _rdm2s(V[1:m-1, i], row_i)
			if i != n
				mpo = mpoleft(blk, i)
				left, err = mpompsmult(mpo, left, alg)
				normalize!(left)
			end
		end		
	end

	return rV
end

function _rdm2s(Us::Vector, x::DoubleLayerSandwichEnv)
	@assert length(x) == length(Us) + 1
	return [_row_rdm2_single(x, pos, Us[pos]) for pos in 1:length(Us)]
end 


function _row_rdm2_single(x::DoubleLayerSandwichEnv, pos::Int, m::Nothing)
	update_left!(x, pos)
	return nothing
end
_row_rdm2_single(x::DoubleLayerSandwichEnv, pos::Int, m::Int) = row_rdm2_single(x, pos)


function rdm2_trivial_operator(shape::Tuple{Int, Int}, periodic::Bool)
	V = Matrix{Union{Int, Nothing}}(nothing, shape)
	H = Matrix{Union{Int, Nothing}}(nothing, shape)
	m, n = shape
	_ncols = periodic ? n : n-1
	for i in 1:m
		for j in 1:_ncols
			H[i, j] = 1
		end
	end
	_nrows = periodic ? m : m-1
	for i in 1:_nrows
		for j in 1:n
			V[i, j] = 1
		end
	end
	return SquareLatticeBonds(V, H)
end


function _merge_rdms!(a::Matrix{Union{M, Nothing}}, b::Matrix{Union{M, Nothing}}) where {M <: AbstractArray}
	@assert size(a) == size(b)
	m, n = size(a)
	for i in 1:m
		for j in 1:n
			bj = b[i, j]
			if !isnothing(bj)
				isnothing(a[i,j]) || error("Element ($i, $j) are notrivial for both matrices.")
				a[i, j] = bj
			end
		end
	end
	return a
end

