using Distributed

@everywhere push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")

@everywhere using QuantumSpins, PeriodicMPS, PEPSTools


function parallel_ground_state!(peps::PEPS, h::PEPSTools.SquareLatticeHamiltonian, alg::ImaginaryTimePEPS)
	isa(alg.stepper, BlockBP) || error("only BlockBP algorithm support parallelization currently")
	isa(alg.measure_alg, BlockBP) || println("Warning! Only BlockBP algorithm support parallelized measurement")

	sh = squeeze(h)
	ntail = 10
	expec_alg = alg.measure_alg
	verbosity = alg.verbosity
	# always measure initial energy and final energy
	energies = [real(parallel_expectation(sh, peps, expec_alg)) / prod(size(peps))]
	(verbosity > 2) && println("initial state energy is $(energies[1])")
	iter = 1
	err = 1.
	for (nsteps, δt, tol) in alg.parameters
		(verbosity > 2) && println("do at most $nsteps sweeps with δt=$(δt) and convergence tolerance $tol")
		(δt < 0.) || error("δt should be negative.")
		U = exponential(sh, δt)
		energies_tmp = Float64[]
		for j in 1:nsteps
			if verbosity > 2
				t = @elapsed parallel_sweep!(peps, U, alg.stepper)
				println("the $iter-th sweep takes $t seconds")
			else
				parallel_sweep!(peps, U, alg.stepper)
			end
			

			if (alg.sweeps_per_measure > 0) && (iter % alg.sweeps_per_measure == 0)
				expec_peps = parallel_expectation(sh, peps, expec_alg) / prod(size(peps))
				if verbosity > 2
					println("energy after $iter-th sweep is $expec_peps")
				end
				push!(energies_tmp, real(expec_peps))
				if length(energies_tmp) >= ntail
					err = QuantumSpins.iterative_error_2(energies_tmp[end-ntail+1:end])
					(verbosity > 2) &&println("error at the $iter-th sweep is $err")
					if err < tol
						(verbosity > 1) && println("early converge in $iter sweeps with error $err.")
						break
					end
				end
				(verbosity > 2) && println()
			end

			iter += 1
		end
		append!(energies, energies_tmp)
	end

	return energies, (iter, err)
end

parallel_expectation(U::PEPSTools.SquareLatticeHamiltonian, peps::PEPS, alg::BlockBP) = parallel_expectation(squeeze(U), peps, alg)
parallel_expectation(U::PEPSTools.SquareLatticeOperator, peps::PEPS, alg::BlockBP) = parallel_expectation(center_splitting(U, alg.block_size), peps, alg)

# this is not parallelized
parallel_expectation(U::PEPSTools.SquareLatticeHamiltonian, peps::PEPS, alg::BoundaryMPS) = parallel_expectation(squeeze(U), peps, alg)
parallel_expectation(U::PEPSTools.SquareLatticeOperator, peps::PEPS, alg::BoundaryMPS) = parallel_expectation(U, borderedpeps(peps), alg)
# use at most two processes
function parallel_expectation(U::PEPSTools.SquareLatticeOperator, blk::BorderedPEPS, alg::BoundaryMPS) 
	PEPSTools.is_nonperiodic(U) || error("BoundaryMPS does not support periodic boundary, using BlockBP instead.")
	@assert size(blk) == size(U)
	mult_alg = PEPSTools.get_mult_alg(alg)
	if nworkers() >= 2
		# println("parallel expec")
		vh = fetch(@spawnat 1 sum(PEPSTools._expect_H(U.H, blk, mult_alg)))
		vv = fetch(@spawnat 2 sum(PEPSTools._expect_V(U.V, blk, mult_alg)))
		return vh + vv
	else
		# println("serial expec")
		return expectation(U, blk, mult_alg)
	end
end


function parallel_expectation(Us::Vector{<:BlockOperator}, peps::PEPS, alg::BlockBP)
	r = 0.
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		r += parallel_expectation(U, blk, alg)
	end
	return r
end

function parallel_expectation(U::BlockOperator, blk::BeliefPEPSBlock, alg::BlockBP)
	@assert blk.partition == U.partition
	parallel_compute_messages!(blk, alg)
	mult_alg = PEPSTools.get_msg_mult_alg(alg)

	index = CartesianIndices((nrows(blk), ncols(blk)))

	n_rank = nworkers()
	# (length(index) % n_rank == 0) || println("total number of jobs $(length(index)) can not be divided by number of workers $(n_rank).")
	nave, resi, jobs = partition_jobs(length(index), n_rank)

	f_per_rank(l::Int) = begin
		r = 0.
		for k in jobs[l]+1:jobs[l+1]
			idx = index[k]
			i, j = idx[1], idx[2]
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = PEPSBlock(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock( U, i, j)
			r += expectation(sU, x, mult_alg)
		end
		return r
	end 

	r = pmap(f_per_rank, 1:n_rank)
	return sum(r)
end

parallel_local_expectations(U::PEPSTools.LocalObservers, peps::SquareTN, alg::BlockBP) = parallel_local_expectations(center_splitting(U, alg.block_size), peps, alg)


function parallel_local_expectations(Us::Vector{<:PEPSTools.BlockLocalOperator}, peps::SquareTN, alg::PEPSTools.AbstractBlockBPPEPSUpdateAlgorithm)
	r = zeros(scalartype(peps), size(peps))
	for U in Us
		blk = PEPSTools.BeliefSquareTNBlock(peps, U.partition)
		r += parallel_local_expectations(U, blk, alg)
	end
	return r
end

function parallel_local_expectations(U::PEPSTools.BlockLocalOperator, blk::PEPSTools.BeliefSquareTNBlock, alg::PEPSTools.AbstractBlockBPPEPSUpdateAlgorithm)
	@assert blk.partition == U.partition
	parallel_compute_messages!(blk, alg)
	mult_alg = PEPSTools.get_msg_mult_alg(alg)

	index = CartesianIndices((nrows(blk), ncols(blk)))

	n_rank = nworkers()
	# (length(index) % n_rank == 0) || println("total number of jobs $(length(index)) can not be divided by number of workers $(n_rank).")
	nave, resi, jobs = partition_jobs(length(index), n_rank)

	f_per_rank(l::Int) = begin
		r = PeriodicArray(zeros(scalartype(blk), size(blk.peps)))
		for k in jobs[l]+1:jobs[l+1]
			idx = index[k]
			i, j = idx[1], idx[2]
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock( U, i, j)
			r[PEPSTools.rowindices(blk, i), PEPSTools.colindices(blk, j)] = PEPSTools.local_expectations(sU, x, mult_alg)
		end
		return r.data
	end 

	r = pmap(f_per_rank, 1:n_rank)

	return sum(r)
end


parallel_sweep!(peps::PEPS, U::PEPSTools.SquareLatticeOperator, alg::BlockBP) = parallel_sweep!(peps, default_splitting(U, alg.block_size), alg)
function parallel_sweep!(peps::PEPS, Us::Vector{<:BlockOperator}, alg::BlockBP)
	for U in Us
		blk = BeliefPEPSBlock(peps, U.partition)
		parallel_sweep!(blk, U, alg)
	end
end
function parallel_sweep!(blk::BeliefPEPSBlock, U::BlockOperator, alg::BlockBP) 
	@assert blk.partition == U.partition
	parallel_compute_messages!(blk, alg)

	index = CartesianIndices((nrows(blk), ncols(blk)))

	n_rank = nworkers()
	# (length(index) % n_rank == 0) || println("total number of jobs $(length(index)) can not be divided by number of workers $(n_rank).")
	nave, resi, jobs = partition_jobs(length(index), n_rank)


	f_per_rank(l::Int) = begin
		r = []
		for k in jobs[l]+1:jobs[l+1]
			idx = index[k]
			i, j = idx[1], idx[2]
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock( U, i, j)
			sweep!(x, sU, alg.update_alg)
			push!(r, x.peps)
		end
		return r
	end 

	r = pmap(f_per_rank, 1:n_rank)
	_collect_peps!(blk, vcat(r...))
end

function _collect_peps!(blk::BeliefPEPSBlock, out)
	index = CartesianIndices((nrows(blk), ncols(blk)))
	@assert length(out) == length(index)

	for l in 1:length(out)
		idx = index[l]
		i, j = idx[1], idx[2]
		blk.peps.data[PEPSTools.rowindices(blk, i), PEPSTools.colindices(blk, j)] = out[l]
	end
end

function parallel_compute_messages!(blk::PEPSTools.AbstractBeliefPEPSBlock, alg::PEPSTools.AbstractBlockBPPEPSUpdateAlgorithm)
	iter = 1
	losses = Float64[]
	maxiter = alg.msg_maxiter
	tol = alg.msg_tol
	mult_alg = PEPSTools.get_msg_mult_alg(alg)
	while iter <= maxiter
		row_msgs, col_msgs = parallel_compute_out_message(blk, mult_alg)
		row_err = sum(PEPSTools.message_distance2(a, b) for (a, b) in zip(row_msgs, blk.row_msgs))
		col_err = sum(PEPSTools.message_distance2(a, b) for (a, b) in zip(col_msgs, blk.col_msgs))
		mse_loss = (row_err + col_err) / (length(row_msgs) + length(col_msgs))
		# println("current loss is $(mse_loss)")
		blk.row_msgs[:] = row_msgs
		blk.col_msgs[:] = col_msgs
		push!(losses, mse_loss)
		if iter > 1
			loss_before = losses[end-1]
			err = abs((loss_before - mse_loss) / losses[1])
			if err < tol
				(alg.verbosity > 2) && println("early converge in $iter sweeps with relative error $err")
				break
			end
		end
		iter += 1
	end
	if (alg.verbosity > 1) && (iter > maxiter)
		println("fail to converge to precision $(tol) in $(maxiter) sweeps.")
	end
	return losses
end

function parallel_compute_out_message(blk::PEPSTools.AbstractBeliefPEPSBlock, alg::MPSCompression)
	index = CartesianIndices((nrows(blk), ncols(blk)))
	n_rank = nworkers()
	# (length(index) % n_rank == 0) || println("total number of jobs $(length(index)) can not be divided by number of workers $(n_rank).")
	nave, resi, jobs = partition_jobs(length(index), n_rank)

	f_per_rank(l::Int) = begin
		r = []
		for k in jobs[l]+1:jobs[l+1]
			idx = index[k]
			i, j = idx[1], idx[2]
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			push!(r, PEPSTools.compute_out_messages(x, alg)) 
		end
		return r
	end 

	r = pmap(f_per_rank, 1:n_rank)
	return _collect_out_message(blk, vcat(r...))
end

function _collect_out_message(blk::PEPSTools.AbstractBeliefPEPSBlock, out)
	index = CartesianIndices((nrows(blk), ncols(blk)))
	row_msgs = copy(blk.row_msgs)
	col_msgs = copy(blk.col_msgs)
	@assert length(out) == length(index)
	for l in 1:length(out)
		left, right, up, down = out[l]
		idx = index[l]
		i, j = idx[1], idx[2]
		row_msgs[i, j] = Message(i=row_msgs[i, j].i, o=up)
		row_msgs[i+1, j] = Message(i=down, o=row_msgs[i+1, j].o)
		col_msgs[i, j] = Message(i=col_msgs[i, j].i, o=left)
		col_msgs[i, j+1] = Message(o=col_msgs[i,  j+1].o, i=right)
	end
	return row_msgs, col_msgs
end

function partition_jobs!(n::Int, nt::Int, r::Vector{Int})
    (length(r) == nt+1) || error("wrong thread vector.")
    nave, resi = div(n, nt), n % nt
    r[1] = 0
    for i in 1:resi
        r[i+1] = r[i] + nave + 1
    end
    for i in resi+1:nt
        r[i+1] = r[i] + nave
    end
    return nave, resi
end

function partition_jobs(n::Int, nt::Int)
    r = Vector{Int}(undef, nt+1)
    nave, resi = partition_jobs!(n, nt, r)
    return nave, resi, r
end
