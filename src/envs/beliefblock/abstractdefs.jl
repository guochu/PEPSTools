
abstract type AbstractBeliefPEPSBlock{T} end

nrows(m::AbstractBeliefPEPSBlock) = nrows(m.partition)
ncols(m::AbstractBeliefPEPSBlock) = ncols(m.partition)
rowindices(blk::AbstractBeliefPEPSBlock, i::Int) = rowindices(blk.partition, i)
colindices(blk::AbstractBeliefPEPSBlock, j::Int) = colindices(blk.partition, j)

scalartype(::Type{<:AbstractBeliefPEPSBlock{T}}) where {T<:Number} = T
scalartype(x::AbstractBeliefPEPSBlock) = scalartype(typeof(x))
Base.size(x::AbstractBeliefPEPSBlock) = size(x.peps)
Base.size(x::AbstractBeliefPEPSBlock, i::Int) = size(x.peps, i)


function subblock(m::AbstractBeliefPEPSBlock, i::Int, j::Int) 
	_peps = subblock(m.peps.data, m.partition, i, j)
	# a, b = nrows(m), ncols(m)
	msgl = m.col_msgs[i, j]
	msgr = m.col_msgs[i, j+1]
	msgu = m.row_msgs[i, j]
	msgd = m.row_msgs[i+1, j]
	return _peps, msgl, msgr, msgu, msgd
end


function compute_out_messages(blk::AbstractBeliefPEPSBlock, alg::MPSCompression)
	row_msgs = copy(blk.row_msgs)
	col_msgs = copy(blk.col_msgs)
	# a, b = nrows(blk), ncols(blk)
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			# println("$i-th row, $j-th column")
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			left, right, up, down = compute_out_messages(x, alg)
			row_msgs[i, j] = Message(i=row_msgs[i, j].i, o=up)
			row_msgs[i+1, j] = Message(i=down, o=row_msgs[i+1, j].o)
			col_msgs[i, j] = Message(i=col_msgs[i, j].i, o=left)
			col_msgs[i, j+1] = Message(o=col_msgs[i, j+1].o, i=right)
		end
	end
	return row_msgs, col_msgs
end

function compute_messages!(blk::AbstractBeliefPEPSBlock, mult_alg::MPSCompression; maxiter::Int, tol::Real, verbosity::Int)
	iter = 1
	losses = Float64[]
	err = -1.
	while iter <= maxiter
		row_msgs, col_msgs = compute_out_messages(blk, mult_alg)
		row_err = sum(message_distance2(a, b) for (a, b) in zip(row_msgs, blk.row_msgs))
		col_err = sum(message_distance2(a, b) for (a, b) in zip(col_msgs, blk.col_msgs))
		mse_loss = (row_err + col_err) / (length(row_msgs) + length(col_msgs))
		# println("current loss is $(mse_loss)")
		blk.row_msgs[:] = row_msgs
		blk.col_msgs[:] = col_msgs
		push!(losses, mse_loss)
		if iter > 1
			loss_before = losses[end-1]
			err = abs((loss_before - mse_loss) / losses[1])
			if err < tol
				(verbosity > 2) && println("early converge in $iter sweeps with relative error $err")
				break
			end
		end
		iter += 1
	end
	if (verbosity > 1) && (iter > maxiter)
		println("fail to converge to precision $(tol) (error=$err) in $(maxiter) sweeps.")
	end
	return losses
end


# function compute_out_messages_serial(blk::AbstractBeliefPEPSBlock, alg::MPSCompression)
# 	row_msgs = blk.row_msgs
# 	col_msgs = blk.col_msgs
# 	# a, b = nrows(blk), ncols(blk)
# 	for i in 1:nrows(blk)
# 		for j in 1:ncols(blk)
# 			# println("$i-th row, $j-th column")
# 			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
# 			x = PEPSBlock(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
# 			left, right, up, down = compute_out_messages(x, alg)
# 			row_msgs[i, j] = Message(i=row_msgs[i, j].i, o=up)
# 			row_msgs[i+1, j] = Message(i=down, o=row_msgs[i+1, j].o)
# 			col_msgs[i, j] = Message(i=col_msgs[i, j].i, o=left)
# 			col_msgs[i, j+1] = Message(o=col_msgs[i, j+1].o, i=right)
# 		end
# 	end
# 	return row_msgs, col_msgs
# end

# function compute_messages_serial!(blk::AbstractBeliefPEPSBlock, mult_alg::MPSCompression; maxiter::Int, tol::Real, verbosity::Int)
# 	iter = 1
# 	losses = Float64[]
# 	err = -1.
# 	while iter <= maxiter
# 		old_row_msgs = copy(blk.row_msgs)
# 		old_col_msgs = copy(blk.col_msgs)
# 		row_msgs, col_msgs = compute_out_messages_serial(blk, mult_alg)
# 		row_err = sum(message_distance2(a, b) for (a, b) in zip(row_msgs, old_row_msgs))
# 		col_err = sum(message_distance2(a, b) for (a, b) in zip(col_msgs, old_col_msgs))
# 		mse_loss = (row_err + col_err) / (length(row_msgs) + length(col_msgs))
# 		# println("current loss is $(mse_loss)")
# 		push!(losses, mse_loss)
# 		if iter > 1
# 			loss_before = losses[end-1]
# 			err = abs((loss_before - mse_loss) / losses[1])
# 			if err < tol
# 				(verbosity > 2) && println("early converge in $iter sweeps with relative error $err")
# 				break
# 			end
# 		end
# 		iter += 1
# 	end
# 	if (verbosity > 1) && (iter > maxiter)
# 		println("fail to converge to precision $(tol) (error=$err) in $(maxiter) sweeps.")
# 	end
# 	return losses
# end
