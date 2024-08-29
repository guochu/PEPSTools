




sweep!(peps::PEPS, U::SquareLatticeOperator, alg::BlockBP) = sweep!(peps, default_splitting(U, alg.block_size), alg)
function sweep!(peps::PEPS, Us::Vector{<:BlockOperator}, alg::BlockBP)
	for U in Us
		blk = peps_partition(peps, U.partition)
		sweep!(blk, U, alg)
	end
end

function sweep!(blk::BlockBPPartitionPEPS, U::BlockOperator, alg::BlockBP) 
	@assert blk.partition == U.partition
	compute_messages!(blk, alg)
	# println("block rows $(nrows(blk)), cols $(ncols(blk))")
	# println("nontrial terms $(n_nontrivial_terms(U))")
	# println(collect(rowindices(blk, 1)))
	# println(collect(colindices(blk, 1)))
	for i in 1:nrows(blk)
		for j in 1:ncols(blk)
			_peps, msgl, msgr, msgu, msgd = subblock(blk, i, j)
			x = borderedpeps(_peps, left=msgl.i, right=msgr.o, up=msgu.i, down=msgd.o)
			sU = subblock( U, i, j)
			# println(isnothing.(sU.H))
			# println(isnothing.(sU.V))
			sweep!(x, sU, alg.update_alg)

			row_pos, row_idx = rowcenterenumerate(blk.partition, i)
			col_pos, col_idx = colcenterenumerate(blk.partition, j)
			# blk.peps.data[rowindices(blk, i), colindices(blk, j)] = x.peps
			blk.peps.data[row_idx, col_idx] = x.peps[row_pos, col_pos]
			# only_set_necessary!(blk.peps.data, rowindices(blk, i), colindices(blk, j), x.peps, sU)
		end
	end
	# println(size.(blk.peps.data))
end


function only_set_necessary!(target, idx1, idx2, source, sU)
	for (pos_i, i) in enumerate(idx1)
		for (pos_j, j) in enumerate(idx2)
			if !isnothing(sU.H[pos_i, pos_j])
				# println("update H $pos_i, $pos_j")
				# @assert (target[i, j] !== source[pos_i, pos_j])
				# @assert (target[i, j+1] !== source[pos_i, pos_j+1])
				target[i, j], target[i, j+1] = source[pos_i, pos_j], source[pos_i, pos_j+1]
			end
			if !isnothing(sU.V[pos_i, pos_j])
				# println("update V $pos_i, $pos_j")
				# @assert (target[i, j] !== source[pos_i, pos_j])
				# @assert (target[i+1, j] !== source[pos_i+1, pos_j])
				target[i, j], target[i+1, j] = source[pos_i, pos_j], source[pos_i+1, pos_j]
			end

			# if (target[i, j] !== source[pos_i, pos_j])
			# 	println("set $pos_i, $pos_j element")
			# 	target[i, j] = source[pos_i, pos_j]
			# else
			# 	println("neglect $pos_i, $pos_j element")
			# end
		end
	end
end
