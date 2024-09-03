

# function fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm, nitr::Int; verbosity::Int=1)
# 	err = 1.
# 	for i in 1:nitr
# 		msg_2 = update_messages(tn, msg)
# 		msg_2 = normalize!(msg_2, alg)
# 		err = average_distance2(msg_2, msg)
# 		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
# 		msg = msg_2
# 	end
# 	(verbosity > 0) && println("final BP error is ", err)
# 	return msg
# end

function fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::MessageNormalizationAlgorithm, nitr::Int, tol::Real; verbosity::Int=1)
	err = 1.
	i = 0
	while i < nitr
		msg_2 = update_messages(tn, msg)
		msg_2 = normalize!(msg_2, alg)
		err = message_distance2(msg_2, msg)
		msg = msg_2
		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
		i += 1
		if (tol > 0) && (err < tol)
			(verbosity > 0) && println("BP converges in $i iterations, error=", err)
			return msg
		end
	end
	if (verbosity > 0)
		if tol > 0
			(i == nitr) && println("BP fails to converge in $(nitr) iterations, error=", err)
		else
			println("final BP error after $(nitr) iterations is ", err)
		end
	end
	return msg
end


function update_messages(tn::SquareTN, msg::SquareLatticeBondMessages)
	msg_new = similar(msg)
	workspace = Vector{scalartype(msg)}()
	for node in 1:length(tn)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		# msg_out = compute_out_messages(aj_data, msg_in)
		msg_out = sl_compute_out_messages(aj_data, msg_in, workspace)
		for (j, n) in enumerate(neighbors(msg, node))
			msg_new[node=>n] = msg_out[j]
		end
	end
	return msg_new
end
function update_messages(tn::PEPS, msg::SquareLatticeBondMessages)
	msg_new = similar(msg)
	for node in 1:length(tn)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		msg_out = dl_compute_out_messages(aj_data, msg_in)
		for (j, n) in enumerate(neighbors(msg, node))
			msg_new[node=>n] = msg_out[j]
		end
	end
	return msg_new
end

get_in_messages(msg::SquareLatticeBondMessages, node::Int) = [msg[n=>node] for n in neighbors(msg, node)]