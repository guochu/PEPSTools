

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

function fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP)
	err = 1.
	i = 0
	nitr = alg.msg_maxiter
	tol = alg.msg_tol
	verbosity = alg.verbosity
	damping = alg.damping
	normalize_alg = alg.normalize_alg
	converged = false
	while i < nitr
		msg_2 = update_messages(tn, msg)
		msg_2 = normalize!(msg_2, normalize_alg)
		msg_2 = mixing!(msg_2, msg, damping)
		err = message_distance2(msg_2, msg)
		msg = msg_2
		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
		i += 1
		if (tol > 0) && (err < tol)
			(verbosity > 0) && println("BP converges in $i iterations, error=", err)
			return msg, true
		end
	end
	if (verbosity > 0)
		if tol > 0
			(i == nitr) && println("BP fails to converge in $(nitr) iterations, error=", err)
		else
			println("final BP error after $(nitr) iterations is ", err)
		end
	end
	return msg, converged
end

Zygote.@adjoint fixedpoint_messages(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP) = begin
	err = 1.
	i = 0
	nitr = alg.msg_maxiter
	tol = alg.msg_tol
	verbosity = alg.verbosity
	damping = alg.damping
	normalize_alg = alg.normalize_alg

	backs = []
	converged = false
	while i < nitr
		msg_2, back = Zygote.pullback(update_messages_normalized, tn, msg, damping, normalize_alg)
		push!(backs, back)

		err = message_distance2(msg_2, msg)
		msg = msg_2
		(verbosity > 1) && println("distance at the $i-th BP iteration ", err)
		i += 1
		if (tol > 0) && (err < tol)
			(verbosity > 0) && println("BP converges in $i iterations, error=", err)
			converged = true
			break
		end
	end

	if (verbosity > 0)
		if tol > 0
			(i == nitr) && println("BP fails to converge in $(nitr) iterations, error=", err)
		else
			println("final BP error after $(nitr) iterations is ", err)
		end
	end
	return (msg, converged), z -> begin
		z_back, z_converged = z
		tn_back = zero.(tn.data)
		for back in reverse(backs)
			tn_back_i, z_back, _, _ = back(z_back)
			axpy!.(true, tn_back_i, tn_back)
		end
		return tn_back, z_back, nothing
	end
end

fixedpoint_messages_n(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP) = fixedpoint_messages_n_util(tn, msg, alg.msg_maxiter, alg.damping, alg.normalize_alg)

function fixedpoint_messages_n_util(tn::Abstract2DTN, msg::SquareLatticeBondMessages, nitr::Int, damping::Real, normalize_alg::MessageNormalizationAlgorithm)
	for i in 1:nitr
		msg = update_messages_normalized(tn, msg, damping, normalize_alg)
	end
	return msg
end
Zygote.@adjoint fixedpoint_messages_n(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP) = begin
	msg_out, back = Zygote.pullback(fixedpoint_messages_n_util, tn, msg, alg.msg_maxiter, alg.damping, alg.normalize_alg)
	return msg_out, z -> begin
		tn_back, msg_back, _, _, _ = back(z)
		return tn_back, msg_back, nothing
	end
end

function fixedpoint_messages_n_v1(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP)
	nitr = alg.msg_maxiter
	damping = alg.damping
	normalize_alg = alg.normalize_alg
	for i in 1:nitr
		msg_2 = update_messages(tn, msg)
		msg_2 = normalize!(msg_2, normalize_alg)
		msg_2 = mixing!(msg_2, msg, damping)
		msg = msg_2
	end
	return msg
end

Zygote.@adjoint fixedpoint_messages_n_v1(tn::Abstract2DTN, msg::SquareLatticeBondMessages, alg::BP) = begin
	nitr = alg.msg_maxiter
	damping = alg.damping
	normalize_alg = alg.normalize_alg
	backs = Vector{Any}(undef, nitr)
	for i in 1:nitr
		msg, backs[i] = Zygote.pullback(update_messages_normalized, tn, msg, damping, normalize_alg)
	end
	return msg, z -> begin
		tn_back = zero.(tn.data)
		for i in nitr:-1:1
			tn_back_i, z, _, _ = backs[i](z)
			axpy!.(true, tn_back_i, tn_back)
		end
		return tn_back, z, nothing
	end
end

function update_messages_normalized(tn::SquareTN, msg::SquareLatticeBondMessages, damping::Real, normalize_alg::MessageNormalizationAlgorithm)
	msg_2 = update_messages(tn, msg)
	msg_2 = normalize(msg_2, normalize_alg)
	msg_2 = mixing(msg_2, msg, damping)
	return msg_2
end

function update_messages(tn::SquareTN, msg::SquareLatticeBondMessages)
	msg_new = similar(msg)
	workspace = Vector{scalartype(msg)}()
	for node in 1:length(tn)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		msg_out = sl_compute_out_messages(aj_data, msg_in, workspace)
		for (j, n) in enumerate(neighbors(msg, node))
			msg_new[node=>n] = msg_out[j]
		end
	end
	return msg_new
end
function update_messages_threaded(tn::SquareTN, msg::SquareLatticeBondMessages)
	msg_new = similar(msg)
	f(mnew, mold, tnet, i) = begin
		msg_out = sl_compute_out_messages_v1(tnet[i], get_in_messages(mold, i))
		for (j, n) in enumerate(neighbors(mold, i))
			mnew[i=>n] = msg_out[j]
		end		
	end
	fetch.([Threads.@spawn f(msg_new, msg, tn, i) for i in 1:length(tn)])
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