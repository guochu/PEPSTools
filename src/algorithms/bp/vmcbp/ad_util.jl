# utility AD adjoints
Zygote.@adjoint SquareTN(data::AbstractMatrix{Array{T, 4}}) where {T<:Number} = SquareTN(data), z -> (z,)


Zygote.@adjoint update_messages(tn::SquareTN, msg::SquareLatticeBondMessages) = begin
	msg_new = similar(msg)
	backs = Vector{Any}(undef, length(tn.data))
	workspace = Vector{scalartype(msg)}()
	for node in 1:length(tn.data)
		msg_in = get_in_messages(msg, node)
		aj_data = tn[node]
		out, backs[node] = Zygote.pullback(sl_compute_out_messages, aj_data, msg_in, workspace)
		for (j, n) in enumerate(neighbors(msg, node))
			msg_new[node=>n] = out[j]
		end
	end
	return msg_new, z -> begin
		tn_back = similar(tn.data) 
		msg_back = similar(msg)
		for node in 1:length(tn.data)
			msg_out = [z[node=>n] for n in neighbors(msg, node)]
			tn_back[node], msg_in = backs[node](msg_out)
			for (j, n) in enumerate(neighbors(msg, node))
				msg_back[n=>node] = msg_in[j]
			end
		end
		return tn_back, msg_back
	end	
end

# Zygote.@adjoint update_messages_threaded(tn::SquareTN, msg::SquareLatticeBondMessages) = begin
# 	msg_new = similar(msg)
# 	backs = Vector{Any}(undef, length(tn.data))
# 	f(mnew, bk, mold, tnet, i) = begin
# 		msg_in = get_in_messages(mold, i)
# 		aj_data = tnet[i]
# 		msg_out, bk[i] = Zygote.pullback(sl_compute_out_messages_v1, aj_data, msg_in)
# 		for (j, n) in enumerate(neighbors(mold, i))
# 			mnew[i=>n] = msg_out[j]
# 		end		
# 	end
# 	f2(tnet_back, mback, zback, bk, i) = begin
# 		msg_out = [zback[i=>n] for n in neighbors(zback, i)]
# 		tnet_back[i], msg_in = bk[i](msg_out)
# 		for (j, n) in enumerate(neighbors(zback, i))
# 			mback[n=>i] = msg_in[j]
# 		end
# 	end
# 	fetch.([Threads.@spawn f(msg_new, backs, msg, tn, i) for i in 1:length(tn)])
# 	return msg_new, z -> begin
# 		tn_back = similar(tn.data) 
# 		msg_back = similar(msg)
# 		fetch.([Threads.@spawn f2(tn_back, msg_back, z, backs, i) for i in 1:length(tn)])
# 		return tn_back, msg_back
# 	end	
# end