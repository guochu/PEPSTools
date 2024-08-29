# these functions are also used for cyclic peps block so they are not typed
function sl_mpoleft_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,4,8,5,9,2,6] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end
	return get_tn.(x.peps[:, i])
end 
function sl_mporight_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[3,7,2,6,5,9,4,8] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[:, i])
end 
function sl_mpoup_util(x, i::Int)
	get_tn(t::AbstractArray{<:Number, 5}) = begin
		@tensor tmp[2,6,5,9,4,8,3,7] := conj(t[2,3,4,5,1]) * t[6,7,8,9,1]
		return tie(tmp, (2,2,2,2))
	end	
	return get_tn.(x.peps[i, :])
end
function sl_mpodown_util(x, i::Int)
	sandwich_single.(x.peps[i, :])
end

dl_mpoleft_util(x, i::Int) = [permute(item, (2,3,4,1)) for item in x.peps[:, i]]
dl_mporight_util(x, i::Int) = [permute(item, (2,1,4,3)) for item in x.peps[:, i]]
dl_mpoup_util(x, i::Int) = [permute(item, (1,4,3,2)) for item in x.peps[i, :]]
dl_mpodown_util(x, i::Int) = x.peps[i, :]