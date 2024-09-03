# compute amplitudes
# _Ψ_util(state::GraphState, x::AbstractVector{Int}, env::AmplitudeEnvironment=environments(state, x, 10, -1)) = get_amplitude(env)

const NNQS_BP_Iterations = 5

# function NNQS._Ψ(state::GraphState, x::AbstractVector{Int})
# 	env = environments(state, _state_to_index(x), NNQS_BP_Iterations, -1, verbosity=0)
# 	return get_amplitude(env)
# end 

_state_to_index(x::AbstractVector{Int}) = [(item == -1) ? 2 : 1 for item in x]


NNQS._Ψ(state::PEPS, x::AbstractVector{Int}) = _Ψ_util(state, dropgrad(_state_to_index(x)), dropgrad(unit_cmessage(state)))

function _Ψ_util(state::PEPS, basis::AbstractVector{Int}, msg::SquareLatticeBondMessages)
	tn = amplitude_tn(state, basis)
	for i in 1:NNQS_BP_Iterations
		msg_2 = update_messages(tn, msg)
		msg = normalize(msg_2, FixedSum())
	end
	msg = canonicalize(msg)
	return bp_contract(tn ,msg)	
end

bp_contract(tn::SquareTN, msg::SquareLatticeBondMessages) = prod(bp_contract_util(tn, msg))
bp_contract_util(tn::SquareTN, msg::SquareLatticeBondMessages) = [sl_contract_node(tn[j], get_in_messages(msg, j)) for j in 1:length(tn)]



Zygote.@adjoint bp_contract_util(tn::SquareTN, msg::SquareLatticeBondMessages) = begin
	vals = zeros(scalartype(msg), length(tn))
	backs = Vector{Any}(undef, length(tn))
	for node in 1:length(tn)
		msgs_in = get_in_messages(msg, node)
		vals[node], backs[node] = Zygote.pullback(sl_contract_node, tn[node], msgs_in)
	end
	return vals, z -> begin
		tc_data_back = similar(tn.data)
		msg_back = similar(msg)
		for node in 1:length(tn)
			tc_data_back[node], msgs_in_back = backs[node](z[node])
			for (j, n) in enumerate(neighbors(msg.graph, node))
				msg_back[n=>node] = msgs_in_back[j]
			end
		end
		return tc_data_back, msg_back
	end
end

function amplitude_tn(state::PEPS, basis::AbstractVector{Int})
	(length(state) == length(basis)) || throw(DimensionMismatch("basis size mismatch"))
	return SquareTN(map((t, x) -> selectdim(t, ndims(t), x), state.data, basis))
end
