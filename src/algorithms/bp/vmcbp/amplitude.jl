# compute amplitudes
# _Ψ_util(state::GraphState, x::AbstractVector{Int}, env::AmplitudeEnvironment=environments(state, x, 10, -1)) = get_amplitude(env)

const NNQS_BP_Iterations = 10

# function NNQS._Ψ(state::GraphState, x::AbstractVector{Int})
# 	env = environments(state, _state_to_index(x), NNQS_BP_Iterations, -1, verbosity=0)
# 	return get_amplitude(env)
# end 

_state_to_index(x::AbstractVector{Int}) = [(item == -1) ? 2 : 1 for item in x]

function _Ψ_threaded(state::PEPS, x::AbstractMatrix{Int}, alg::BP=BP(FixedSum()))
	f(r, s, bj, j) = begin
		r[j] = NNQS._Ψ(s, bj, alg)
	end 
	amps = zeros(scalartype(state), size(x, 2))
	fetch.([Threads.@spawn f(amps, state, x[:, j], j) for j in 1:size(x, 2)])
	return transpose(amps)
end

Zygote.@adjoint _Ψ_threaded(state::PEPS, x::AbstractMatrix{Int}, alg::BP=BP(FixedSum())) = begin
	f(r, bks, s, bj, j) = begin
		r[j], bks[j] = Zygote.pullback(NNQS._Ψ, s, bj, alg)
	end 
	amps = zeros(scalartype(state), size(x, 2))
	backs = Vector{Any}(undef, size(x, 2))
	fetch.([Threads.@spawn f(amps, backs, state, x[:, j], j) for j in 1:size(x, 2)])
	return transpose(amps), z -> begin
		state_back = zero.(state.data)
		tmps = [Threads.@spawn backs[j](z[j]) for j in 1:size(x, 2)]
		for j in 1:size(x, 2)
			state_back_j, _, _ = fetch(tmps[j])
			axpy!.(true, state_back_j, state_back)
		end
		return state_back, nothing
	end
end

NNQS.Ψ(state::PEPS, x::AbstractMatrix{Int}, alg::BP=BP(FixedSum())) = (Threads.nthreads() == 1) ? NNQS._Ψ(state, x, alg) : _Ψ_threaded(state, x, alg)
Zygote.@adjoint NNQS.Ψ(state::PEPS, x::AbstractMatrix{Int}, alg::BP=BP(FixedSum())) = (Threads.nthreads() == 1) ? Zygote.pullback(
																					NNQS._Ψ, state, x, alg) : Zygote.pullback(_Ψ_threaded, state, x, alg)

NNQS._Ψ(state::PEPS, x::AbstractVector{Int}, alg::BP=BP(FixedSum())) = _Ψ_util(state, dropgrad(_state_to_index(x)), 
												dropgrad(_init_c_bondmessages(state, alg.initguess, alg.seed)), dropgrad(alg))

function _Ψ_util(state::PEPS, basis::AbstractVector{Int}, msg_init::SquareLatticeBondMessages, alg::BP)
	# println("here... ", alg.seed)
	tn = amplitude_tn(state, basis)
	msg = fixedpoint_messages(tn, msg_init, alg)
	msg = canonicalize(msg)
	return bp_contract(tn, msg)	
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
			for (j, n) in enumerate(neighbors(msg, node))
				msg_back[n=>node] = msgs_in_back[j]
			end
		end
		return tc_data_back, msg_back
	end
end

function amplitude_tn(state::PEPS, basis::AbstractVector{Int})
	(length(state) == length(basis)) || throw(DimensionMismatch("basis size mismatch"))
	return SquareTN(reshape(map((t, x) -> convert(Array, selectdim(t, ndims(t), x)), get_data(state), basis), size(state)))
end

get_data(state::PEPS) = state.data
Zygote.@adjoint get_data(state::PEPS) = get_data(state), z -> (z,)
