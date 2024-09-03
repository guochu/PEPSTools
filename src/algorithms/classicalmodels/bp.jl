function magnetizations(x::Classical2DModel, alg::BP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	env = bp_environments(peps, alg.msg_maxiter, alg.msg_tol, verbosity=alg.verbosity)
	return [_local_expectation(magnetization_tensor(x, i, j, β=β), i, j, peps, env, β=β) for i in 1:size(x, 1), j in 1:size(x, 2)]
end

function magnetization(x::Classical2DModel, i::Int, j::Int, alg::BP; β::Real)
	peps = SquareTN(site_tensors(x, β=β))
	mT = magnetization_tensor(x, i, j, β=β)
	env = bp_environments(peps, alg.msg_maxiter, alg.msg_tol, verbosity=alg.verbosity)
	return _local_expectation(mT, i, j, peps, env)
end

function _local_expectation(mT::AbstractArray{<:Number, 4}, i::Int, j::Int, peps::SquareTN, env::ClassicalBPEnv)
	@assert size(mT) == size(peps[i, j])
	messages = env.messages
	index = LinearIndices(size(peps))
	msgs_in = get_in_messages(messages, index[i, j])
	return sl_contract_node(mT, msgs_in) / sl_contract_node(peps[i, j], msgs_in)
end