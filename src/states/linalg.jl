# linalg
LinearAlgebra.axpy!(α, x::PEPS, y::PEPS) = axpy!.(α, x.data, y.data)


function contract(state::SquareTN; kwargs...)
	# compute labels
	labels = Vector{Vector{Int}}(undef, length(state))
	start_index = 0
	for node in 1:length(state)
		nbs = neighbors(state, node)
		idxes = Int[]
		for n in nbs
			if n < node
				pos = findfirst(x->x==node, neighbors(state, n))
				isnothing(pos) && error("node $(n) is not a neighbour of node $(node).")
				push!(idxes, labels[n][pos])
			else
				start_index += 1
				push!(idxes, start_index)
			end
		end
		# start_index += 1
		# push!(idxes, start_index)
		labels[node] = idxes
	end
	return only(ncon(reshape(state.data, length(state)), labels; kwargs...))
end

exact_amplitude(state::PEPS, basis::AbstractVector{Int}; kwargs...) = contract(amplitude_tn(state, basis); kwargs...)