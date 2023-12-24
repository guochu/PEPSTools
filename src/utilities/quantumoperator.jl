

function QuantumSpins.QuantumOperator(h::SquareLatticeHamiltonianBase)
	index = LinearIndices(size(h))
	terms = []
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			if !isnothing(h.H[i, j])
				if (j == size(h,2))
					for (a, b) in h.H[i, j]
						push!(terms, QTerm(index[i, j]=>a, index[i, 1]=>b))
					end	
				else
					for (a, b) in h.H[i, j]
						push!(terms, QTerm(index[i, j]=>a, index[i, j+1]=>b))
					end						
				end
			end
		end
	end
	for i in 1:size(h,1)
		for j in 1:size(h,2)
			if !isnothing(h.V[i, j])
				if i == size(h,1) 
					for (a, b) in h.V[i, j]
						push!(terms, QTerm(index[i, j]=>a, index[1, j]=>b))
					end					
				else
					for (a, b) in h.V[i, j]
						push!(terms, QTerm(index[i, j]=>a, index[i+1, j]=>b))
					end					
				end
			end
		end
	end
	return QuantumOperator([terms...])
end

