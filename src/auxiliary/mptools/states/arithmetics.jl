LinearAlgebra.dot(a::MPS, b::MPS) = _dot(a, b) * (scaling(a) * scaling(b))^length(a)

function _dot(psiA::MPS, psiB::MPS) 
	(length(psiA) == length(psiB)) || throw(ArgumentError("dimension mismatch."))
    hold = l_LL(psiA)
    for i in 1:length(psiA)
        hold = updateleft(hold, psiA[i], psiB[i])
    end
    return tr(hold) 
end

function LinearAlgebra.norm(psi::MPS; iscanonical::Bool=false) 
	r = iscanonical ? norm(psi[1]) : sqrt(real(_dot(psi, psi)))
    return r * scaling(psi)^(length(psi))
end


"""
    distance(a, b)
Square of Euclidean distance between a and b.
"""
distance2(a::MPS, b::MPS) = _distance2(a, b)

"""
    distance(a, b)
Euclidean distance between a and b.
"""
distance(a::MPS, b::MPS) = _distance(a, b)