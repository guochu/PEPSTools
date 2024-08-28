function Base.:*(h::MPO, psi::MPS)
    @assert !isempty(h)
    (length(h) == length(psi)) || throw(DimensionMismatch("mpo mps size mismatch"))
    r = [@tensor tmp[-1 -2; -3 -4 -5] := a[-1, -3, -4, 1] * b[-2, 1, -5] for (a, b) in zip(h.data, psi.data)]
    return MPS([tie(item,(2,1,2)) for item in r], scaling=scaling(h) * scaling(psi))
end