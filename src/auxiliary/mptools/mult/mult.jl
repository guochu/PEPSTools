include("svdmult.jl")
include("iterativemult.jl")

# wrapper
mult(x::MPO, y::MPS, alg::SVDCompression) = svdmult(x, y, alg.trunc)
mult(x::MPO, y::MPS, alg::IterativeCompression) = iterativemult(x, y, alg)
mpompsmult(x::MPO, y::MPS, alg::MPSCompression) = mult(x, y, alg)

