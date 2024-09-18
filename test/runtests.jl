include("../src/includes.jl")

using Test, Zygote


include("util.jl")

# # auxiliary
# include("tensorops.jl")
# include("bputility.jl")
# include("mptools.jl")

# # state definitions
# include("definitions.jl")

# amplitudes
include("amplitudes.jl")

# # algorithms
# include("algorithms/bpfixedpoint.jl")
# include("algorithms/classic_ising.jl")
# include("algorithms/blockbp.jl")
# include("algorithms/expectation.jl")