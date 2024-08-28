include("../src/includes.jl")

using Test


include("util.jl")

# auxiliary
include("tensorops.jl")
include("mptools.jl")

# state definitions
include("definitions.jl")

include("algorithms/classic_ising.jl")
include("algorithms/blockbp.jl")
include("algorithms/expectation.jl")