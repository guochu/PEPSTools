abstract type AbstractBPEnvironment end

Base.size(m::AbstractBPEnvironment) = size(m.peps)


include("graphinterface.jl")
include("computemessages/computemessages.jl")
include("bondmessages.jl")
include("updatemessages.jl")
include("doublelayer.jl")
include("classical.jl")