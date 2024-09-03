abstract type AbstractBPEnvironment end

Base.size(m::AbstractBPEnvironment) = size(m.peps)


include("graphinterface.jl")
include("messageinitializers.jl")
include("normalizemessage.jl")
include("updatemessages.jl")
include("doublelayer.jl")
include("classical.jl")