abstract type Classical2DModel end

# MagnetizationTensors(data::AbstractMatrix{<:AbstractArray{T, 4}}) where T = SquareLatticeSites(data)
# MagnetizationTensors(::Type{T}, shape::Tuple{Int, Int}) where {T <: Number} = SquareLatticeSites(Array{T, 4}, shape)
# function MagnetizationTensors(shape::Tuple{Int, Int}, op::Array{T,4}) where {T <: Number}
# 	r = PeriodicArray{Union{Array{T,4}, Nothing}, 2}(undef, shape)
# 	for i in 1:length(r)
# 		r[i] = op
# 	end
# 	return SquareLatticeSites(r)
# end
# MagnetizationTensors(m::Int, n::Int, args...) = MagnetizationTensors((m, n), args...)


include("ising.jl")
include("bmps.jl")
include("blockbp.jl")