abstract type Abstract2DTN{T<:Number} end
scalartype(::Type{<:Abstract2DTN{T}}) where {T<:Number} = T
scalartype(x::Abstract2DTN) = scalartype(typeof(x))
raw_data(x::Abstract2DTN) = x.data.data

const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

Base.size(x::Abstract2DTN) = size(x.data)
Base.size(x::Abstract2DTN, i::Int) = size(x.data, i)
Base.getindex(x::Abstract2DTN, i::ValidIndices, j::ValidIndices) = getindex(x.data, i, j)
Base.setindex!(x::Abstract2DTN, v, i::ValidIndices, j::ValidIndices) = setindex!(x.data, v, i, j)
Base.isempty(x::Abstract2DTN) = isempty(x.data)

abstract type AbstractPEPS{T} <: Abstract2DTN{T} end

include("peps.jl")
include("singlelayer.jl")