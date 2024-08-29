

struct SquareLatticeSites{M}
	data::PeriodicArray{Union{M, Nothing}, 2}
end
SquareLatticeSites{M}(shape::Tuple{Int, Int}) where {M <: AbstractArray} = SquareLatticeSites{M}(PeriodicArray{Union{M, Nothing}, 2}(nothing, shape))
SquareLatticeSites(data::AbstractMatrix{Union{M, Nothing}}) where {M<:AbstractArray} = SquareLatticeSites(PeriodicArray(data))
SquareLatticeSites(data::AbstractMatrix{M}) where {M<:AbstractArray} = SquareLatticeSites(convert(Matrix{Union{M, Nothing}}, data))

Base.size(x::SquareLatticeSites) = size(x.data)
Base.size(x::SquareLatticeSites, i::Int) = size(x.data, i)
Base.repeat(x::SquareLatticeSites, i::Int...) = SquareLatticeSites(repeat(x.data.data, i...))

Base.getindex(x::SquareLatticeSites, i::ValidIndices, j::ValidIndices) = getindex(x.data, i, j)
Base.setindex!(x::SquareLatticeSites, v, i::ValidIndices, j::ValidIndices) = setindex!(x.data, v, i, j)

Base.fill!(x::SquareLatticeSites{M}, v::Union{M, Nothing}) where M = fill!(x.data, v)

n_nontrivial_terms(x::SquareLatticeSites) = n_nontrivial_terms(x.data)


const LocalCObservers{T} = SquareLatticeSites{Array{T, 4}}
const LocalQObservers{T} = SquareLatticeSites{Array{T, 2}}
