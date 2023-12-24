


struct LocalObservers{M}
	data::PeriodicArray{Union{M, Nothing}, 2}
end

LocalObservers{M}(shape::Tuple{Int, Int}) where {M <: AbstractArray} = LocalObservers{M}(PeriodicArray{Union{M, Nothing}, 2}(nothing, shape))
LocalObservers(::Type{M}, shape::Tuple{Int, Int}) where {M<:AbstractArray} = LocalObservers{M}(shape)
LocalObservers(data::AbstractMatrix{Union{M, Nothing}}) where {M<:AbstractArray} = LocalObservers(PeriodicArray(data))
LocalObservers(data::AbstractMatrix{M}) where {M<:AbstractArray} = LocalObservers(convert(Matrix{Union{M, Nothing}}, data))


Base.size(x::LocalObservers) = size(x.data)
Base.size(x::LocalObservers, i::Int) = size(x.data, i)
Base.repeat(x::LocalObservers, i::Int...) = LocalObservers(repeat(x.data.data, i...))

Base.getindex(x::LocalObservers, i::ValidIndices, j::ValidIndices) = getindex(x.data, i, j)
Base.setindex!(x::LocalObservers, v, i::ValidIndices, j::ValidIndices) = setindex!(x.data, v, i, j)

nontrivial_terms(x::LocalObservers) = nontrivial_terms(x.data)


const LocalClassicalObservers{T} = LocalObservers{Array{T, 4}}
const LocalQuantumObservers{T} = LocalObservers{Array{T, 2}}

MagnetizationTensors(data::AbstractMatrix{<:AbstractArray{T, 4}}) where T = LocalObservers(data)
MagnetizationTensors(::Type{T}, shape::Tuple{Int, Int}) where {T <: Number} = LocalObservers(Array{T, 4}, shape)
function MagnetizationTensors(shape::Tuple{Int, Int}, op::Array{T,4}) where {T <: Number}
	r = PeriodicArray{Union{Array{T,4}, Nothing}, 2}(undef, shape)
	for i in 1:length(r)
		r[i] = op
	end
	return LocalObservers(r)
end
MagnetizationTensors(m::Int, n::Int, args...) = MagnetizationTensors((m, n), args...)

