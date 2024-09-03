"""
    struct PEPS{T<:Number}
"""
struct PEPS{T<:Number} <: AbstractPEPS{T}
	data::PeriodicArray{Array{T, 5}, 2}
end

# PEPS{T}(m::Int, n::Int) where {T<:Number} = PEPS(PeriodicArray{Array{T, 5}, 2}(undef, m, n))
# PEPS(::Type{T}, m::Int, n::Int) where {T<:Number} = PEPS{T}(m, n)
PEPS(data::AbstractMatrix{Array{T, 5}}) where {T <: Number} = PEPS(PeriodicArray(data))

Base.complex(x::PEPS) = PEPS(complex.(x.data))
Base.conj(x::PEPS) = PEPS(conj(x.data))
Base.copy(x::PEPS) = PEPS(copy(x.data))

Base.repeat(x::PEPS, i::Int...) = PEPS(repeat(x, i...))

"""
    physical_dimensions(mps::AbstractPEPS)
Return physical dimensions of mps
"""
physical_dimensions(peps::PEPS) = size.(peps.data, 5)


Flux.@functor PEPS (data,)


# initializers
function prodpeps(::Type{T}, ds::AbstractMatrix{Int}, physectors::AbstractMatrix{Int}) where {T<:Number}
    @assert size(ds) == size(physectors)
    m, n = size(ds)
    r = Matrix{Array{T, 5}}(undef, m, n)
    for j in 1:n
        for i in 1:m
            dj = ds[i, j]
           r[i, j] = reshape(onehot(T, dj, physectors[i, j]), (1, 1, 1, 1, dj))
        end
    end
    return PEPS(r)
end
prodpeps(::Type{T}, physectors::AbstractMatrix{Int}; d::Int) where {T<:Number} = prodpeps(T, ones(Int, size(physectors)) .* d, physectors)
prodpeps(ds::AbstractMatrix{Int}, physectors::AbstractMatrix{Int}) = prodpeps(Float64, ds, physectors)
prodpeps(physectors::AbstractMatrix{Int}; kwargs...) = prodpeps(Float64, physectors; kwargs...)

function randompeps(f, ::Type{T}, ds::AbstractMatrix{Int}; periodic::Bool=false, D::Int) where {T<:Number}
    m, n = size(ds)
    r = Matrix{Array{T, 5}}(undef, m, n)
    for j in 1:n
        for i in 1:m
            s1 = ds[i, j]
            s2 = s3 = s4 = s5 = D
            if !periodic
                s2 = j==1 ? 1 : D
                s3 = i==1 ? 1 : D
                s4 = j==n ? 1 : D
                s5 = i==m ? 1 : D
            end
            tmp = f(T, s2, s3, s4, s5, s1)
            if D==1
                r[i, j] = normalize!(tmp)
            else
                r[i, j] = tmp / D
            end
        end
    end
    return PEPS(r)
end

randompeps(f, T::Type{<:Number}, m::Int, n::Int; d::Int, kwargs...) = randompeps(f, T, ones(Int, m, n) .* d; kwargs...)
randompeps(T::Type{<:Number}, m::Int, n::Int; kwargs...) = randompeps(rand, T, m, n; kwargs...)
randompeps(ds::AbstractMatrix{Int}; kwargs...) = randompeps(Float64, ds; kwargs...)
randompeps(m::Int, n::Int; kwargs...) = randompeps(Float64, m, n; kwargs...)



sandwich_single(a::AbstractArray{<:Number, 5}) = sandwich_single(a, a)
function sandwich_single(a::AbstractArray{<:Number, 5}, b::AbstractArray{<:Number, 5})
    @tensor tmp[2,6,3,7,4,8,5,9] := conj(a[2,3,4,5,1]) * b[6,7,8,9,1] 
    return tie(tmp, (2,2,2,2))
end

function sandwich_single(a::AbstractArray{<:Number, 5}, op::AbstractMatrix, b::AbstractArray{<:Number, 5})
    @tensor tmp[2,6,3,7,4,8,5,9] := conj(a[2,3,4,5,1]) * op[1,10] * b[6,7,8,9,10] 
    return tie(tmp, (2,2,2,2))
end
sandwich_single(a::AbstractArray{<:Number, 5}, op::Nothing, b::AbstractArray{<:Number, 5}) = sandwich_single(a, b)

function sandwich(peps::AbstractMatrix{Array{T, 5}}, op::AbstractDict{Tuple{Int, Int}, <:AbstractMatrix} = Dict{Tuple{Int, Int}, Matrix{Float64}}()) where {T}
    m, n = size(peps)
    r = Matrix{Array{T, 4}}(undef, m, n)
    for i in 1:m
        for j in 1:n
            ts = peps[i, j]
            mj = get(op, (i, j), nothing)
            r[i, j] = sandwich_single(ts, mj, ts)
        end
    end
    return r
end
