"""
    struct SquareTN{T<:Number}
"""
struct SquareTN{T<:Number} <: Abstract2DTN{T}
	data::PeriodicArray{Array{T, 4}, 2}
end


# SquareTN{T}(m::Int, n::Int) where {T<:Number} = SquareTN(PeriodicArray{Array{T, 4}, 2}(undef, m, n))
# SquareTN(::Type{T}, m::Int, n::Int) where {T<:Number} = SquareTN{T}(m, n)
SquareTN(data::AbstractMatrix{Array{T, 4}}) where {T <: Number} = SquareTN(PeriodicArray(data))


Base.complex(x::SquareTN) = SquareTN(complex.(x.data))
Base.conj(x::SquareTN) = SquareTN(conj(x.data))
Base.copy(x::SquareTN) = SquareTN(copy(x.data))

Base.repeat(x::SquareTN, i::Int...) = SquareTN(repeat(x.data, i...))

SquareTN(x::PEPS) = SquareTN(sandwich(x.data))


function randomsquaretn(f, ::Type{T}, m::Int, n::Int; periodic::Bool=false, D::Int) where {T<:Number}
    r = Matrix{Array{T, 4}}(undef, m, n)
    for j in 1:n
        for i in 1:m
            s2 = s3 = s4 = s5 = D
            if !periodic
                s2 = j==1 ? 1 : D
                s3 = i==1 ? 1 : D
                s4 = j==n ? 1 : D
                s5 = i==m ? 1 : D
            end
            tmp = f(T, s2, s3, s4, s5)
            r[i, j] = normalize!(tmp)
        end
    end
    return SquareTN(r)
end
randomsquaretn(T::Type{<:Number}, m::Int, n::Int; kwargs...) = randomsquaretn(rand, T, m, n; kwargs...)
randomsquaretn(m::Int, n::Int; kwargs...) = randomsquaretn(Float64, m, n; kwargs...)



