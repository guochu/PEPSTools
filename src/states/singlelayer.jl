"""
    struct SquareTN{T<:Number}
"""
struct SquareTN{T<:Number} <: Abstract2DTN{T}
	data::PeriodicArray{Array{T, 4}, 2}
end


SquareTN{T}(m::Int, n::Int) where {T<:Number} = SquareTN{PeriodicArray{Array{T, 4}, 2}}(undef, m, n)
SquareTN(::Type{T}, m::Int, n::Int) where {T<:Number} = SquareTN{T}(m, n)
SquareTN(data::AbstractMatrix{Array{T, 4}}) where {T <: Number} = SquareTN(PeriodicArray(data))


Base.complex(x::SquareTN) = SquareTN(complex.(x.data))
Base.conj(x::SquareTN) = SquareTN(conj(x.data))
Base.copy(x::SquareTN) = SquareTN(copy(x.data))

Base.repeat(x::SquareTN, i::Int...) = SquareTN(repeat(x.data, i...))

SquareTN(x::PEPS) = SquareTN(sandwich(x.data))


"""
    bond_dimensions(peps::SquareTN)
Return bond dimensions of peps

two dimensional tensor network.
The index convention for single site
peps: ------------2-------------
--------------1-------3---------
------------------4-------------
"""
function bond_dimensions(peps::SquareTN)
    m, n = size(peps)
    Vs = Matrix{Int}(undef, m, n) 
    Hs = Matrix{Int}(undef, m, n) 
    for j in 1:n
        for i in 1:m
            Vs[i, j] = size(peps[i, j], 4)
            Hs[i, j] = size(peps[i, j], 3)
        end
    end
    return Dict("H"=>Hs, "V"=>Vs)
end


function is_h_nonperiodic(peps::SquareTN)
    m, n = size(peps)
    for i in 1:m
        ((size(peps[i, 1], 1) == 1) && (size(peps[i, n], 3) == 1)) || return false
    end
    return true
end
function is_v_nonperiodic(peps::SquareTN)
    m, n = size(peps)
    for i in 1:n
        (size(peps[1, i], 2) == 1) && (size(peps[m, i], 4) == 1) || return false
    end
    return true
end
is_nonperiodic(peps::SquareTN) = is_h_nonperiodic(peps) && is_v_nonperiodic(peps)

check(peps::SquareTN; periodic::Bool=false) = periodic ? check_periodic(peps) : check_non_periodic(peps)
function check_non_periodic(peps::SquareTN)
    isempty(peps) && error("peps is empty.")
    m, n = size(peps)
    (size(peps[1,1], 1)==1 && size(peps[1,1], 2)==1) || return false
    for i in 1:n
        (size(peps[1, i], 2)==1) || return false
        (size(peps[m, i], 4)==1) || return false
    end
    for i in 1:m
        (size(peps[i, 1], 1)==1) || return false
        (size(peps[i, n], 3)==1) || return false
    end
    for j in 1:n-1
        for i in 1:m-1
            (size(peps[i, j], 3) == size(peps[i, j+1], 1)) || return false
            (size(peps[i, j], 4) == size(peps[i+1, j], 2)) || return false
        end
    end
    return true
end

function check_periodic(peps::SquareTN)
    isempty(peps) && error("peps is empty.")
    m, n = size(peps)
    for j in 1:n-1
        for i in 1:m-1
            (size(peps[i, j], 3) == size(peps[i, j+1], 1)) || return false
            (size(peps[i, j], 4) == size(peps[i+1, j], 2)) || return false
        end
    end
    for j in 1:n
        (size(peps[m, j], 4) == size(peps[1, j], 2)) || return false
    end
    for i in 1:m
        (size(peps[i, n], 4) == size(peps[i, 1], 1)) || return false
    end
    return true
end


function randomsquaretn(f, ::Type{T}, m::Int, n::Int; periodic::Bool=false, D::Int) where {T<:Number}
    r = SquareTN(T, m, n)
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
            if D==1
                r[i, j] = normalize!(tmp)
            else
                r[i, j] = tmp / D
            end
        end
    end
    return r
end
randomsquaretn(::Type{T}, m::Int, n::Int; kwargs...) where {T<:Number} = randomsquaretn(randn, T, m, n; kwargs...)
randomsquaretn(m::Int, n::Int; kwargs...) = randomsquaretn(Float64, m, n; kwargs...)



