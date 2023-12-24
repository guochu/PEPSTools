abstract type AbstractPEPS{T} end
const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

Base.eltype(x::AbstractPEPS{T}) where {T<:Number} = T
Base.size(x::AbstractPEPS) = size(x.data)
Base.size(x::AbstractPEPS, i::Int) = size(x.data, i)
Base.getindex(x::AbstractPEPS, i::ValidIndices, j::ValidIndices) = getindex(x.data, i, j)
Base.setindex!(x::AbstractPEPS, v, i::ValidIndices, j::ValidIndices) = setindex!(x.data, v, i, j)
Base.isempty(x::AbstractPEPS) = isempty(x.data)

raw_data(x::AbstractPEPS) = x.data.data

struct PEPS{T<:Number} <: AbstractPEPS{T}
	data::PeriodicArray{Array{T, 5}, 2}
end

PEPS{T}(m::Int, n::Int) where {T<:Number} = PEPS(PeriodicArray{Array{T, 5}, 2}(undef, m, n))
PEPS(::Type{T}, m::Int, n::Int) where {T<:Number} = PEPS{T}(m, n)
PEPS(data::AbstractMatrix{Array{T, 5}}) where {T <: Number} = PEPS(PeriodicArray(data))

Base.complex(x::PEPS) = PEPS(complex.(x.data))
Base.conj(x::PEPS) = PEPS(conj(x.data))
Base.copy(x::PEPS) = PEPS(copy(x.data))

Base.repeat(x::PEPS, i::Int...) = PEPS(repeat(raw_data(x), i...))

"""
    bond_dimensions(peps::AbstractPEPS)
Return bond dimensions of peps

two dimensional tensor network.
The index convention for single site
peps: ------------3-------------
--------------2---1---4---------
------------------5-------------
"""
function QuantumSpins.bond_dimensions(peps::PEPS)
    m, n = size(peps)
    Vs = Matrix{Int}(undef, m, n) 
    Hs = Matrix{Int}(undef, m, n) 
    for j in 1:n
        for i in 1:m
            Vs[i, j] = size(peps[i, j], 5)
            Hs[i, j] = size(peps[i, j], 4)
        end
    end
    return Dict("H"=>Hs, "V"=>Vs)
end


"""
    physical_dimensions(mps::AbstractPEPS)
Return physical dimensions of mps
"""
function QuantumSpins.physical_dimensions(peps::PEPS)
    m, n = size(peps)
    r = zeros(Int, m, n)
    for j in 1:n
        for i in 1:m
            r[i, j] = size(peps[i, j], 1)
        end
    end
    return r
end

"""
    bond_dimension(mps::AbstractPEPS)
Return maximum bond dimension of mps
"""
function QuantumSpins.bond_dimension(peps::PEPS)
	D = 0
	for item in peps.data
		D = max(D, size(item)...)
	end
	return D
end 

function is_h_nonperiodic(peps::PEPS)
    m, n = size(peps)
    for i in 1:m
        ((size(peps[i, 1], 2) == 1) && (size(peps[i, n], 4) == 1)) || return false
    end
    return true
end
function is_v_nonperiodic(peps::PEPS)
    m, n = size(peps)
    for i in 1:n
        (size(peps[1, i], 3) == 1) && (size(peps[m, i], 5) == 1) || return false
    end
    return true
end
is_nonperiodic(peps::PEPS) = is_h_nonperiodic(peps) && is_v_nonperiodic(peps)

check(peps::PEPS; periodic::Bool=false) = periodic ? check_periodic(peps) : check_non_periodic(peps)
function check_non_periodic(peps::PEPS)
    isempty(peps.data) && error("peps is empty.")
    m, n = size(peps)
    (size(peps[1,1], 2)==1 && size(peps[1,1], 3)==1) || return false
    for i in 1:n
        (size(peps[1, i], 3)==1) || return false
        (size(peps[m, i], 5)==1) || return false
    end
    for i in 1:m
        (size(peps[i, 1], 2)==1) || return false
        (size(peps[i, n], 4)==1) || return false
    end
    for j in 1:n-1
        for i in 1:m-1
            (size(peps[i, j], 4) == size(peps[i, j+1], 2)) || return false
            (size(peps[i, j], 5) == size(peps[i+1, j], 3)) || return false
        end
    end
    return true
end

function check_periodic(peps::PEPS)
    isempty(peps.data) && error("peps is empty.")
    m, n = size(peps)
    for j in 1:n-1
        for i in 1:m-1
            (size(peps[i, j], 4) == size(peps[i, j+1], 2)) || return false
            (size(peps[i, j], 5) == size(peps[i+1, j], 3)) || return false
        end
    end
    for j in 1:n
        (size(peps[m, j], 5) == size(peps[1, j], 3)) || return false
    end
    for i in 1:m
        (size(peps[i, n], 4) == size(peps[i, 1], 2)) || return false
    end
    return true
end


function prodpeps(::Type{T}, ds::AbstractMatrix{Int}, physectors::AbstractMatrix{Int}) where {T<:Number}
    @assert size(ds) == size(physectors)
    m, n = size(ds)
    r = PEPS(T, m, n)
    for j in 1:n
        for i in 1:m
            dj = ds[i, j]
           r[i, j] = reshape(QuantumSpins.onehot(T, dj, physectors[i, j]), (dj, 1, 1, 1, 1))
        end
    end
    return r
end
prodpeps(::Type{T}, physectors::AbstractMatrix{Int}; d::Int) where {T<:Number} = prodpeps(T, ones(Int, size(physectors)) .* d, physectors)
prodpeps(ds::AbstractMatrix{Int}, physectors::AbstractMatrix{Int}) = prodpeps(Float64, ds, physectors)
prodpeps(physectors::AbstractMatrix{Int}; kwargs...) = prodpeps(Float64, physectors; kwargs...)

function randompeps(::Type{T}, ds::AbstractMatrix{Int}; periodic::Bool=false, D::Int) where {T<:Number}
    m, n = size(ds)
    r = PEPS(T, m, n)
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
            tmp = randn(T, s1, s2, s3, s4, s5)
            if D==1
                r[i, j] = normalize!(tmp)
            else
                r[i, j] = tmp / D
            end
        end
    end
    return r
end

randompeps(::Type{T}, m::Int, n::Int; d::Int, kwargs...) where {T<:Number} = randompeps(T, ones(Int, m, n) .* d; kwargs...)
randompeps(ds::AbstractMatrix{Int}; kwargs...) = randompeps(Float64, ds; kwargs...)
randompeps(m::Int, n::Int; kwargs...) = randompeps(Float64, m, n; kwargs...)



sandwich_single(a::AbstractArray{<:Number, 5}) = sandwich_single(a, a)
function sandwich_single(a::AbstractArray{<:Number, 5}, b::AbstractArray{<:Number, 5})
    @tensor tmp[2,6,3,7,4,8,5,9] := conj(a[1,2,3,4,5]) * b[1,6,7,8,9] 
    return tie(tmp, (2,2,2,2))
end

function sandwich_single(a::AbstractArray{<:Number, 5}, op::AbstractMatrix, b::AbstractArray{<:Number, 5})
    @tensor tmp[2,6,3,7,4,8,5,9] := conj(a[1,2,3,4,5]) * op[1,10] * b[10,6,7,8,9] 
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