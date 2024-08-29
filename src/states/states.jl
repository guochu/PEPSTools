abstract type Abstract2DTN{T<:Number} end
scalartype(::Type{<:Abstract2DTN{T}}) where {T<:Number} = T
scalartype(x::Abstract2DTN) = scalartype(typeof(x))
raw_data(x::Abstract2DTN) = x.data.data

# const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

Base.size(x::Abstract2DTN) = size(x.data)
Base.size(x::Abstract2DTN, i::Int) = size(x.data, i)
Base.length(x::Abstract2DTN) = length(x.data)
Base.getindex(x::Abstract2DTN, i::ValidIndices, j::ValidIndices) = getindex(x.data, i, j)
Base.setindex!(x::Abstract2DTN, v, i::ValidIndices, j::ValidIndices) = setindex!(x.data, v, i, j)
Base.isempty(x::Abstract2DTN) = isempty(x.data)



"""
    bond_dimensions(peps::Abstract2DTN)
Return bond dimensions of peps

two dimensional tensor network.
The index convention for single site
peps: ------------3-------------
--------------2---1---4---------
------------------5-------------
"""
function bond_dimensions(peps::Abstract2DTN)
    m, n = size(peps)
    Vs = Matrix{Int}(undef, m, n) 
    Hs = Matrix{Int}(undef, m, n) 
    for j in 1:n
        for i in 1:m
            Vs[i, j] = size(peps[i, j], 4)
            Hs[i, j] = size(peps[i, j], 3)
        end
    end
    return SquareLatticeBonds(H=Hs, V=Vs)
end

"""
    bond_dimension(mps::Abstract2DTN)
Return maximum bond dimension of mps
"""
function bond_dimension(peps::Abstract2DTN)
	D = 0
	for item in peps.data
		D = max(D, maximum(size(item)[1:4]))
	end
	return D
end 

function is_h_nonperiodic(peps::Abstract2DTN)
    m, n = size(peps)
    for i in 1:m
        ((size(peps[i, 1], 1) == 1) && (size(peps[i, n], 3) == 1)) || return false
    end
    return true
end
function is_v_nonperiodic(peps::Abstract2DTN)
    m, n = size(peps)
    for i in 1:n
        (size(peps[1, i], 2) == 1) && (size(peps[m, i], 4) == 1) || return false
    end
    return true
end
is_nonperiodic(peps::Abstract2DTN) = is_h_nonperiodic(peps) && is_v_nonperiodic(peps)
is_periodic(peps::Abstract2DTN) = (!is_h_nonperiodic(peps)) && (!is_v_nonperiodic(peps))



check_consistency(peps::Abstract2DTN; periodic::Bool=false) = periodic ? _check_periodic(peps) : _check_non_periodic(peps)
function _check_non_periodic(peps::Abstract2DTN)
    isempty(peps) && error("peps is empty")
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

function _check_periodic(peps::Abstract2DTN)
    isempty(peps) && error("peps is empty")
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


abstract type AbstractPEPS{T} <: Abstract2DTN{T} end

include("peps.jl")
include("singlelayer.jl")