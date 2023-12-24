
densevector(peps::PEPS) = is_nonperiodic(peps) ? _non_periodic_peps_todense_impl(peps) : _periodic_peps_todense_impl(peps)

function contract_peps(x::SquareTN)
    r = densevector(_convert_to_PEPS(x))
    @assert length(r)==1
    return only(r)
end

function _tr(a::AbstractArray, i::Int, j::Int)
    (i == j) && error("trace dimensions can not be the same.")
    (i >= 1 && i <= ndims(a)) || error("dim $i out of bound.")
    (j >= 1 && j <= ndims(a)) || error("dim $j out of bound.")
    (size(a, i)==size(a, j)) || error("dim $i and $j size mismatch.")
    (ndims(a) == 2) && return tr(a)
    ranges = Any[1:size(a, l) for l in 1:ndims(a)]
    rest_sizes = [size(a, l) for l in 1:ndims(a) if ((l != i) && (l != j))]
    r = zeros(eltype(a), rest_sizes...)
    for l in 1:size(a, i)
        ranges[i] = l
        ranges[j] = l
        r .+= view(a, ranges...)
    end
    return r
end


function _non_periodic_peps_todense_impl(peps::PEPS)
    isempty(peps) && error("input peps is empty.")
    m, n = size(peps)

    # the first column
    if m == 1
        v = dropdims(peps[1,1], dims=(2,3,5))
    else
        v = dropdims(peps[1,1], dims=(2,3))
        for i in 2:m-1
            # s1, s2, s3, s4, s5 = shape(peps[i, 1])    
            # tmp = reshape(peps[i, 1], s1, s3, s4, s5)
            tmp = dropdims(peps[i, 1], dims=2)
            v = contract(v, tmp, ((ndims(v),), (2,)))
        end
        tmp = dropdims(peps[m, 1], dims=(2,5))
        v = contract(v, tmp, ((ndims(v),), (2,)))
    end
    
    (ndims(v) == 2*m) || error("something wrong.")
    ss = [0 for i in 1:2*m]
    for i in 1:m
        ss[i] = 2*i-1
        ss[m+i] = 2*i
    end
    v = permute(v, ss)
    ss = size(v)
    v = reshape(v, prod(ss[1:m]), ss[m+1:2*m]...)
    (ndims(v) == m+1) || error("something wrong.")
    if n == 1
        return reshape(v, length(v))
    end

    # the middle columns
    for j in 2:n-1
        if m == 1
            tmp = dropdims(peps[1, j], dims=(3, 5) )
            v = contract(v, tmp, ((2,), (2,)))
        else
            tmp = dropdims(peps[1, j], dims=3)
            v = contract(v, tmp, ((2,), (2,)))
            for i in 2:m-1
                v = contract(v, peps[i, j], ((2, ndims(v)), (2, 3)))
            end
            tmp = dropdims(peps[m, j], dims=5)
            v = contract(v, tmp, ((2, ndims(v)), (2, 3)))
        end
        
        (ndims(v) == 2*m+1) || error("something, wrong.")
        ss = [0 for i in 1:2*m+1]
        ss[1] = 1
        for i in 1:m
            ss[i+1] = 2*i
            ss[m+i+1] = 2*i+1
        end
        v = permute(v, ss)
        ss = size(v)
        v = reshape(v, prod(ss[1:m+1]), ss[m+2:end]...)
        (ndims(v) == m+1) || error("something wrong.")
    end

    # the last column
    if m == 1
        tmp = dropdims(peps[1, n], dims=(3,4,5))
        v = contract(v, tmp, ((2,), (2,)))
    else
        tmp = dropdims(peps[1, n], dims=(3,4))
        v = contract(v, tmp, ((2,), (2,)))
        for i in 2:m-1

            tmp = dropdims(peps[i, n], dims=4)
            v = contract(v, tmp, ((2, ndims(v)), (2, 3)))
        end
        tmp = dropdims(peps[m, n], dims=(4,5))
        v = contract(v, tmp, ((2, ndims(v)), (2, 3)))
    end

    return reshape(v, length(v))
end

function _periodic_peps_todense_impl(peps::PEPS)
    isempty(peps) && error("input peps is empty.")
    m, n = size(peps)

    # the first column
    if m == 1
        v = _tr(peps[1, 1], 3, 5)
    else
        v = peps[1,1]
        for i in 2:m-1
            v = contract(v, peps[i, 1], ((ndims(v),), (3,)))
        end
        v = contract(v, peps[m, 1], ((2, ndims(v)), (5, 3)))
    end

     @assert ndims(v) == 3 * m

    if n == 1
        k = 1
        for i in 1:m
            v = _tr(v, k+1, k+2)
            k += 1
        end
        return reshape(v, length(v))
    end

    v = permute(v, vcat(2:3:3*m, 1:3:3*m, 3:3:3*m))
    size_v = size(v)
    v = reshape(v, (size_v[1:m]..., prod(size_v[m+1:2*m]), size_v[2*m+1:3*m]...))

    # the middle columns
    for j in 2:n-1
        @assert ndims(v) == 2*m+1
        if m == 1
            tmp = _tr(peps[1, j], 3, 5)
            v = contract(v, tmp, ((m+2,), (2,)))
        else
            v = contract(v, peps[1, j], ((m+2,), (2,)))
            for i in 2:m-1
                v = contract(v, peps[i, j], ((m+2, ndims(v)), (2, 3)) )
            end
            v = contract(v, peps[m, j], ((m+2, m+4, ndims(v)), (2, 5, 3)))
        end
        
        @assert ndims(v) == 3*m + 1
        v = permute(v, vcat(1:m+1, m+2:2:3*m+1, m+3:2:3*m+1))
        size_v = size(v)
        v = reshape(v, (size_v[1:m]..., prod(size_v[m+1:2*m+1]), size_v[2*m+2:3*m+1]...))
    end

    # the last column
    if m == 1
        tmp = _tr(peps[1, n], 3, 5)
        v = contract(v, tmp, ((1, m+2), (3, 2)))
    else
        v = contract(v, peps[1, n], ((1, m+2), (4, 2)))
        k = m+1
        for i in 2:m-1
            v = contract(v, peps[i, n], ((1,k,ndims(v)), (4,2,3)))
            k -= 1
        end
        v = contract(v, peps[m, n], ((1, k, 5, ndims(v)), (4,2,5,3)))
    end

    return reshape(v, length(v))
end


_convert_to_PEPS(x::SquareTN) = PEPS([_insert_dim(item) for item in x.data])

