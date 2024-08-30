"""
	struct ClassicalSandwichEnv{T, _MPS}
"""
struct ClassicalSandwichEnv{T, _MPS} <: AbstractSquareTNSandwichEnv
	up::_MPS
	middle::Vector{Array{T, 4}}
	down::_MPS
	hstorage::Vector{Array{T, 3}}
end

function ClassicalSandwichEnv(up::MPS, middle::Vector{<:AbstractArray}, down::MPS, left::AbstractArray{<:Number, 3}, right::AbstractArray{<:Number, 3})
	T = scalartype(up)
	middle = convert(Vector{Array{T, 4}}, middle)
	return ClassicalSandwichEnv(up, middle, down, compute_hstorage_right(up, middle, down, left, right))
end 

row_environments(up::MPS, middle::Vector{<:AbstractArray{<:Number, 4}}, down::MPS, left::AbstractArray{<:Number, 3}, 
	right::AbstractArray{<:Number, 3}) = ClassicalSandwichEnv(up, middle, down, left, right)


# expectation values of all the sites
function expectation_sites(x::ClassicalSandwichEnv{T}, obs::Vector) where T
	@assert length(x.middle) == length(obs)
	return [unsafe_expectation_site(x, pos, obs[pos]) for pos in 1:length(obs)]
end 
function unsafe_expectation_site(x::ClassicalSandwichEnv, pos::Int, ob::AbstractArray{<:Number, 4})
	tmp = bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], x.middle[pos])
	@tensor n = tmp[1,2,3] * x.hstorage[pos+1][1,2,3]
	# println("norm of n is $(norm(n))")
	tmp = bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], ob)
	@tensor r = tmp[1,2,3] * x.hstorage[pos+1][1,2,3]
	# println("norm of r is $(norm(r))")

	update_storage_left!(x, pos)
	return r / n
end
function expectation_site(x::ClassicalSandwichEnv, pos::Int, ob::AbstractArray{<:Number, 4})
	left = compute_hleft(x, pos)
	tmp = bm_update_left(left, x.up[pos+1], x.down[pos+1], x.middle[pos])
	@tensor n = tmp[1,2,3] * x.hstorage[pos+1][1,2,3]
	tmp = bm_update_left(left, x.up[pos+1], x.down[pos+1], ob)
	@tensor r = tmp[1,2,3] * x.hstorage[pos+1][1,2,3]

	return r / n
end

# expectation values of a bond
function expectation_bond(x::ClassicalSandwichEnv, pos::Int, ob1::AbstractArray{<:Number, 4}, ob2::AbstractArray{<:Number, 4})
	@assert pos < length(x.middle)
	left = compute_hleft(x, pos)
	tmp = bm_update_left(left, x.up[pos+1], x.down[pos+1], x.middle[pos])
	tmp = bm_update_left(tmp, x.up[pos+2], x.down[pos+2], x.middle[pos+1])	
	@tensor n = tmp[1,2,3] * x.hstorage[pos+2][1,2,3]
	tmp = bm_update_left(left, x.up[pos+1], x.down[pos+1], ob1)	
	tmp = bm_update_left(tmp, x.up[pos+2], x.down[pos+2], ob2)	
	@tensor r = tmp[1,2,3] * x.hstorage[pos+2][1,2,3]

	return r / n
end


function update_storage_left!(x::ClassicalSandwichEnv, pos::Int)
    x.hstorage[pos+1] = normalize!(bm_update_left(x.hstorage[pos], x.up[pos+1], x.down[pos+1], x.middle[pos]))
end

function compute_hleft(x::ClassicalSandwichEnv, pos::Int)
	hleft = x.hstorage[1]
	for i in 1:pos-1
		hleft = normalize!(bm_update_left(hleft, x.up[i+1], x.down[i+1], x.middle[i]))
	end
	return hleft
end
