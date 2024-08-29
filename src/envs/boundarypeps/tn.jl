


struct SquareTNBlock{T, _MPS} <: AbstractSquareTNBlock{T}
	peps::Matrix{Array{T, 4}}
	left::_MPS
	right::_MPS
	up::_MPS
	down::_MPS

function SquareTNBlock(peps::Matrix{Array{T, 4}}, left::MPS, right::MPS, up::MPS, down::MPS) where {T<:Number}
	m, n = size(peps)
	@assert length(left) == length(right) == m
	@assert length(up) == length(down) == n
	for i in 1:m
		@assert size(left[i], 2) == size(peps[i, 1], 1)
		@assert size(right[i], 2) == size(peps[i, n], 3)
	end
	for i in 1:n
		@assert size(up[i], 2) == size(peps[1, i], 2)
		@assert size(down[i], 2) == size(peps[m, i], 4)
	end
	new{T, MPS{T, real(T)}}(peps, left, right, up, down)
end

end


function SquareTNBlock(peps::AbstractMatrix{Array{T, 4}}; left::MPS=trivial_mps(T, size(peps, 1)), right::MPS=trivial_mps(T, size(peps, 1)), 
	up::MPS=trivial_mps(T, size(peps, 2)), down::MPS=trivial_mps(T, size(peps, 2))) where {T <: Number}
	# @assert check(peps)
	return SquareTNBlock(peps, left, right, up, down)
end
SquareTNBlock(peps::SquareTN, args...; kwargs...) = SquareTNBlock(raw_data(peps), args...; kwargs...)
PEPSBlock(peps::AbstractMatrix{Array{T, 4}}, args...; kwargs...) where {T<:Number} = SquareTNBlock(peps, args...; kwargs...)

