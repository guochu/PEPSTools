abstract type Classical2DModel end



struct ClassicalIsing2D <: Classical2DModel
	hz::Matrix{Float64}
	Jh::PeriodicArray{Float64, 2}
	Jv::PeriodicArray{Float64, 2}
	isperiodic::Bool
end

Base.size(x::ClassicalIsing2D, args...) = size(x.hz, args...)
isperiodic(x::ClassicalIsing2D) = x.isperiodic

function ClassicalIsing2D(hz::AbstractMatrix{<:Real}, Jh::AbstractMatrix{<:Real}, Jv::AbstractMatrix{<:Real}; periodic::Bool=true)
	m, n = size(hz)
	if periodic
		@assert size(hz) == size(Jh) == size(Jv)
	else
		@assert m == size(Jh, 1) == size(Jv, 1) + 1
		@assert n == size(Jh, 2) + 1  == size(Jv, 2)
	end
	return ClassicalIsing2D(convert(Matrix{Float64}, hz), PeriodicArray(convert(Matrix{Float64}, Jh)), PeriodicArray(convert(Matrix{Float64}, Jv)), periodic)
end

function ClassicalIsing2D(hz::AbstractMatrix{<:Real}; J::Real, periodic::Bool=true)
	if periodic
		r = ClassicalIsing2D(hz, J .* ones(size(hz)), J .* ones(size(hz)); periodic=periodic)
	else
		m, n = size(hz)
		r = ClassicalIsing2D(hz, J .* ones(m, n-1), J .* ones(m-1, n); periodic=periodic)
	end
	return r
end

ClassicalIsing2D(shape::Tuple{Int, Int}; J::Real=1, hz::Real=1.0e-3, periodic::Bool=true) = ClassicalIsing2D(hz .* ones(shape); J=J, periodic=periodic)


ClassicalIsing2D(m::Int, n::Int; kwargs...) = ClassicalIsing2D((m, n); kwargs...)

function site_tensor_util(x::ClassicalIsing2D, i::Int, j::Int, β::Real)
	Jl, Jr, Ju, Jd = x.Jh[i, j-1], x.Jh[i, j], x.Jv[i-1, j], x.Jv[i, j]
	hz = x.hz[i, j]
	T, T0, T1 = ising_tensor_2D(Jl, Jr, Ju, Jd, hz, β, size(x,1), size(x,2), i, j, periodic=isperiodic(x))
	# return T
	return T, T0, T1
end
site_tensor(x::ClassicalIsing2D, i::Int, j::Int; β::Real) = site_tensor_util(x, i, j, β)[1]
function site_tensors(x::ClassicalIsing2D; β::Real)
	r = Matrix{Array{Float64, 4}}(undef, size(x))
	for i in 1:size(r, 1)
		for j in 1:size(r, 2)
			r[i, j] = site_tensor(x, i, j, β=β)
		end
	end
	return r
end
function magnetization_tensor(x::ClassicalIsing2D, i::Int, j::Int; β::Real)
	T, T0, T1 = site_tensor_util(x, i, j, β)
	return T0 - T1
end

function magnetization_tensors(x::ClassicalIsing2D; β::Real)
	r = Matrix{Array{Float64, 4}}(undef, size(x))
	for i in 1:size(r, 1)
		for j in 1:size(r, 2)
			r[i, j] = magnetization_tensor(x, i, j, β=β)
		end
	end
	return r
end

function ising_tensor_2D(Jl::Real, Jr::Real, Ju::Real, Jd::Real,  h::Real, beta::Real, Nx::Int, Ny::Int, i::Int, j::Int; periodic::Bool=false)
	DL = 2
	DR = 2
	DU = 2
	DD = 2
	if !periodic
		if j == 1
			DL = 1
		end
		if j == Ny
			DR = 1
		end
		if i == 1
			DU = 1
		end
		if i == Nx
			DD = 1
		end
	end
	# println("here i=$i, j=$j, DL=$DL, DU=$DU, DR=$DR, DD=$DD")
	T0 = zeros(DL, DU, DR, DD)
	T1 = zeros(DL, DU, DR, DD)

	T0[1,1,1,1] = exp(beta*h)
	T1[end, end, end, end] = exp(-beta*h)

	# x=exp(beta*J)
	# int_mat = [x 1/x; 1/x x]

	# r = cholesky(int_mat)
	# C = Matrix(r.L)

	if DL > 1
		int_mat = gen_bond_matrix(Jl, beta)
		r = cholesky(int_mat)
		C = Matrix(r.L)
		T0 = @tensor tmp[1,3,4,5] := C[2,1] * T0[2,3,4,5]
		T1 = @tensor tmp[1,3,4,5] := C[2,1] * T1[2,3,4,5]
	end
	if DU > 1
		int_mat = gen_bond_matrix(Ju, beta)
		r = cholesky(int_mat)
		C = Matrix(r.L)
		T0 = @tensor tmp[3,1,4,5] := C[2,1] * T0[3,2,4,5]
		T1 = @tensor tmp[3,1,4,5] := C[2,1] * T1[3,2,4,5]
	end
	if DR > 1
		int_mat = gen_bond_matrix(Jr, beta)
		r = cholesky(int_mat)
		C = Matrix(r.L)
		T0 = @tensor tmp[3,4,1,5] := C[2,1] * T0[3,4,2,5]
		T1 = @tensor tmp[3,4,1,5] := C[2,1] * T1[3,4,2,5]
	end
	if DD > 1
		int_mat = gen_bond_matrix(Jd, beta)
		r = cholesky(int_mat)
		C = Matrix(r.L)
		T0 = @tensor tmp[3,4,5,1] := C[2,1] * T0[3,4,5,2]
		T1 = @tensor tmp[3,4,5,1] := C[2,1] * T1[3,4,5,2]
	end

	T = T0 + T1
	nr = norm(T)

	return T/nr, T0/nr, T1/nr
end

function gen_bond_matrix(J::Real, beta::Real)
	x = exp(beta*J)
	return [x 1/x; 1/x x]
end

function decompose_bond(m::AbstractMatrix)
	r = cholesky(m)
	C = Matrix(r.L)
	return C, C'
end
