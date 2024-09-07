push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON


function gen_peps_path(m::Int, D2::Int)
	block_size = div(m, 2)
	return "data/Table2_bp_m$(m)_n$(m)_blocksize$(block_size)_$(block_size)_D2$(D2).peps"
end 

fidelity(x::AbstractMatrix, y::AbstractMatrix) = real(tr(sqrt(x) * sqrt(y)))

trace_distance(x::AbstractMatrix, y::AbstractMatrix) = 0.5*tr(abs.(x - y))

function trivial_operator(m::Int, n::Int)
	iden4 = reshape(one(zeros(4, 4)), (2,2,2,2))
	H = Matrix{Union{Array{Float64, 4}, Nothing}}(nothing, m, n)
	V = Matrix{Union{Array{Float64, 4}, Nothing}}(nothing, m, n)
	for i in 1:m
		for j in 1:n-1
			H[i, j] = iden4
		end
	end

	for i in 1:m-1
		for j in 1:n
			V[i, j] = iden4
		end
	end
	return PEPSTools.SquareLattice(V=V, H=H)
end

function get_su_peps(peps::PEPS, D2::Int, nsweeps::Int=1000)
	m, n = size(peps)
	U = trivial_operator(m, n)

	alg = SimpleUpdate(D2=D2)
	cpeps = CanonicalPEPS(peps)

	for i in 1:nsweeps
		sweep!(cpeps, U, alg)
	end

	println("the middle singular vector $(cpeps.Hbonds[10, 10])")

	return PEPS(cpeps)
end


function compute_rmd2s(m::Int, D2::Int)
	peps_path = gen_peps_path(m, D2)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPO(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=(7,7), update_alg=BoundaryMPS(D1=D1, D2=D2))
	alg3 = SimpleUpdate(D2=D2)


	rho1s = PEPSTools.rdm2sH(peps, alg1)
	rho2s = rdm2s(peps, alg2)

	su_nsweeps = 100
	println("trivial iteration for $su_nsweeps sweeps")
	@time su_peps = get_su_peps(peps, D2, su_nsweeps)
	rho3s = rdm2s(su_peps, alg3)


	rho2s = rho2s["H"]
	rho3s = rho3s["H"]

	println(size(rho1s))
	println(size(rho2s))
	println(size(rho3s))

	rho1s = [reshape(item, (4, 4)) for item in rho1s]
	rho2s = [reshape(item, (4, 4)) for item in rho2s]
	rho3s = [reshape(item, (4, 4)) for item in rho3s]


	f1s = [fidelity(a, b) for (a, b) in zip(rho1s, rho2s)]
	f2s = [fidelity(a, b) for (a, b) in zip(rho1s, rho3s)]

	d1s = [trace_distance(a, b) for (a, b) in zip(rho1s, rho2s)]
	d2s = [trace_distance(a, b) for (a, b) in zip(rho1s, rho3s)]

	# println("quantum fidelity on site ($i, $j) is $f1 and $f2")

	results = Dict("f_bp"=>f1s, "f_su"=>f2s, "d_bp"=>d1s, "d_su"=>d2s)

	result_path = "result/periodic_xxx_rdm2_m_$(m)_D2$(D2)_D1$(D1).json"

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end

	return rho1s, rho2s, rho3s
end


