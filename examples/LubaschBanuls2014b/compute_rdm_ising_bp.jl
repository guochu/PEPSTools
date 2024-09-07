push!(LOAD_PATH, "../../src")

using Random, LinearAlgebra
using PEPSTools
using Serialization, JSON


gen_peps_path(D2::Int, B::Real, block_size::Tuple{Int, Int}) = "data/Table_vii_bp_blocksize$(block_size[1])_$(block_size[2])_D2$(D2)_B$(B).peps"

fidelity(x::AbstractMatrix, y::AbstractMatrix) = real(tr(sqrt(x) * sqrt(y)))


function compute_rmd2s(D2::Int, B::Real)
	block_size = (7, 7)

	peps_path = gen_peps_path(D2, B, block_size)  
	println("read initial peps from path $peps_path")
	peps = Serialization.deserialize(peps_path)

	D1 = 2*D2^2 + 10
	alg1 = BoundaryMPS(D1=D1, D2=D2)
	alg2 = BlockBP(block_size=block_size, update_alg=alg1)
	alg3 = SimpleUpdate(D2=D2)

	rho1s = rdm2sH(peps, alg1)
	rho2s = similar(rho1s)
	for i in 1:size(rho2s, 1)
		for j in 1:size(rho2s, 2)
			rho2s[i, j] = rdm2H(peps, i, j, alg2)
		end
	end
	rho3s = rdm2sH(peps, alg3; periodic=false)

	rho1s = [reshape(item, (4, 4)) for item in rho1s]
	rho2s = [reshape(item, (4, 4)) for item in rho2s]
	rho3s = [reshape(item, (4, 4)) for item in rho3s]


	f1s = [fidelity(a, b) for (a, b) in zip(rho1s, rho2s)]
	f2s = [fidelity(a, b) for (a, b) in zip(rho1s, rho3s)]

	# println("quantum fidelity on site ($i, $j) is $f1 and $f2")

	results = Dict("f_bp"=>f1s, "f_su"=>f2s)

	result_path = "result/ising_bp_rdm2_D2$(D2)_D1$(D1)_B$(B).json"

	open(result_path, "w") do f
		write(f, JSON.json(results))
	end

	return f1s, f2s
end


