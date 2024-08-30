println("------------------------------------")
println("|      classical Ising model       |")
println("------------------------------------")

function _convert_upper(mpsu::Vector)
    n = length(mpsu)
    r = Vector{Any}(undef, n)
    for i in 1:n
        r[i] = dropdims(mpsu[i], dims=2)
        r[i] = permute(r[i], (1,3,2))
    end
    return r
end

function _convert_lower(mpsd::Vector)
    n = length(mpsd)
    r = Vector{Any}(undef, n)
    for i in 1:n
        # r[i], tmp = fusion(mpsd[i], nothing, ((1,4),(1,1)))
        r[i] = dropdims(mpsd[i], dims=4)
    end
    return r
end

function _brute_force_update(r::AbstractArray, m::Vector)
    L = length(m)
    (ndims(r) == L) || error("wrong input tensor rank.")
    tmp = dropdims(m[L], dims=3)
    s = contract(tmp, r, ((2,), (L,)))
    for i in L-1:-1:2
        s = contract(m[i], s, ((2,3), (L+1, 1)))
    end
    tmp = dropdims(m[1], dims=1)
    s = contract(tmp, s, ((1,2), (L+1, 1)))
    return s
end

function contract_2D_nonperiodic(ms::AbstractMatrix)
    m, n = size(ms)
    (m <= 1) && error("the size of PEPS is less than 1.")
    mpsu = _convert_upper(ms[1, :])
    r = dropdims(mpsu[1], dims=1)
    for i in 2:n
        r = contract(r, mpsu[i], ((i,), (1,)))
    end
    r = dropdims(r, dims=n+1)
    # r is an n-dimensional tensor
    for i in 2:m-1
        r = _brute_force_update(r, ms[i, :])
    end

    mpsd = ms[m, :]
    tmp = dropdims(mpsd[n], dims=(3,4))
    r = contract(r, tmp, ((n,), (2,)))
    for i in n-1:-1:2
        tmp = dropdims(mpsd[i], dims=4)
        r = contract(r, tmp, ((ndims(r)-1, ndims(r)), (2,3)))
    end
    tmp = dropdims(mpsd[1], dims=(1,4))
    return vdot(r, tmp)
end

vdot(x::AbstractArray, y::AbstractArray) = dot(conj(x), y)

function exact_magnetization(model::ClassicalIsing2D, i::Int, j::Int; β::Real)
	ms = site_tensors(model, β=β)
	mj = magnetization_tensor(model, i, j, β=β)
	ms_2 = copy(ms)
	ms_2[i, j] = mj
	return contract_2D_nonperiodic(ms_2) / contract_2D_nonperiodic(ms)
end


function exact_bond_energy(model::ClassicalIsing2D, i::Int, j::Int; β::Real)
	ms = site_tensors(model, β=β)
	mj1 = magnetization_tensor(model, i, j, β=β)
	mj2 = magnetization_tensor(model, i, j+1, β=β)
	ms_2 = copy(ms)
	ms_2[i, j] = mj1
	ms_2[i, j+1] = mj2
	return contract_2D_nonperiodic(ms_2) / contract_2D_nonperiodic(ms)
end

function Onsager_m(beta)
	
	beta_C =  log(1+sqrt(2.0))/2
	
	if beta<= beta_C
		exact=0.0
	else
		exact = ( 1-1/(sinh(2*beta))^4)^(1.0/8)
	end
	return exact

end

@testset "Classical Ising model, OBC, BMPS against exact" begin
	m = 4

	D2 = 2
	D1 = 20

	update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=20, mult_alg=IterativeCompression(D=D1,maxiter=20))


	J = 1.
	h = 0.23
	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=false)
	i = j = div(m, 2)

	beta_min = 0.2
	beta_max = 1.5
	beta_n = 40
	del_beta = (beta_max-beta_min)/(beta_n-1)

	betas_list = [beta_min + i*del_beta for i in 0:beta_n-1]

	tol = 1.0e-4
	for i in (2,3)
		for j in (1,2)
			for beta in betas_list
				obs1 = magnetization(model, i, j, update_alg, β=beta)
				obs2 = exact_magnetization(model, i, j, β=beta)

				@test abs(obs1 - obs2) < tol

				obs1 = bond_energy(model, i, j, update_alg, β=beta)
				obs2 = exact_bond_energy(model, i, j, β=beta)

				@test abs(obs1 - obs2) < tol
			end	
		end
	end
end


@testset "Classical Ising model, PBC, BlockBP against analytic" begin
	
	m = 5
	i = j = div(m+1, 2)

	J = 1.0

	D2 = 2
	# D1 = 2*D2^2 + 10
	D1 = 8
	update_alg = BoundaryMPS(D2=D2,D1=D1,als_maxiter=50,mult_alg=IterativeCompression(D=D1,maxiter=50))
	exp_alg = BlockBP(msg_maxiter=50, msg_tol=1.0e-11, msg_D=6, block_size=(m, m), update_alg=update_alg)

	h = 0.
	model = ClassicalIsing2D(m, m, J=J, hz=h, periodic=true)


	betas_list = vcat(collect(0.2:0.05:0.35), collect(0.49:0.05:1))

	tol = 1.0e-2

	for beta in betas_list
		obs1 = abs(magnetization(model, i, j, exp_alg, β=beta))
		obs2 = Onsager_m(beta)

		@test abs(obs1 - obs2) < tol
		# println("beta=$beta, BP=$(obs1), Onsager=$(obs2)")
	end
end


