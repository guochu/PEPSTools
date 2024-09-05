include("../src/includes.jl")


using Random
using Serialization, JSON


function test_gs_energy(m, n; hz=1, J=1, D=2)

    h1 = Ising2D((m, n), h=hz, J=J)

    T = Float64

    state = randompeps(T, m, n, d=2, D=D)

    sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)

    @time E1 = energy(h1, state, sampler)

    println("NNQS energy ", E1/length(state))


    h2 = squeeze(ising2D(m, n, hz=hz, J=J))

    @time E2 = energy(h2, state, bp_environments(state, 10, 1.0e-8, verbosity=1))

    println("BP energy ", E2/length(state))

    alg = BoundaryMPS(D1=20, D2=D, als_maxiter=10)

    @time E3 = energy(h2, state, alg)

    println("BMPS energy ", E3/length(state))

end

function test_gs_energy_2(m ,n; hz=1, J=1, D=2)


    T = Float64

    state = randompeps(T, m, n, d=2, D=D)

    h2 = squeeze(ising2D(m, n, hz=hz, J=J))

    @time E2 = expectation(h2, state, bp_environments(state, 10, 1.0e-8, verbosity=1))

    println("BP energy ", E2 )

    alg = BoundaryMPS(D1=20, D2=D, als_maxiter=10)

    @time E3 = expectation(h2, state, alg)

    println("BMPS energy ", E3 )

    bp_alg = BlockBP(block_size=(div(m, 2), div(n, 2)), msg_D=8, update_alg=alg)

    @time E4 = expectation(h2, state, bp_alg)

    println("BlockBP energy ", E4 )


    return E2, E3, E4
end

function test_gs_energy_3(m ,n; hz=1, J=1, D=2)


    T = Float64

    state = randompeps(T, m, n, d=2, D=D, periodic=true)

    h2 = squeeze(ising2D(m, n, hz=hz, J=J, periodic=true))

    @time E2 = expectation(h2, state, bp_environments(state, 10, 1.0e-8, verbosity=1))

    println("BP energy ", E2 )

    alg = BoundaryMPS(D1=20, D2=D, als_maxiter=10)

    bp_alg = BlockBP(block_size=(div(m, 2), div(n, 2)), msg_D=8, update_alg=alg)

    @time E4 = expectation(h2, state, bp_alg)

    println("BlockBP energy ", E4 )


    return E2, E4
end
