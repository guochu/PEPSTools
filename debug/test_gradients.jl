include("../src/includes.jl")


using Random
using Serialization, JSON
using Flux, Flux.Optimise


function test_amplitude(L; hz=1, J=1, D=2)

    T = Float64

    state = randompeps(T, L, L, d=2, D=D)

    sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)

    @time Ψ(state, rand((-1, 1), length(state)))

end

function test_amplitude_grad(L; D=2)

    T = Float64
    state = randompeps(T, L, L, d=2, D=D)

    basis = rand((-1, 1), length(state))
    @time grad = gradient(Ψ, state, basis)
end


function test_amplitudes(L, n; hz=1, J=1, D=2)

    T = Float64

    state = randompeps(T, L, L, d=2, D=D)

    sampler = MetropolisLocal(length(state), n_thermal=100, n_sample_per_chain=500, n_discard=10)

    @time Ψ(state, rand((-1, 1), length(state), n))

end

function test_amplitudes_grad(L, n; D=2)

    T = Float64
    state = randompeps(T, L, L, d=2, D=D)

    basis = rand((-1, 1), length(state), n)
    @time begin
        r, back = Zygote.pullback(Ψ_threaded, state, basis)
        back(r)
    end
end