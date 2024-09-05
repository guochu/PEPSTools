include("../src/includes.jl")



function main4()

	# Random.seed!(1234)

	m = 4
	n = 4

	D = 2
	T = ComplexF64

	state = randompeps(T, m, n, d=2, D=D)	
	basis = rand(1:2, length(state))

	msg = random_c_bondmessages(state)
	# msg = unit_cmessage(T, graph, D=D)

	normalize_alg = FixedSum()

	normalize!(msg, normalize_alg)
	# println(msg[1])

	tn = amplitude_tn(state, basis)
	# tn = amplitude_tn(state, ones(Int, length(state)))


	# for i in 1:20
	# 	msg_2 = update_messages(tn, copy(msg))
	# 	msg_2 = normalize!(msg_2)
	# 	println("distance at the $i-th step ", distance2(msg_2, msg))
	# 	msg = msg_2
	# end

	msg = fixedpoint_messages(tn, msg, normalize_alg, 20, -1, verbosity=2)

	println("exact value ", contract(tn))

	return bp_contract(tn, canonicalize(msg))
end

