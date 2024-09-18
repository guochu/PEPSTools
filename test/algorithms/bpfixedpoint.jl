println("------------------------------------")
println("|          BP fixedpoint           |")
println("------------------------------------")


@testset "BP fixedpoint adjoint" begin
	D = 2
	
	tol = 1.0e-8
	for T in (Float64, ComplexF64)
		for damping in (0, 0.5)
			bp_alg = BP(FixedSum(), msg_maxiter=5, msg_tol=-1, damping=damping, verbosity=0)

			tn = randomsquaretn(T, 3, 3, D=D)
			msg = random_c_bondmessages(tn)
			msg_out1, back1 = Zygote.pullback(fixedpoint_messages_n, tn, msg, bp_alg)
			msg_out2, back2 = Zygote.pullback(fixedpoint_messages_n_v1, tn, msg, bp_alg)
			@test message_distance2(msg_out1, msg_out2) < tol
			tn_back1, msg_back1, _ = back1(msg_out1)
			tn_back2, msg_back2, _ = back2(msg_out2)
			@test message_distance2(msg_back1, msg_back2) < tol
			for (a, b) in zip(tn_back1, tn_back2)
				@test norm(a-b) / norm(a) < tol
			end

			msg_out3_all, back3 = Zygote.pullback(fixedpoint_messages, tn, msg, bp_alg)
			msg_out3, converged = msg_out3_all
			tn_back3, msg_back3, _ = back3(msg_out3_all)
			@test message_distance2(msg_out1, msg_out3) < tol
			@test message_distance2(msg_back1, msg_back3) < tol
			for (a, b) in zip(tn_back1, tn_back3)
				@test norm(a-b) / norm(a) < tol
			end		

		end
	end
	
end