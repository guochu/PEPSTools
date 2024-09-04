println("------------------------------------")
println("|             BP Utilities         |")
println("------------------------------------")


@testset "BP single layer" begin
	D = 5
	tol = 1.0e-8

	for T in (Float64, ComplexF64)
		all_ts = [rand(T, D), rand(T, D, D), rand(T, D, D, D), rand(T, D, D, D, D)]
		for t in all_ts

			msg_in = [rand(T, D) for i in 1:ndims(t)]
			msg_out1, back1 = Zygote.pullback(sl_compute_out_messages, t, msg_in)
			msg_out2, back2 = Zygote.pullback(sl_compute_out_messages_v1, t, msg_in)
			for (a, b) in zip(msg_out1, msg_out2)
				@test norm(a-b) / norm(a) < tol
			end
			z_back = [rand(T, D) for i in 1:ndims(t)]
			t_back1, msg_in_back1 = back1(z_back)
			t_back2, msg_in_back2 = back2(z_back)
			@test norm(t_back1-t_back2) / norm(t_back1) < tol
			if ndims(t) > 1
				for (a, b) in zip(msg_in_back1, msg_in_back2)
					@test norm(a-b) / norm(a) < tol
				end	
			end

			out1, back1 = Zygote.pullback(sl_contract_node, t, msg_in)
			out2, back2 = Zygote.pullback(sl_contract_node_v2, t, msg_in)
			@test abs(out1 - out2) / abs(out1) < tol
			t_back1, msg_in_back1 = back1(one(out1))
			t_back2, msg_in_back2 = back2(one(out1))
			@test norm(t_back1-t_back2) / norm(t_back1) < tol
			for (a, b) in zip(msg_in_back1, msg_in_back2)
				@test norm(a-b) / norm(a) < tol
			end		

		end
	end
end