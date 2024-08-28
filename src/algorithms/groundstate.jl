

struct ImaginaryTimePEPS{A<:AbstractPEPSUpdateAlgorithm, B<:AbstractPEPSUpdateAlgorithm} <: AbstractPEPSGroundStateAlgorithm
	stepper::A
	measure_alg::B
	parameters::Vector{Tuple{Int, Float64, Float64}}
	sweeps_per_measure::Int
	verbosity::Int
end


ImaginaryTimePEPS(paras::Vector{<:Tuple{Int, <:Real, <:Real}}; verbosity::Int=1, 
	stepper::AbstractPEPSUpdateAlgorithm=BoundaryMPS(verbosity=verbosity), measure_alg::AbstractPEPSUpdateAlgorithm=stepper, 
	sweeps_per_measure::Int=10) = ImaginaryTimePEPS(
	stepper, measure_alg, convert(Vector{Tuple{Int, Float64, Float64}}, paras), sweeps_per_measure, verbosity)

ImaginaryTimePEPS(paras::Vector{<:Tuple{Int, <:Real}}=[(1000, -0.01), (1000, -0.001), (500, -0.0001), (200, -0.00001)]; 
	tol::Real=1.0e-6, kwargs...) = ImaginaryTimePEPS([(item[1], item[2], tol) for item in paras]; kwargs...)


function ground_state!(peps::AbstractPEPS, h::SquareLatticeHamiltonianBase, alg::ImaginaryTimePEPS)
	sh = squeeze(h)
	ntail = 10
	expec_alg = alg.measure_alg
	verbosity = alg.verbosity
	# always measure initial energy and final energy
	energies = [real(expectation(sh, peps, expec_alg)) / prod(size(peps))]
	(verbosity > 2) && println("initial state energy is $(energies[1])")
	iter = 1
	err = 1.
	for (nsteps, δt, tol) in alg.parameters
		(verbosity > 2) && println("do at most $nsteps sweeps with δt=$(δt) and convergence tolerance $tol")
		(δt < 0.) || error("δt should be negative.")
		U = exponential(sh, δt)
		energies_tmp = Float64[]
		for j in 1:nsteps
			if verbosity > 2
				t = @elapsed sweep!(peps, U, alg.stepper)
				println("the $iter-th sweep takes $t seconds")
			else
				sweep!(peps, U, alg.stepper)
			end
			

			if (alg.sweeps_per_measure > 0) && (iter % alg.sweeps_per_measure == 0)
				expec_peps = expectation(sh, peps, expec_alg) / prod(size(peps))
				if verbosity >= 2
					println("energy after $iter-th sweep is $expec_peps")
				end
				push!(energies_tmp, real(expec_peps))
				if length(energies_tmp) >= ntail
					err = iterative_error_2(energies_tmp[end-ntail+1:end])
					(verbosity > 2) &&println("error at the $iter-th sweep is $err")
					if err < tol
						(verbosity > 1) && println("early converge in $iter sweeps with error $err.")
						break
					end
				end
				(verbosity > 2) && println()
			end

			iter += 1

		end
		append!(energies, energies_tmp)
	end

	return energies, (iter, err)
end

# function ground_state!(peps::PEPS, h::SquareLatticeHamiltonianBase, alg::ImaginaryTimePEPS{PEPSSimpleUpdate})
# 	cpeps = CanonicalPEPS(peps)
# 	res = ground_state!(cpeps, h, alg)
# 	r = PEPS(cmps)
# 	peps.data[:] = r.data[:]
# 	return res
# end
