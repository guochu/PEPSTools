function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM)
end

function heisenberg2D(m::Int, n::Int; J::Real=1, Jzz::Real=1, periodic::Bool=false)
	p = spin_half_matrices()
	h = SquareLatticeHamiltonian(Float64, m, n)
	L = periodic ? size(h,2) : size(h,2)-1
	for i in 1:size(h,1)
		for j in 1:L
			push!(h.H[i, j], (0.5 * J * p["+"], p["-"]))
			push!(h.H[i, j], (0.5 * J * p["-"], p["+"]))
			push!(h.H[i, j], (0.25 * Jzz * p["z"], p["z"]))
		end
	end

	L = periodic ? size(h,1) : size(h,1) - 1
	for i in 1:L
		for j in 1:size(h,2)
			push!(h.V[i, j], (0.5 * J * p["+"], p["-"]))
			push!(h.V[i, j], (0.5 * J * p["-"], p["+"]))
			push!(h.V[i, j], (0.25 * Jzz * p["z"], p["z"]))
		end
	end	
	return h	
end

function ising2D(m::Int, n::Int; J::Real=1, hz::Real=1, periodic::Bool=false)
	p = spin_half_matrices()
	h = SquareLatticeHamiltonian(Float64, m, n)
	I2 = one(zeros(2,2))
	sx, sz = p["x"], p["z"]
	# one body terms are contained in the horizontal term
	if periodic
		for i in 1:size(h,1)
			for j in 1:size(h,2)
				push!(h.H[i, j], (J * sz, sz))
				push!(h.H[i, j], (0.25*hz*sx, I2))
				push!(h.H[i, j], (I2, 0.25*hz*sx))
			end
		end	
	else
		for i in 1:size(h,1)
			for j in 1:size(h,2)-1
				push!(h.H[i, j], (J * sz, sz))
				if j == 1
					push!(h.H[i, j], (0.5 * hz * sx, I2))
					push!(h.H[i, j], (I2, 0.25*hz*sx))
				elseif j == size(h,2)-1
					push!(h.H[i, j], (0.25*hz*sx, I2))
					push!(h.H[i, j], (I2, 0.5*hz*sx))
				else
					push!(h.H[i, j], (0.25*hz*sx, I2))
					push!(h.H[i, j], (I2, 0.25*hz*sx))
				end
			end
		end		
	end

	if periodic
		for i in 1:size(h,1)
			for j in 1:size(h,2)
				push!(h.V[i, j], (J * sz, sz))
				push!(h.V[i, j], (0.25*hz*sx, I2))
				push!(h.V[i, j], (I2, 0.25*hz*sx))
			end
		end	
	else
		for i in 1:size(h,1)-1
			for j in 1:size(h,2)
				push!(h.V[i, j], (J * sz, sz))
				if i == 1
					push!(h.V[i, j], (0.5 * hz * sx, I2))
					push!(h.V[i, j], (I2, 0.25*hz*sx))
				elseif i == size(h,2)-1
					push!(h.V[i, j], (0.25*hz*sx, I2))
					push!(h.V[i, j], (I2, 0.5 * hz*sx))
				else
					push!(h.V[i, j], (0.25*hz*sx, I2))
					push!(h.V[i, j], (I2, 0.25*hz*sx))
				end
			end
		end	
	end

	return h	
end
