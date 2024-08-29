abstract type AbstractRowEnv end

abstract type AbstractPEPSRowEnv <: AbstractRowEnv end
abstract type AbstractSquareTNRowEnv <: AbstractRowEnv end


Base.length(x::AbstractRowEnv) = length(x.middle)
compute_center(x::AbstractPEPSRowEnv, center::Int; kwargs...) = error("compute_center not implemented for env type $(typeof(x))")

update_left!(x::AbstractPEPSRowEnv, pos::Int) = error("update_left! not implemented for env type $(typeof(x))")

function update_single!(x::AbstractPEPSRowEnv, pos::Int, U::AbstractArray{<:Number, 4}; kwargs...)
	Qleft, aL, Qright, aR, Xt = compute_center(x, pos)
	aLn, aRn = center_minimization(aL, aR, Xt, U; kwargs...)	

    @tensor cl[1,2,6,3,5] := Qleft[1,2,3,4] * aLn[4,5,6]
    @tensor cr[1,4,5,6,2] := aRn[1,2,3] * Qright[3,4,5,6]
    x.middle[pos] = cl
    x.middle[pos+1] = cr
    update_left!(x, pos)
end
update_single!(x::AbstractPEPSRowEnv, pos::Int, U::Nothing; kwargs...) = update_left!(x, pos)

function row_expectation_single(x::AbstractPEPSRowEnv, pos::Int, U::AbstractArray{<:Number, 4})
	Qleft, aL, Qright, aR, Xt = compute_center(x, pos)
	@tensor down[2,4,6] := aL[1,2,3] * Xt[1,4,5] * aR[3,6,5] 
	@tensor r = conj(down[1,2,3]) * U[1,3,4,5] * down[4,2,5]
	update_left!(x, pos)
	return r / dot(down, down)
end

function row_expectation_single(x::AbstractPEPSRowEnv, pos::Int, U::Nothing) 
	update_left!(x, pos)
	return 0.
end 

function row_rdm2_single(x::AbstractPEPSRowEnv, pos::Int)
	Qleft, aL, Qright, aR, Xt = compute_center(x, pos)
	@tensor down[2,4,6] := aL[1,2,3] * Xt[1,4,5] * aR[3,6,5] 
	@tensor rho[4,5,1,3] := conj(down[1,2,3]) * down[4,2,5]
	update_left!(x, pos)

	return normalize_rho!(rho)
end

