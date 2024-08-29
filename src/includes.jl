using LinearAlgebra, TensorOperations


abstract type AbstractPEPSUpdateAlgorithm end
abstract type AbstractPEPSGroundStateAlgorithm end


# definition of periodic array
include("auxiliary/periodicarray.jl")
include("auxiliary/distance.jl")
include("auxiliary/mptools/mptools.jl")

# definition of lattice Hamiltonian and (first order) evolutionary operator
include("operators/operators.jl")


# definition of PEPS
include("states/states.jl")


# environments
include("envs/boundarypeps/boundarypeps.jl")
include("envs/rowenvs/rowenvs.jl")
include("envs/beliefblock/beliefblock.jl")


# algorithms
# simple update
include("algorithms/simpleupdate/simpleupdate.jl")
include("algorithms/simpleupdate/rdms.jl")

# definition of classical models
include("algorithms/classicalmodels.jl")
include("algorithms/magnetizations.jl")

# boundarymps algorithm of nonperiodic peps
include("algorithms/boundarymps/update.jl")
include("algorithms/boundarymps/expecs.jl")
include("algorithms/boundarymps/classical_expecs.jl")


# block belief propagation update
include("algorithms/blockbp/blockoperator.jl")
include("algorithms/blockbp/operator_splitting.jl")
include("algorithms/blockbp/observer_splitting.jl")
include("algorithms/blockbp/blockbp.jl")
include("algorithms/blockbp/update.jl")
include("algorithms/blockbp/expecs.jl")

# classical
include("algorithms/blockbp/classical_expecs.jl")

# a high level interface
include("algorithms/groundstate.jl")

# utilities
include("utilities/utilities.jl")