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


## environments
# border PEPS
include("borderedpeps/borderedpeps.jl")
# BMPS environment
include("bmpsenv/bmpsenv.jl")
# BlockBP partition
include("blockbppartition/blockbppartition.jl")


# algorithms
# simple update
include("algorithms/simpleupdate/simpleupdate.jl")
include("algorithms/simpleupdate/rdms.jl")


# boundarymps algorithm of nonperiodic peps
include("algorithms/boundarymps/update.jl")
include("algorithms/boundarymps/expecs.jl")


# block belief propagation update
include("algorithms/blockbp/blockoperator.jl")
include("algorithms/blockbp/operator_splitting.jl")
include("algorithms/blockbp/observer_splitting.jl")
include("algorithms/blockbp/blockbp.jl")
include("algorithms/blockbp/update.jl")
include("algorithms/blockbp/expecs.jl")

# classical models
include("algorithms/classicalmodels/classicalmodels.jl")


# a high level interface
include("algorithms/groundstate.jl")

# utilities
include("utilities/utilities.jl")
