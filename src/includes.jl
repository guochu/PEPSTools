using Base: @boundscheck, front, tail
using Random, LinearAlgebra, TensorOperations, Strided
using Zygote
using Flux, NNQS
import NNQS: energy


abstract type PEPSUpdateAlgorithm end
abstract type ImaginaryTimePEPSUpdateAlgorithm <: PEPSUpdateAlgorithm end
abstract type PEPSGroundStateAlgorithm end


# definition of periodic array
include("auxiliary/periodicarray.jl")
include("auxiliary/distance.jl")
include("auxiliary/mptools/mptools.jl")
include("auxiliary/bputility/bputility.jl")

# simple definition of message
include("message.jl")

# definition of lattice Hamiltonian and (first order) evolutionary operator
include("operators/operators.jl")


# definition of PEPS
include("states/states.jl")


## environments
# border PEPS
include("borderedpeps/borderedpeps.jl")
# BMPS environment
include("bmpsenv/bmpsenv.jl")
# BlockBP environment
include("blockbpenv/blockbpenv.jl")
# BP environment
include("bpenv/bpenv.jl")


# algorithms
# simple update
include("algorithms/simpleupdate/simpleupdate.jl")
include("algorithms/simpleupdate/rdms.jl")


# boundarymps algorithm of nonperiodic peps
include("algorithms/boundarymps/update.jl")
include("algorithms/boundarymps/expecs.jl")
include("algorithms/boundarymps/rdms.jl")


# block belief propagation update
include("algorithms/blockbp/blockoperator.jl")
include("algorithms/blockbp/operator_splitting.jl")
include("algorithms/blockbp/observer_splitting.jl")
include("algorithms/blockbp/blockbp.jl")
include("algorithms/blockbp/update.jl")
include("algorithms/blockbp/expecs.jl")
include("algorithms/blockbp/rdms.jl")

# bp update
include("algorithms/bp/variationalbp/variationalbp.jl")
include("algorithms/bp/vmcbp/vmcbp.jl")

# classical models
include("algorithms/classicalmodels/classicalmodels.jl")


# a high level interface
include("algorithms/groundstate.jl")

# utilities
include("utilities/utilities.jl")