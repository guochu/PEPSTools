module PEPSTools

# Message
export Message, FixedNorm, FixedSum

# # some basic definitions
# simple definition of Hamiltonian and operators
export scalartype
export SquareLatticeBonds, SquareLatticeHamiltonian, SquareLatticeOperator, is_periodic, is_nonperiodic
export squeeze, exponential
export SquareLatticeSites, LocalCObservers, LocalQObservers
# matrix product tools
export IterativeCompression, SVDCompression


# definition of PEPS
export PEPS, randompeps, prodpeps, randomsquaretn, CanonicalPEPS
export bond_dimensions, physical_dimensions, bond_dimension, increase_bond!


# environments
export borderedpeps
export SquareLatticePartition, lattice_partition, blockbp_environments

# blockbp partition
export nrows, ncols, blockbp_environments


## algorithms
export BoundaryMPS, sweep!, expectation, energy
export SimpleUpdate
# reduced density matrices
export rdm1s, rdm2s
# block bp algorithm
export BlockBP, BlockOperator, BlockLocalOperator, default_splitting, center_splitting, subblock, subblocks
# BP algorithm
export BP

# high-level
export ImaginaryTimePEPS, ground_state!

# classical tensor network
export Classical2DModel, ClassicalIsing2D, magnetizations, magnetization, bond_energy



## utility
# predefined models
export heisenberg2D, ising2D


using Base: @boundscheck, front, tail
using Logging: @warn
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

end

