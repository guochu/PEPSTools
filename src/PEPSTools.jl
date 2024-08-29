module PEPSTools

# # some basic definitions
# simple definition of Hamiltonian and operators
export scalartype
export SquareLatticeBonds, SquareLatticeHamiltonian, SquareLatticeOperator, is_periodic, is_nonperiodic
export squeeze, exponential
export SquareLatticeSites, LocalCObservers, LocalQObservers


# definition of PEPS
export PEPS, randompeps, prodpeps, randomsquaretn
export bond_dimensions, physical_dimensions, bond_dimension


# environments

# export SquareLatticePartition, block_partition, peps_partition

# # algorithms

# # classical tensor network
# export LocalObservers, MagnetizationTensors, Classical2DModel, ClassicalIsing2D, magnetizations, magnetization, interactionH


# export PEPSBlock, nrows, ncols
# export BoundaryMPS, sweep!, expectation, row_expectations, local_expectations, expectationfull
# export CanonicalPEPS, PEPSSimpleUpdate
# # reduced density matrices
# export rdm1s, rdm2s



# # block bp algorithm
# export Message, BlockBP, BeliefPEPSBlock, BlockOperator, BlockLocalOperator, default_splitting, center_splitting, subblock, subblocks

# # high-level
# export ImaginaryTimePEPS, ground_state!

# export densevector

# export heisenberg2D, ising2D



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

end

