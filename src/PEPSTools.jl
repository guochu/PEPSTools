module PEPSTools

# # some basic definitions
# simple definition of Hamiltonian and operators
export scalartype
export SquareLatticeBonds, SquareLatticeHamiltonian, SquareLatticeOperator, is_periodic, is_nonperiodic
export squeeze, exponential
export SquareLatticeSites, LocalCObservers, LocalQObservers


# definition of PEPS
export PEPS, randompeps, prodpeps, randomsquaretn, CanonicalPEPS
export bond_dimensions, physical_dimensions, bond_dimension


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
export Message, BlockBP, BlockOperator, BlockLocalOperator, default_splitting, center_splitting, subblock, subblocks

# high-level
export ImaginaryTimePEPS, ground_state!

# classical tensor network
export Classical2DModel, ClassicalIsing2D, magnetizations, magnetization, bond_energy



## utility
# predefined models
export heisenberg2D, ising2D



using LinearAlgebra, TensorOperations


abstract type PEPSUpdateAlgorithm end
abstract type ImaginaryTimePEPSUpdateAlgorithm <: PEPSUpdateAlgorithm end
abstract type PEPSGroundStateAlgorithm end


# definition of periodic array
include("auxiliary/periodicarray.jl")
include("auxiliary/distance.jl")
include("auxiliary/mptools/mptools.jl")

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
# BlockBP partition
include("blockbpenv/blockbpenv.jl")


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

# classical models
include("algorithms/classicalmodels/classicalmodels.jl")


# a high level interface
include("algorithms/groundstate.jl")

# utilities
include("utilities/utilities.jl")

end

