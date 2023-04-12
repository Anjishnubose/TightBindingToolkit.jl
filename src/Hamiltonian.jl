include("UnitCell.jl")
include("BZ.jl")

using LinearAlgebra

"""
Function to fill the Hamiltonian at a given momentum point, k,  given a unit cell object
# """
function FillHamiltonian(uc::UnitCell, k::Vector{Float64})
    dims    =   uc.localDim * length(uc.basis)
    H       =   zeros(ComplexF64, dims, dims)

    for site in 1:length(uc.basis)
        b1  =   uc.localDim * (site - 1) + 1
        b2  =   uc.localDim * (site - 1) + 1
        H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   sum(uc.fields[site] .* SpinVec) 
    end

    for bond in uc.bonds
        b1  =   uc.localDim * (bond.base - 1) + 1
        b2  =   uc.localDim * (bond.target - 1) + 1
        H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  +=   exp( im * dot(k, sum(bond.offset .* uc.primitives))) * bond.mat
        H[b2 : b2 + uc.localDim - 1, b1 : b1 + uc.localDim - 1]  +=   exp(-im * dot(k, sum(bond.offset .* uc.primitives))) * bond.mat'  
    end

    return H
end

"""
Function to get H over all k points
"""
function FullHamiltonian(uc::UnitCell, bz::BZ)
    return FillHamiltonian.(Ref(uc), bz.ks)
end

"""
The Hamiltonian struct contains the following
    1. H : a matrix of matrix s.t. H[i, j] corresponds to the Hamiltonian at k=ks[i, j] where ks=bz.ks in the brillouin Zone
    2. bands : a matrix of vectors s.t. bands[i, j] are the energies at k=ks[i, j]
    3. states : a matrix of unitary matrices s.t. states[i, j] are all the eigenvectors at k=ks[i, j]
    4. Gk :  a matrix of matrix s.t. Gk[i, j] are all the correlations <c†_k . c_k> b/w the sublattices at k=ks[i, j] ----> momentum space correlations
    5. Gr : a matrix of matrix s.t. Gr[i, j] = <c†_0 c_δ> where δ = i * a1 + j * a2  ----> real space correlations
"""
mutable struct Hamiltonian
    H       ::  Matrix{Matrix{ComplexF64}}
    bands   ::  Matrix{Vector{Float64}}
    states  ::  Matrix{Matrix{ComplexF64}}
    filling ::  Float64
    Gk      ::  Matrix{Matrix{ComplexF64}}
    Gr      ::  Matrix{Matrix{ComplexF64}}

    Hamiltonian(uc::UnitCell, bz::BZ, filling::Float64=0.5) = new{}(FullHamiltonian(uc, bz), Array{Vector{Float64}}(undef, 0, 0), 
                                                Array{Matrix{Float64}}(undef, 0, 0), filling, Array{Matrix{Float64}}(undef, 0, 0), Array{Matrix{Float64}}(undef, 0, 0))
end

"""
Diagonalize a Hamiltonian at all momentum points simultaneuosly!
"""
function Diagonalize(H::Matrix{Matrix{ComplexF64}})
    sols    =   eigen.(Hermitian.(H))
    states  =   getfield.(sols, :vectors)
    bands   =   getfield.(sols, :values)
    return bands, states
end

"""
Finding the Greens functions in momentum space at some filling
"""
function getGk(states::Matrix{Matrix{ComplexF64}}, filling::Float64)
    N           =   size(states[1, 1])[1]  ##### each Gk is going to be a N x N complex matrix
    quasiCount 	=	diagm(cat(ones(round(Int64, N * filling)), zeros(round(Int64, N - (N * filling))), dims=1))   ##### Matrix with 1s and 0s along the diagonal s.t. the trace is as close to N * filling as possible.
    Gk          =   states .* Ref(quasiCount) .* adjoint.(states)
    return tranpose.(Gk)
end

function SolveHamiltonian!(Ham::Hamiltonian)
    Ham.bands, Ham.states   =   Diagonalize(Ham.H)
    Ham.Gk      =   getGk(Ham.states, Ham.filling)
end

"""
An example!
"""

if abspath(PROGRAM_FILE) == @__FILE__
    ##### Square lattice unit cell with enlarged unit cell along one direction
    a1      =   [2.0, 0.0]
    a2      =   [0.0, 1.0]
    uc      =   UnitCell(a1, a2, 2)
    addBasisSite!(uc, [0.0, 0.0])
    addBasisSite!(uc, [1.0, 0.0])
    addIsotropicBonds!(uc, 1.0, SpinVec[3], "1NN", 1)
    ##### BZ with size 100
    kSize   =   100
    bz      =   BZ(kSize)
    fillBZ!(bz, uc)
    ##### timing stuff
    @time Hamiltonian(uc, bz)
    @time Hamiltonian(uc, bz)
    @time Hamiltonian(uc, bz)

end
