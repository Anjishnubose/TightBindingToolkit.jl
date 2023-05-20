include("UnitCell.jl")
include("BZ.jl")
include("SpinMats.jl")
using LinearAlgebra

"""
Function to fill the Hamiltonian at a given momentum point, k,  given a unit cell object
# """
function FillHamiltonian(uc::UnitCell, k::Vector{Float64} ; SpinMatrices::Vector{Matrix{ComplexF64}})
    dims    =   uc.localDim * length(uc.basis)
    H       =   zeros(ComplexF64, dims, dims)

    for site in 1:length(uc.basis)
        b1  =   uc.localDim * (site - 1) + 1
        b2  =   uc.localDim * (site - 1) + 1
        H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   sum(uc.fields[site][1:3] .* SpinMatrices[1:3]) 
    end

    for bond in uc.bonds
        b1  =   uc.localDim * (bond.base - 1) + 1
        b2  =   uc.localDim * (bond.target - 1) + 1
        
        if b1==b2 && bond.offset==zeros(length(uc.basis))
            H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   bond.mat
        else
            H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   exp( im * dot(k, sum(bond.offset .* uc.primitives))) * bond.mat
            H[b2 : b2 + uc.localDim - 1, b1 : b1 + uc.localDim - 1]  .+=   exp(-im * dot(k, sum(bond.offset .* uc.primitives))) * bond.mat' 
        end 
    end

    return H
end

"""
Function to get H over all k points
"""
function FullHamiltonian(uc::UnitCell, bz::BZ)
    SpinMatrices    =   SpinMats((uc.localDim-1)//2)
    return FillHamiltonian.(Ref(uc), bz.ks ; SpinMatrices = SpinMatrices)
end

"""
The Hamiltonian struct contains the following
    1. H : a matrix of matrix s.t. H[i, j] corresponds to the Hamiltonian at k=ks[i, j] where ks=bz.ks in the brillouin Zone
    2. bands : a matrix of vectors s.t. bands[i, j] are the energies at k=ks[i, j]
    3. states : a matrix of unitary matrices s.t. states[i, j] are all the eigenvectors at k=ks[i, j]

"""
mutable struct Hamiltonian
    H           ::  Matrix{Matrix{ComplexF64}}
    bands       ::  Matrix{Vector{Float64}}
    states      ::  Matrix{Matrix{ComplexF64}}
    bandwidth   ::  Tuple{Float64, Float64}

    Hamiltonian(uc::UnitCell, bz::BZ) = new{}(FullHamiltonian(uc, bz), Array{Vector{Float64}}(undef, 0, 0), Array{Matrix{ComplexF64}}(undef, 0, 0), (0.0, 0.0))
end

"""
Diagonalize a Hamiltonian at all momentum points simultaneuosly!
"""
function DiagonalizeHamiltonian!(Ham::Hamiltonian)
    sols            =   eigen.(Hermitian.(Ham.H))
    Ham.bands       =   getfield.(sols, :values)
    Ham.states      =   getfield.(sols, :vectors)
    Ham.bandwidth   =   extrema(reduce(vcat, Ham.bands))
    println("Hamiltonian Diagonalized")
end

"""
Density of state calculator
"""
function DOS(Omega::Float64, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) :: Float64

    energies    =   reduce(vcat, Ham.bands)
    n_bands     =   length(Ham.bands[1, 1])
    energies    =   energies[filter(i -> (i-1) % n_bands + 1 <= till_band, 1:length(energies))]
    dos         =   @. 1 / (Omega - energies + im * spread)
    dos         =   imag(sum(dos))
    return dos
end

function DOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) :: Vector{Float64}

    dos     =   DOS.(Omegas, Ref(Ham); till_band=till_band, spread=spread)
    dos     =   @. dos / sum(dos)
    return dos
end