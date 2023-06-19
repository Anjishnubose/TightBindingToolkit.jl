module BdG

using ..TightBindingToolkit.UCell:UnitCell, ModifyIsotropicFields!
using ..TightBindingToolkit.BZone:BZ
using ..TightBindingToolkit.Hams:Hamiltonian, ModifyHamiltonianField!
using ..TightBindingToolkit.TBModel:dist, distDer, BinarySearch

import ..TightBindingToolkit.TBModel:findFilling, getMu!, getFilling!, getGk!, SolveModel!

export BdGModel, findFilling, getMu!, getFilling!, getGk!, SolveModel!

using LinearAlgebra, Tullio, Statistics, TensorCast

@doc """
    `BdGModel` is a data type representing a general Tight Binding system with pairing.

    # Attributes
    - `uc_hop  ::  UnitCell`: the Unit cell of the lattice with the hoppings.
    - `uc_pair ::  UnitCell`: the Unit cell of the lattice with the pairings.
    - `bz      ::  BZ`: The discretized Brillouin Zone.
    - `Ham     ::  Hamiltonian`: the Hamiltonian at all momentum-points.
    - `T       ::  Float64`: the temperature of the system.
    - `filling ::  Float64`: The filling of the system.
    - `mu      ::  Float64`: The chemical potential of the system.
    - `stat    ::  Int64` : ±1 for bosons and fermions.
    - `gap     ::  Float64` : the energy gap of excitations at the given filling.
    - `Gk      ::  Array{Matrix{ComplexF64}}` : An Array (corresponding to the grid of k-points in `BZ`) of Greens functions.
    - `Fk      ::  Array{Matrix{ComplexF64}}` : An Array (corresponding to the grid of k-points in `BZ`) of anomalous Greens functions.

    Initialize this structure using 
    ```julia
    BdGModel(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=0.0, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)
    ```
    You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
    """
mutable struct BdGModel
    uc_hop  ::  UnitCell
    uc_pair ::  UnitCell
    bz      ::  BZ
    Ham     ::  Hamiltonian
    T       ::  Float64         ##### Temperature
    filling ::  Float64         ##### Filling fraction
    mu      ::  Float64         ##### Chemical potential
    gap     ::  Float64
    stat    ::  Int64           ##### +1 for Boson, -1 for Fermions
    Gk      ::  Array{Matrix{ComplexF64}}
    Fk      ::  Array{Matrix{ComplexF64}}
    Gr      ::  Array{Matrix{ComplexF64}}
    Fr      ::  Array{Matrix{ComplexF64}}
    
    BdGModel(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=0.0, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1) = new{}(uc_hop, uc_pair, bz, Ham, T, filling, mu, 0.0, stat ,Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...),Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...))
    ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
end


@doc """
```julia
findFilling(M::BdGModel) --> Float64
findFilling(mu::Float64, M::BdGModel)
```
Find filling at given temperature and chemical potential that takes BdGModel object as argument instead of Hamiltonian, since for a BdG case, the filling depends on the wavefunctions also.
Because of this, if you want to calculate the filling at a different chemical potential, you have to modify the Hamiltonian, the UnitCells, rediagonalize, and recalculate eveyrthing.

"""
function findFilling(M::BdGModel) :: Float64
    @assert M.Ham.is_BdG==true "Use other format for pure hopping Hamiltonian"

    N       =   length(M.Ham.bands[begin])
    Eks     =   reshape(getindex.(M.Ham.bands, Ref(1:N÷2)), prod(M.bz.gridSize))
    U11s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(1:N÷2)), prod(M.bz.gridSize))
    U21s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(1:N÷2)), prod(M.bz.gridSize))

    nFs     =   dist.(Eks; T=M.T, mu=0.0, stat=M.stat)

    @tullio filling := ((conj(U11s[k1][i, j]) * U11s[k1][i, j] * nFs[k1][j]) 
                      + (conj(U21s[k1][i, j]) * U21s[k1][i, j] * (1 - nFs[k1][j])))

    return real(filling) / (prod(M.bz.gridSize) * M.uc_hop.localDim * length(M.uc_hop.basis))
end

function findFilling(mu::Float64, M::BdGModel) :: Float64
    ModifyHamiltonianField!(M.Ham, M.uc_hop, repeat([mu], outer=length(M.uc_hop.basis)) ; dim=4)
    return findFilling(M)
end


"""
```julia
getFilling!(M::BdGModel)
```
Function to get filling given a chemical potential!

"""
function getFilling!(M::BdGModel)
    M.filling   =   findFilling(M.mu, M)

end


"""
```julia
getMu!(M::BdGModel ;  tol::Float64=0.001)
```
Function to get the correct chemical potential given a filling.
"""
function getMu!(M::BdGModel ;  tol::Float64=0.001)
    M.mu    =   BinarySearch(M.filling, M.Ham.bandwidth, findFilling, (M,) ; tol=tol)
    println("Found chemical potential for given filling.")
    

end


"""
```julia
getGk!(M::BdGModel)
```
Finding the Greens functions, and anomalous greens functions in momentum space at some chemical potential.
"""
function getGk!(M::BdGModel)

    N       =   length(M.Ham.bands[begin])
    Eks     =   reshape(getindex.(M.Ham.bands, Ref(1:N÷2)) , prod(M.bz.gridSize))   ##### Only the negative energies from the bdG spectrum

    U11s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(1:N÷2)) , prod(M.bz.gridSize))          ##### The 4 different quadrants in the Unitary relating the BdG quasiparticles and the nambu basis
    U21s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(1:N÷2)) , prod(M.bz.gridSize)) 
    U12s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(N÷2 + 1: N)) , prod(M.bz.gridSize)) 
    U22s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(N÷2 + 1: N)) , prod(M.bz.gridSize))

    nFs     =   dist.(Eks; T=M.T, mu=0.0, stat=M.stat)
    #//TODO CHECK!!!!
    @reduce Gk[k1][i, j] |= sum(l) ((conj(U11s[k1][i, l]) * U11s[k1][j, l] * nFs[k1][l]) 
                                      + (conj(U12s[k1][i, l]) * U12s[k1][j, l] * (1 - nFs[k1][l])))

    @reduce Fk[k1][i, j] |= sum(l) ((conj(U11s[k1][i, l]) * U21s[k1][j, l] * nFs[k1][l]) 
                                      + (conj(U12s[k1][i, l]) * U22s[k1][j, l] * (1 - nFs[k1][l])))
    M.Gk    =   reshape(Gk, M.bz.gridSize...)
    M.Fk    =   reshape(Fk, M.bz.gridSize...)

end

"""
```julia
SolveModel!(M::BdGModel)
```
one-step function to find all the attributes in  BdGModel after it has been initialized.
"""
function SolveModel!(M::BdGModel)
    @assert M.Ham.is_BdG==true "Use other format for pure hopping Hamiltonian"

    if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
        getFilling!(M)
    else
        getMu!(M)
    end

    energies    =   sort(reduce(vcat, M.Ham.bands))
    M.gap     =   energies[floor(Int64, length(energies)*M.filling) + 1] - energies[floor(Int64, length(energies)*M.filling)]

    getGk!(M)
    println("System filled!")
end


end