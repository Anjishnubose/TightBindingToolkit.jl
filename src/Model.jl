include("UnitCell.jl")
include("BZ.jl")
include("Hamiltonian.jl")

using LinearAlgebra


@doc """
`Model` is a data type representing a general Tight Binding system.

# Attributes
 - `uc      ::  UnitCell`: the Unit cell of the lattice.
 - `bz      ::  BZ`: The discretized Brillouin Zone.
 - `Ham     ::  Hamiltonian`: the Hamiltonian at all momentum-points.
 - `T       ::  Float64`: the temperature of the system.
 - `filling ::  Float64`: The filling of the system.
 - `mu      ::  Float64`: The chemical potential of the system.
 - `stat    ::  Int64` : ±1 for bosons and fermions.
 - `gap     ::  Float64` : the energy gap of excitations at the given filling.
 - `Gk      ::  Matrix{Matrix{ComplexF64}}` : A matrix (corresponding to the matrix of k-points in `BZ`) of Greens functions.

Initialize this structure using 
```julia
Model(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=0.0, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)
```
You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
"""
mutable struct Model
    uc      ::  UnitCell
    bz      ::  BZ
    Ham     ::  Hamiltonian
    T       ::  Float64         ##### Temperature
    filling ::  Float64         ##### Filling fraction
    mu      ::  Float64         ##### Chemical potential
    stat    ::  Int64           ##### +1 for bosons, -1 for fermions
    gap     ::  Float64
    Gk      ::  Matrix{Matrix{ComplexF64}}
    Gr      ::  Matrix{Matrix{ComplexF64}}
    
    Model(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=0.0, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1) = new{}(uc, bz, Ham, T, filling, mu, stat, 0.0, Array{Matrix{ComplexF64}}(undef, 0, 0), Array{Matrix{ComplexF64}}(undef, 0, 0))
    ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
end


@doc """
```julia
dist(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
```
Distribution function. `stat`=1 --> Bose-Einstein distribution, and `stat`=-1 --> Fermi distribution.

"""
function dist(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
    return @. 1 / (exp((E - mu) / T) - stat)
end


@doc """
```julia
distDer(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
```
derivative of [`dist`](@ref) w.r.t the energy. `stat`=1 --> Bose-Einstein distribution, and `stat`=-1 --> Fermi distribution.

"""
function distDer(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
    return @. - (1/T) * exp((E - mu) / T) / ((exp((E - mu) / T) - stat) ^ 2)
end


@doc """
```julia
findFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1) --> Float64
```
Find filling at given `T`=temperature and `mu`=chemical potential, for a given `bands`.

"""
function findFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1) :: Float64
    return sum(dist(bands; T=T, mu=mu, stat=stat)) / length(bands)
end


@doc """
```julia
getMu!(M::Model, tol::Float64=0.001)
```
Function to get chemical potential for a given `Model`, within a tolerance.

"""
function getMu!(M::Model, tol::Float64=0.001)
    energies    =   sort(reduce(vcat, M.Ham.bands))
    if M.T≈0.0
        M.mu      =   (energies[floor(Int64, length(energies)*M.filling)] + energies[floor(Int64, length(energies)*M.filling)+1])/2
    else
        muExt       =   collect(extrema(energies))
        guess       =   Nothing
        steps       =   Int(ceil(log2((muExt[end]-muExt[1])/(tol))))
        for i in 1:steps
            guess   =   mean(muExt)
            check   =   findFilling(M.Ham.bands, guess, M.T, M.stat)
            if check>M.filling
                muExt[end]  =   guess
            elseif check<M.filling
                muExt[1]    =   guess
            else
                break
            end
        end
        M.mu  =   guess
    end
end


@doc """
```julia
getFilling!(M::Model)
```
Find filling for a given `Model`.

"""
function getFilling!(M::Model)
    if M.T≈0.0
        energies    =   sort(reduce(vcat, M.Ham.bands))
        M.filling   =   searchsortedlast(energies, M.mu) / length(energies)
    else
        M.filling   =   findFilling(reduce(vcat, M.Ham.bands), M.mu, M.T, M.stat)
    end
end


@doc """
```julia
getCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) --> Matrix{Float64}
```
Function to return a diagonal matrix whose entries are M[i, i] = θ(-(E^i(k)-μ)) ----> 1 if the energy is below the chemical potential, otherwise 0.

"""
function getCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) :: Matrix{Float64}
    return diagm(dist(Es; T=T, mu=mu, stat=stat))
end


@doc """
```julia
getGk!(M::Model)
```
Finding the equal-time Greens functions in momentum space of a `Model`.

"""
function getGk!(M::Model)
    quasiCount 	=	getCount.(M.Ham.bands, Ref(M.mu), Ref(M.T), Ref(M.stat))   ##### Matrix with 1s and 0s along the diagonal. The count of the quasiparticles at each k point determined by the bands and the chemical potential
    M.Gk      =   transpose.(M.Ham.states .* quasiCount .* adjoint.(M.Ham.states))
end


@doc """
```julia
SolveModel!(M::Model)
```
one-step function to find all the attributes in Model after it has been initialized.

"""
function SolveModel!(M::Model)
    energies    =   sort(reduce(vcat, M.Ham.bands))
    if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
        getFilling!(M)
    else
        getMu!(M)
    end
    M.gap     =   energies[floor(Int64, length(energies)*M.filling) + 1] - energies[floor(Int64, length(energies)*M.filling)]
    getGk!(M)
    println("System filled!")
end

#//TODO : Add real space Greens function calculation using FFTW.