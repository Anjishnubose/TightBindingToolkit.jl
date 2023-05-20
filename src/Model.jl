include("UnitCell.jl")
include("BZ.jl")
include("Hamiltonian.jl")

using LinearAlgebra


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


"""
Distribution Function : E can be a number, vector, array etc
"""
function dist(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
    return @. 1 / (exp((E - mu) / T) - stat)
end

function distDer(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
    return @. - (1/T) * exp((E - mu) / T) / ((exp((E - mu) / T) - stat) ^ 2)
end

"""
Find filling at given temperature and chemical potential
"""
function findFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1)
    return sum(dist(bands; T=T, mu=mu, stat=stat)) / length(bands)
end

"""
Function to get chemical potential given a filling!
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
"""
Function to get filling given a chemical potential!
"""
function getFilling!(M::Model)
    if M.T≈0.0
        energies    =   sort(reduce(vcat, M.Ham.bands))
        M.filling   =   searchsortedlast(energies, M.mu) / length(energies)
    else
        M.filling   =   findFilling(reduce(vcat, M.Ham.bands), M.mu, M.T, M.stat)
    end
end
"""
Function to return a diagonal matrix whose entries are M[i, i] = θ(-(E^i(k)-μ)) ----> 1 if the energy is below the chemical potential, otherwise 0.
"""
function getCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64)
    return diagm(dist(Es; T=T, mu=mu, stat=stat))
end
"""
Finding the Greens functions in momentum space at some chemical potential
"""
function getGk!(M::Model)
    quasiCount 	=	getCount.(M.Ham.bands, Ref(M.mu), Ref(M.T), Ref(M.stat))   ##### Matrix with 1s and 0s along the diagonal. The count of the quasiparticles at each k point determined by the bands and the chemical potential
    M.Gk      =   transpose.(M.Ham.states .* quasiCount .* adjoint.(M.Ham.states))
end


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