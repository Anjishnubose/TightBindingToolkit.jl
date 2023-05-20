include("Model.jl")
include("SpinMats.jl")

using LinearAlgebra
using TensorCast
using MappedArrays

directions = ["x" , "y" , "z"]

"""
Structure containing all the relevant info regarding magnetic susceptibility
"""
mutable struct susceptibility
    M       ::  Model   
    Qs      ::  Vector{Vector{Float64}}
    Omegas  ::  Vector{Float64}
    Spread  ::  Float64
    chis    ::  Dict

    susceptibility(M::Model , Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, [], Omegas, eta, Dict())
    susceptibility(M::Model , Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, Qs, Omegas, eta, Dict())
end


"""
function to calculate (nF(E(k)) - nF(E(k+Q))) / (Ω + i * η - (E(k+Q) - E(k))) for a general multi-band system, at all k-points
"""
function nF_factor(Q_index::Vector{Int64}, Omega::Float64 , M::Model; eta::Float64=1e-2) :: Matrix{Matrix{ComplexF64}}
    
    Es_k    =   M.Ham.bands
    Es_kpQ  =   circshift(Es_k, -Q_index)  ##### The bands with momentum indices rolled to k+Q instead of k.

    ##### The fermi distribution functions for each band at each k-point
    nf_k    =   dist.(Es_k ; T=M.T, mu=M.mu, stat=M.stat) 
    nf_kpQ  =   circshift(nf_k, -Q_index)

    @cast nF[k1, k2][i, j] |= (nf_k[k1, k2][i] - nf_kpQ[k1, k2][j]) / (Omega + im * eta - (Es_kpQ[k1, k2][j] - Es_k[k1, k2][i]));
    ##### Returning a matrix of matrix type
    return nF
end 

function FindChi(Q::Vector{Float64}, Omega::Float64 , M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2) ::ComplexF64

    SpinVec =   SpinMats((M.uc.localDim-1)//2)
    expQ    =   diagm(exp.(-im .* dot.(Ref(Q) , M.uc.basis)))
    Qa      =   kron(expQ , SpinVec[a])
    Qb      =   kron(expQ , SpinVec[b])

    Q_index =   GetQIndex(Q, M.bz) .- [(M.bz.gridSize + 1)÷2 , (M.bz.gridSize + 1)÷2]
    nF      =   nF_factor(Q_index, Omega , M; eta=eta)

    Uk      =   M.Ham.states
    UkQ     =   circshift(M.Ham.states, -Q_index)

    if a!=b
        Ma      =   adjoint.(UkQ) .* Ref(Qa) .* Uk
        Mb      =   adjoint.(UkQ) .* Ref(Qb) .* Uk
        chi     =   -tr(sum(map( .* , Ma, mappedarray(conj , Mb)) .* nF))
    else
        Ma      =   adjoint.(UkQ) .* Ref(Qa) .* Uk
        chi     =   -tr(sum(map( .* , Ma, mappedarray(conj , Ma)) .* nF))
    end

    return chi 
end

function FindChi_FullBZ(Omega::Float64 , M::Model ; a::Int64=3, b::Int64=3, eta::Float64=1e-2) :: Matrix{ComplexF64}
    return FindChi.(M.bz.ks, Ref(Omega), Ref(M) ; a=a, b=b, eta=eta)
end

function FindChi_Path(Omega::Float64 , M::Model, path::Vector{Vector{Float64}} ; a::Int64=3, b::Int64=3, eta::Float64=1e-2) :: Vector{ComplexF64}
    return FindChi.(path, Ref(Omega), Ref(M) ; a=a, b=b, eta=eta)
end


function FillChis!(chi::susceptibility; fill_BZ::Bool=false, a::Int64=3, b::Int64=3)
    @assert chi.M.bz.gridSize%2==1  "Only works if the Gamma point is included in the grid ---> grid size must be odd!"

    if fill_BZ
        chis    =   FindChi_FullBZ.(chi.Omegas, Ref(chi.M) ; a=a, b=b, eta=chi.Spread)
        @cast chisNew[i, j, k] |= chis[i][j, k] / length(bz.ks);
    else
        chis    =   FindChi_Path.(chi.Omegas, Ref(chi.M), Ref(chi.Qs) ; a=a, b=b, eta=chi.Spread)
        @cast chisNew[i, j] |= chis[i][j] / length(bz.ks);
    end

    chi.chis[directions[a] * directions[b]]     =   chisNew
    println("Chis filled along " * directions[a] * directions[b])
end



