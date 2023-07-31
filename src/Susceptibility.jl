module suscep
    export Susceptibility , FindChi , FillChis!

    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.Useful: DistFunction
    using ..TightBindingToolkit.BZone:BZ, GetQIndex
    using ..TightBindingToolkit.TBModel:Model


    using LinearAlgebra, TensorCast, Tullio, Logging

    directions = ["x" , "y" , "z"]


@doc """
`Susceptibility` is a data type representing the magnetic response, ``χ^{ab}(Q , Ω)`` for a general tight-binding `Model`.

# Attributes
- `M       ::  Model `: the given model.
- `Qs      ::  Vector{Vector{Float64}}`: the set of momentum points over which ``χ^{ab}(Q , Ω)`` is calculated.
- `Omegas  ::  Vector{Float64}`: the set of energies over which ``χ^{ab}(Q , Ω)`` is calculated.
- `Spread  ::  Float64` : the finite spread when summing over delta functions.
- `chis    ::  Dict`: a dictionary containing ``χ^{ab}(Q , Ω)`` for the different directions e.g. `chis["xx"]` etc.

Initialize this structure using 
```julia
Susceptibility(M::Model , Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, [], Omegas, eta, Dict())
Susceptibility(M::Model , Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, Qs, Omegas, eta, Dict())
```
"""
    mutable struct Susceptibility
        M       ::  Model   
        Qs      ::  Vector{Vector{Float64}}
        Omegas  ::  Vector{Float64}
        Spread  ::  Float64
        chis    ::  Dict

        Susceptibility(M::Model , Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, [], Omegas, eta, Dict())
        Susceptibility(M::Model , Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, Qs, Omegas, eta, Dict())
    end


@doc raw"""
```julia
nF_factor(Q_index::Vector{Int64}, Omega::Float64 , M::Model; eta::Float64=1e-2) --> Matrix{Matrix{ComplexF64}}
```
function to calculate ``(nF(E(k)) - nF(E(k+Q))) / (Ω + i * η - (E(k+Q) - E(k)))`` for a general multi-band system, at all k-points.

"""
    function nF_factor(Q_index::Vector{Int64}, Omega::Float64 , M::Model; eta::Float64=1e-2) :: Matrix{Matrix{ComplexF64}}
        
        Es_k    =   M.Ham.bands
        Es_kpQ  =   circshift(Es_k, -Q_index)  ##### The bands with momentum indices rolled to k+Q instead of k.

        ##### The fermi distribution functions for each band at each k-point
        nf_k    =   DistFunction.(Es_k ; T=M.T, mu=M.mu, stat=M.stat) 
        nf_kpQ  =   circshift(nf_k, -Q_index)

        @cast nF[k1, k2][i, j] |= (nf_k[k1, k2][i] - nf_kpQ[k1, k2][j]) / (Omega + im * eta - (Es_kpQ[k1, k2][j] - Es_k[k1, k2][i]));
        ##### Returning a matrix of matrix type
        return nF
    end 


@doc """
```julia
FindChi(Q::Vector{Float64}, Omega::Float64 , M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2) --> ComplexF64
```
function to calculate susceptibility at a fixed Ω=`Omega`, and `Q`, and along a fixed direction given by `a` and `b`.

"""
    function FindChi(Q::Vector{Float64}, Omega::Float64 , M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2, SpinMatrix::Vector{Matrix{ComplexF64}} = SpinMats((M.uc.localDim-1)//2)) ::ComplexF64

        expQ    =   diagm(exp.(-im .* dot.(Ref(Q) , M.uc.basis)))
        Qa      =   kron(expQ , SpinMatrix[a])
        Qb      =   kron(expQ , SpinMatrix[b])
        @assert size(Qa) == size(Qb) == (length(M.uc.basis) * M.uc.localDim, length(M.uc.basis) * M.uc.localDim)

        Q_index =   GetQIndex(Q, M.bz) .- ((M.bz.gridSize .+ 1) .÷ 2)
        nF      =   nF_factor(Q_index, Omega , M; eta=eta)

        Uk      =   M.Ham.states
        UkQ     =   circshift(M.Ham.states, -Q_index)

        if a!=b
            Ma      =   adjoint.(UkQ) .* Ref(Qa) .* Uk
            Mb      =   adjoint.(UkQ) .* Ref(Qb) .* Uk
            @tullio chi     :=   - Ma[k1, k2][i, j] * conj(Mb[k1, k2][i, j]) * nF[k1, k2][j, i]
        else
            Ma      =   adjoint.(UkQ) .* Ref(Qa) .* Uk
            @tullio chi     :=   - Ma[k1, k2][i, j] * conj(Ma[k1, k2][i, j]) * nF[k1, k2][j, i]
        end

        return chi 
    end


@doc """
```julia
FindChi_FullBZ(Omega::Float64 , M::Model ; a::Int64=3, b::Int64=3, eta::Float64=1e-2) --> Matrix{ComplexF64}
```
function to calculate susceptibility at a fixed Ω=`Omega`, but for all `Q` present in the `BZ`, and along a fixed direction given by `a` and `b`.

"""
    function FindChi_FullBZ(Omega::Float64 , M::Model ; a::Int64=3, b::Int64=3, eta::Float64=1e-2, SpinMatrix::Vector{Matrix{ComplexF64}} = SpinMats((M.uc.localDim-1)//2)) :: Matrix{ComplexF64}
        return FindChi.(M.bz.ks, Ref(Omega), Ref(M) ; a=a, b=b, eta=eta, SpinMatrix = SpinMatrix)
    end


@doc """
```julia
FindChi_Path(Omega::Float64 , M::Model, path::Vector{Vector{Float64}} ; a::Int64=3, b::Int64=3, eta::Float64=1e-2) --> Vector{ComplexF64}
```
function to calculate susceptibility at a fixed Ω=`Omega`, but for all `Q` present in the given path, and along a fixed direction given by `a` and `b`.

"""
    function FindChi_Path(Omega::Float64 , M::Model, path::Vector{Vector{Float64}} ; a::Int64=3, b::Int64=3, eta::Float64=1e-2, SpinMatrix::Vector{Matrix{ComplexF64}} = SpinMats((M.uc.localDim-1)//2)) :: Vector{ComplexF64}
        return FindChi.(path, Ref(Omega), Ref(M) ; a=a, b=b, eta=eta, SpinMatrix = SpinMatrix)
    end


@doc """
```julia
FillChis!(chi::Susceptibility; fill_BZ::Bool=false, a::Int64=3, b::Int64=3)
```
function to calculate susceptibility at a all given Ω=`Omegas`, but for all `Q` present in the given path, and along a fixed direction given by `a` and `b`.

"""
    function FillChis!(chi::Susceptibility; fill_BZ::Bool=false, a::Int64=3, b::Int64=3, SpinMatrix::Vector{Matrix{ComplexF64}} = SpinMats((chi.M.uc.localDim-1)//2))
        @assert prod(chi.M.bz.gridSize .% 2) == 1  "Only works if the Gamma point is included in the grid ---> grid size must be odd!"

        if fill_BZ
            chis    =   FindChi_FullBZ.(chi.Omegas, Ref(chi.M) ; a=a, b=b, eta=chi.Spread, SpinMatrix = SpinMatrix)
            @cast chisNew[i, j, k] |= chis[i][j, k] / length(chi.M.bz.ks);
        else
            chis    =   FindChi_Path.(chi.Omegas, Ref(chi.M), Ref(chi.Qs) ; a=a, b=b, eta=chi.Spread, SpinMatrix = SpinMatrix)
            @cast chisNew[i, j] |= chis[i][j] / length(chi.M.bz.ks);
        end

        chi.chis[directions[a] * directions[b]]     =   chisNew
        @info "Chis filled along $(directions[a])-$(directions[b])"
    end

end


