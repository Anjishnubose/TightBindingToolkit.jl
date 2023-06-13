module Hams
    export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , DiagonalizeHamiltonian! , DOS, ModifyHamiltonianMu!

    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.UCell:UnitCell
    using ..TightBindingToolkit.BZone:BZ

    using LinearAlgebra

    @doc """
    ```julia
    FillHamiltonian(uc::UnitCell, k::Vector{Float64} ; SpinMatrices::Vector{Matrix{ComplexF64}} = SpinMats((uc.localDim-1)//2)) :: Matrix{ComplexF64}
    ```
    Returns the Hamiltonian at momentum point `k`, corresponding to the bonds present in `UnitCell`.

    """
    function FillHamiltonian(uc::UnitCell, k::Vector{Float64} ; SpinMatrices::Vector{Matrix{ComplexF64}} = SpinMats((uc.localDim-1)//2)) :: Matrix{ComplexF64}
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


    @doc raw"""
    ```julia
    FullHamiltonian(uc::UnitCell, bz::BZ) --> Matrix{Matrix{ComplexF64}}
    ```
    Returns the full Hamiltonian at all momentum points in `BZ`, corresponding to the bonds present in `UnitCell`.

    """
    function FullHamiltonian(uc::UnitCell, bz::BZ) :: Matrix{Matrix{ComplexF64}}
        SpinMatrices    =   SpinMats((uc.localDim-1)//2)
        return FillHamiltonian.(Ref(uc), bz.ks ; SpinMatrices = SpinMatrices)
    end


    @doc """
    `Hamiltonian` is a data type representing a general momentum-space Hamiltonian corresponding to the given `UnitCell` and `BZ`.

    # Attributes
    - `H           :: Matrix{Matrix{ComplexF64}}`: A matrix (corresponding to the matrix of k-points in `BZ`) of Hamiltonian matrices.
    - `bands       :: Matrix{Vector{Float64}}`: A matrix (corresponding to the matrix of k-points in `BZ`) of band spectrums.
    - `states      :: Matrix{Matrix{ComplexF64}}`: A matrix (corresponding to the matrix of k-points in `BZ`) of band wavefunctions.
    - `bandwidth   :: Tuple{Float64, Float64}` : the tuple of minimum and maximum energies in the band structure.

    Initialize this structure using
    ```julia
    Hamiltonian(uc::UnitCell, bz::BZ)
    ```
    """
    mutable struct Hamiltonian
        H           ::  Matrix{Matrix{ComplexF64}}
        bands       ::  Matrix{Vector{Float64}}
        states      ::  Matrix{Matrix{ComplexF64}}
        bandwidth   ::  Tuple{Float64, Float64}

        Hamiltonian(uc::UnitCell, bz::BZ) = new{}(FullHamiltonian(uc, bz), Array{Vector{Float64}}(undef, 0, 0), Array{Matrix{ComplexF64}}(undef, 0, 0), (0.0, 0.0))
    end


    @doc """
    ```julia
    DiagonalizeHamiltonian!(Ham::Hamiltonian)
    ```
    Diagonalize the `Hamiltonian` at all momentum points in the `BZ`.

    """
    function DiagonalizeHamiltonian!(Ham::Hamiltonian)
        sols            =   eigen.(Hermitian.(Ham.H))
        Ham.bands       =   getfield.(sols, :values)
        Ham.states      =   getfield.(sols, :vectors)
        Ham.bandwidth   =   extrema(reduce(vcat, Ham.bands))
        println("Hamiltonian Diagonalized")
    end


    @doc """
    ```julia
    DOS(Omega::Float64, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Float64
    DOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Vector{Float64}
    ```
    Calculate the Density of State correspondingto the given energies in `Omegas`, for the lowest bands upto `till_band`.
    The calculation is done at a finite `spread` of the delta-function sum. 
    """
    function DOS(Omega::Float64, energies::Vector{Float64}; spread::Float64=1e-3) :: Float64

        dos         =   @. 1 / ((Omega - energies) + im * spread)
        return sum(imag(dos))
    end

    function DOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) :: Vector{Float64}

        energies    =   reduce(vcat, Ham.bands)
        n_bands     =   length(Ham.bands[1, 1])
        energies    =   energies[filter(i -> (i-1) % n_bands + 1 <= till_band, 1:length(energies))]

        dos     =   DOS.(Omegas, Ref(energies); spread=spread)

        dOmega  =   Omegas[2:end] .- Omegas[1:end-1]
        norm    =   sum(dos[1:end-1] .* dOmega)

        dos     .=  dos ./ norm
        return dos
    end
    
end