module LatHam
    export FillHamiltonian, LatticeHamiltonian, DiagonalizeHamiltonian!, Slater, SingleParticleFidelity, SlaterOverlap, GaugeFix!

    using ..TightBindingToolkit.LatticeStruct: Lattice
    import ..TightBindingToolkit.Hams: FillHamiltonian, DiagonalizeHamiltonian!

    using LinearAlgebra


@doc """
```julia
FillHamiltonian( lattice::Lattice{2} )
```
Fills the Hamiltonian matrix for a given `Lattice` object using the bonds stored in the `Lattice` object.
"""
    function FillHamiltonian( lattice::Lattice{2} )

        localDim    = lattice.uc.localDim
        Hamsize     = lattice.length * localDim
        Ham         = zeros( ComplexF64 , Hamsize , Hamsize )

        for site in 1:lattice.length
            ############# On-Site terms
            Ham[localDim * (site - 1) + 1 : localDim * (site), localDim * (site - 1) + 1 : localDim * (site)] += sum( lattice.fields[site] .* lattice.uc.OnSiteMats )

            ############# Bond terms
            for neighbour in 1:length(lattice.bondSites[site, :])

                if lattice.bondSites[site , neighbour] != 0
                    Ham[ localDim*(site-1) + 1 : localDim*(site) , localDim*(lattice.bondSites[site, neighbour]-1) + 1 : localDim*(lattice.bondSites[site, neighbour]) ] += lattice.bondMats[site, neighbour]
                    
                    if lattice.bondDists[site, neighbour] > 0.0
                        Ham[ localDim*(lattice.bondSites[site, neighbour]-1) + 1 : localDim*(lattice.bondSites[site, neighbour])  , localDim*(site-1) + 1 : localDim*(site)] += adjoint(lattice.bondMats[site, neighbour])
                    end
                end
            end 
        end 

        return Ham
    end

@doc """
`LatticeHamiltonian` is a data type representing a general real-space Hamiltonian constructed using a `Lattice` object.

# Attributes
- `H           ::      Matrix{ComplexF64}`: Hamiltonian matrix.
- `bands       ::      Vector{Float64}`: eigenvalues of the Hamiltonian.
- `states      ::      Matrix{ComplexF64}`: eigenvectors of the Hamiltonian.
- `is_BdG      ::      Bool`: whether the Hamiltonian is a BdG Hamiltonian or not.
- `bandwidth   ::      Tuple{Float64, Float64}`: minimum and maximum eigenvalues of the Hamiltonian.

Initialize this structure using 
```julia
LatticeHamiltonian( lattice::Lattice{2} )
```
"""
    mutable struct LatticeHamiltonian

        H           ::      Matrix{ComplexF64}
        bands       ::      Vector{Float64}
        states      ::      Matrix{ComplexF64}
        is_BdG      ::      Bool
        bandwidth   ::      Tuple{Float64, Float64}

        function LatticeHamiltonian(lattice::Lattice{2})

            localDim    = lattice.uc.localDim
            Hamsize     = lattice.length * localDim

            return new{}(FillHamiltonian(lattice), Vector{Float64}(undef, Hamsize), Matrix{ComplexF64}(undef, Hamsize, Hamsize), false, (0.0, 0.0))
        end
    end


@doc """
```julia
DiagonalizeHamiltonian!(H::LatticeHamiltonian)
```
Diagonalizes the Hamiltonian stored in `H` and stores the eigenvalues and eigenvectors in `H.bands` and `H.states` respectively.
Calculates the bandwidth too.
"""
    function DiagonalizeHamiltonian!(H::LatticeHamiltonian)

        sol        =   eigen(Hermitian(H.H))
        H.bands    =   sol.values
        H.states   =   sol.vectors
        H.bandwidth     =   (H.bands[begin], H.bands[end])

    end


    function GaugeFix!(H::LatticeHamiltonian, gauge::Tuple{Int64, Float64} = (1, 0.0))

        site, DesiredPhase = gauge[1], gauge[2]

        NormCheck = (abs2.(H.states[site:end, :]) .> 0.0)
        NonTrivialNormSites = findfirst.(==(1), eachcol(NormCheck)) .+ (site - 1)     

        phases = angle.(getindex.(Ref(H.states), NonTrivialNormSites, 1:length(H.bands)))
        phaseShift = exp.(im .* (DesiredPhase .- phases))

        H.states = H.states .* phaseShift'

        return (NonTrivialNormSites, phaseShift)

    end


    function SingleParticleFidelity(H1::LatticeHamiltonian, H2::LatticeHamiltonian, states::Vector{Int64} = collect(1:length(H1.bands)) ; fixGauge::Bool = true, gauge::Tuple{Int64, Float64} = (1, 0.0)) :: Matrix{ComplexF64}

        @assert size(H1.H) == size(H2.H) "The two hamiltonians must be the same size!"

        if fixGauge
            GaugeFix!(H1, gauge)
            GaugeFix!(H2, gauge)
        end

        fidelity = adjoint(H1.states[:, states]) * H2.states[:, states]

        return fidelity
    end


    function Slater(H::LatticeHamiltonian, positions::Vector{Int64}, states::Vector{Int64} = collect(1:length(positions))) :: ComplexF64
        
        @assert length(positions) == length(states) "Number of particles must be equal to the number of positions to be filled!"

        SlaterDet = H.states[positions, states]
        SlaterDet = det(SlaterDet)

        return SlaterDet
    end


    function SlaterOverlap(H1::LatticeHamiltonian, H2::LatticeHamiltonian, states::Vector{Int64} = collect(1:Int64(length(H1.bands)//2)); fixGauge::Bool = true, gauge::Tuple{Int64, Float64} = (1, 0.0)) :: ComplexF64

        return det((SingleParticleFidelity(H1, H2, states ; fixGauge = fixGauge, gauge = gauge)))
    end

end