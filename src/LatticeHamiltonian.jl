module LatHam
    export FillHamiltonian, LatticeHamiltonian, DiagonalizeHamiltonian!

    using ..TightBindingToolkit.LatticeStruct: Lattice
    import ..TightBindingToolkit.Hams: FillHamiltonian, DiagonalizeHamiltonian!

    using LinearAlgebra


    function FillHamiltonian( lattice::Lattice{2} )

        localDim    = lattice.uc.localDim
        Hamsize     = lattice.length * localDim
        Ham         = zeros( ComplexF64 , Hamsize , Hamsize )

        for site in 1:lattice.length
            ############# On-Site terms
            Ham[localDim * (site - 1) + 1 : localDim * (site), localDim * (site - 1) + 1 : localDim * (site)] += sum( lattice.fields[site] .* lattice.uc.OnSiteMats )

            ############# Bond terms
            for neighbour in 1:length(lattice.BondSites[site, :])

                if lattice.BondSites[site , neighbour] != 0
                    Ham[ localDim*(site-1) + 1 : localDim*(site) , localDim*(lattice.BondSites[site, neighbour]-1) + 1 : localDim*(lattice.BondSites[site, neighbour]) ] += lattice.BondMats[site, neighbour]
                    
                    if lattice.BondDists[site, neighbour] > 0.0
                        Ham[ localDim*(lattice.BondSites[site, neighbour]-1) + 1 : localDim*(lattice.BondSites[site, neighbour])  , localDim*(site-1) + 1 : localDim*(site)] += adjoint(lattice.BondMats[site, neighbour])
                    end
                end
            end 
        end 

        return Ham
    end


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


    function DiagonalizeHamiltonian!(H::LatticeHamiltonian)

        sol        =   eigen(Hermitian(H.H))
        H.bands    =   sol.values
        H.states   =   sol.vectors
        H.bandwidth     =   (H.bands[begin], H.bands[end])

    end


















































end