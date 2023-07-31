module LatticeStruct

    export Lattice, FillSites!, FillBonds!, FillLattice!, GetBCPhase, ApplyBCToSite

    using LinearAlgebra

    using ..TightBindingToolkit.Useful: Meshgrid
    using ..TightBindingToolkit.UCell: Bond, UnitCell


@doc """
```julia
GetMaxCoordinationNumber(uc::UnitCell) --> Int64
```
Returns the maximum coordination number amongs all the sublattices of the `UnitCell`.

"""
    function GetMaxCoordinationNumber(uc::UnitCell) :: Int64

        coords  =   zeros(Int64, length(uc.basis))
        for bond in uc.bonds

            if bond.base == bond.target && bond.offset == zeros(Int64, length(uc.primitives))
                coords[bond.base]    +=  1
            else
                coords[bond.base]    +=  1
                coords[bond.target]  +=  1
            end   
            
        end

        return max(coords...)
    end

    
@doc """
`Lattice{T}` is a data type representing a general real-space lattice constructed out of a `UnitCell{T}`.

# Attributes
- `uc          ::  UnitCell{T}`: Unit Cell of the lattice.
- `size        ::  Vector{Int64}`: size of the lattice in units of the UnitCell `primitives`.
- `length      ::  Int64`: total number of sites.
- `sites       ::  Dict{ Tuple{Int64, Vector{Int64}}, Int64}` : a dictionary with key = (sublattice, offset), and the corresponding value being the site number on the lattice.
- `positions   ::  Vector{ Vector{Float64}}`: real-space positions of all the lattice sites.
- `fields      ::  Vector{ Vector{Float64}}`: the fields on all the lattice sites.
- `BondSites   ::  Matrix{Int64}`: a matrix with `lattice.length` rows such that `BondSites[s, i]` is the site number of the ith neighbour of site-s on the lattice.
- `BondDists   ::  Matrix{Float64}`: a matrix with `lattice.length` rows such that `BondDists[s, i]` is the distance to the ith neighbour of site-s on the lattice.
- `BondLabels  ::  Matrix{String}`: a matrix with `lattice.length` rows such that `BondLabels[s, i]` is the label of the bond to the ith neighbour of site-s on the lattice.
- `BondMats    ::  Matrix{Array{ComplexF64, T}}`: a matrix with `lattice.length` rows such that `BondMats[s, i]` is the `Array{ComplexF64, T}` of the bond connecting to the ith neighbour of site-s on the lattice.

Initialize this structure using 
```julia
Lattice( uc::UnitCell{T}, size::Vector{Int64} ; null_dist::Float64 = -1.0, null_label::String = "-")
```

where `null_dist` and `null_label` are used for bonds which are not allowed due to given boundary conditions, but still tracked in the code.
"""
    mutable struct Lattice{T}

        uc          ::  UnitCell{T}
        size        ::  Vector{Int64}
        length      ::  Int64
        """
        All sites and fields
        """
        sites       ::  Dict{ Tuple{Int64, Vector{Int64}}, Int64}
        positions   ::  Vector{ Vector{Float64}}
        fields      ::  Vector{ Vector{Float64}}
        """
        Bond information
        """
        BondSites   ::  Matrix{Int64}
        BondDists   ::  Matrix{Float64}
        BondLabels  ::  Matrix{String}
        BondMats    ::  Matrix{Array{ComplexF64, T}}

        function Lattice( uc::UnitCell{T}, size::Vector{Int64} ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}

            Ncoords         =   GetMaxCoordinationNumber(uc)
            L               =   length(uc.basis) * prod(size)

            return new{T}( uc, size, L, Dict{ Tuple{Int64, Vector{Int64}}, Int64}(), Vector{Float64}[], Vector{Float64}[], zeros( Int64, L, Ncoords), null_dist .* ones( Float64, L, Ncoords), repeat([null_label], L, Ncoords), Matrix{Array{ComplexF64, T}}( undef, L, Ncoords))
        end
    end


@doc """
```julia
FillSites!(lattice::Lattice{T})
```
Fills all the sites, positions, and fields of the lattice using information in the `UnitCell`.

"""
    function FillSites!(lattice::Lattice{T}) where {T}

        offsets     =   collect.(Meshgrid(lattice.size .- Ref(1) ; starts = zeros(Int64, length(lattice.size))))
        site        =   1

        for offset in offsets
            for (b, basis) in enumerate(lattice.uc.basis)

                position    =   sum(offset .* lattice.uc.primitives) + basis

                lattice.sites[(b, offset)]  =   site
                site        =   site + 1

                push!(lattice.positions , position)
                push!(lattice.fields, lattice.uc.fields[b])

            end
        end
    
    end


@doc """
```julia
ApplyBCToSite(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) -->  Vector{Int64}
```
Returns the effective real-space offset of `site` given a lattice size `L`, and boundary conditions `BC`.

"""
    function ApplyBCToSite(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) :: Vector{Int64}

        BCStrength  =   Bool.(abs.(BC))     ##### strength of the boundary condition refers to it being 1 (for a periodic system with arbitrary phase), or 0 for open systems.
        BCPhase     =   angle.(BC)

        LEffective  =   L .* BCStrength + .!(BCStrength)    ##### Leff_i = L if system is closed, or 1 if system is open along that dimension
        return mod.(site, LEffective) + (.!(BCStrength)) .* (site)  ##### returns site_i mod L_i if system is closed, or site_i if system is open along that dimension.
    end


@doc """
```julia
GetBCPhase(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) --> ComplexF64
```
Returns the effective phase of `site` given a lattice size `L`, and boundary conditions `BC`.

"""
    function GetBCPhase(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) :: ComplexF64

        BCPhase     =   angle.(BC)  ##### the phase, ϕ_i corresponding to the boundary condition along each dimension
        phase       =   sum((div.(site, L, Ref(RoundDown))) .* BCPhase)     ##### a vector of [θ_i] s.t. θ_i = n_i * ϕ_i, where n_i  = site_i ÷ L_i
        return exp(im * phase)

    end


@doc """
```julia
FillBonds!(lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-" ) where {T}
```
Fills all the bond information in the lattice using the bonds in `UnitCell`.

"""
    function FillBonds!(lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-" ) where {T}

        bases           =   getproperty.(lattice.uc.bonds, :base)
        targets         =   getproperty.(lattice.uc.bonds, :target)
        coord           =   ones(Int64, lattice.length)

        flippedIndices 	=	collect(T:-1:1)

        for site in keys(lattice.sites)

            sub, offset     =   site
            s               =   lattice.sites[site]
            pos             =   lattice.positions[s]

            OutgoingBonds   =   lattice.uc.bonds[findall(==(sub), bases)]

            for bond in OutgoingBonds

                targetOffset    =   ApplyBCToSite(offset + bond.offset, lattice.size, lattice.uc.BC)
                targetPhase     =   GetBCPhase(   offset + bond.offset, lattice.size, lattice.uc.BC)

                target          =   get(lattice.sites, (bond.target, targetOffset), 0)

                lattice.BondSites[ s, coord[s]]    =   target
                lattice.BondDists[ s, coord[s]]    =   !iszero(target) * bond.dist + iszero(target) * null_dist
                lattice.BondLabels[s, coord[s]]    =   iszero(target) ? null_label : bond.label
                lattice.BondMats[  s, coord[s]]    =   !iszero(target) * targetPhase * bond.mat
                
                coord[s]   =   coord[s] + 1

            end

            IncomingBonds   =   lattice.uc.bonds[findall(==(sub), targets)]

            for bond in IncomingBonds

                if bond.base == bond.target && bond.offset == zeros(Int64, length(lattice.size))

                    continue
                else

                    targetOffset    =   ApplyBCToSite(offset - bond.offset, lattice.size, lattice.uc.BC)
                    targetPhase     =   GetBCPhase(   offset - bond.offset, lattice.size, lattice.uc.BC)

                    target          =   get(lattice.sites, (bond.base, targetOffset), 0)

                    lattice.BondSites[ s, coord[s]]    =   target
                    lattice.BondDists[ s, coord[s]]    =   !iszero(target) * bond.dist + iszero(target) * null_dist
                    lattice.BondLabels[s, coord[s]]    =   iszero(target) ? null_label : bond.label
                    lattice.BondMats[  s, coord[s]]    =   !iszero(target) * targetPhase * collect(conj.(permutedims(bond.mat, flippedIndices)))
                    
                    coord[s]   =   coord[s] + 1

                end
            end
        end
    end


@doc """
```julia
FillLattice!( lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}
```
Wrapper function to fill both site and bond information in the lattice.

"""
    function FillLattice!( lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}

        FillSites!(lattice)
        FillBonds!(lattice ; null_dist = null_dist, null_label = null_label)

    end





end