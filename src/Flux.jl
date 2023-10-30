module Flux

    export GetBondPhase, CheckGaugeValidity, ModifyGauge!, GetStringPhase, InsertMonopolePair!

    using LinearAlgebra

    using ..TightBindingToolkit.Useful: SegmentIntersection
    using ..TightBindingToolkit.LatticeStruct: Lattice, FillLattice!


@doc """
```julia
GetBondPhase(A::Function, r1::Vector{Float64}, r2::Vector{Float64} ; n::Int64 = 100, _kwargs::Dict = Dict()) --> ComplexF64
```
Returns the tight binding phase on a bond between two sites `r1` and `r2` on a lattice with a given gauge vector potential `A`. 
The phase is calculated by integrating the vector potential along the bond.
"""
    function GetBondPhase(A::Function, r1::Vector{Float64}, r2::Vector{Float64} ; n::Int64 = 100, _kwargs::Dict{Symbol, <:Any} = Dict{Symbol, Any}()) :: ComplexF64

        ts  =   collect(range(0.0, 1.0, length = n + 1))
        rs  =   Ref(r1) .+ (Ref(r2 - r1) .* ts)
        As  =   A.(rs ; _kwargs...)

        dr  =   (r2 - r1) ./ n
        ProjectedAs =   dot.(As, Ref(dr))

        return sum(ProjectedAs) - 1/2 * (ProjectedAs[begin] + ProjectedAs[end])
    end


@doc """
```julia
CheckGaugeValidity(lat::Lattice{T}, A::Function ; _kwargs::Dict = Dict(), accuracy::Int64 = 5) --> Bool 
```
Checks if the gauge vector potential `A` is valid for the lattice `lat` under periodic boundary conditions.
"""
    function CheckGaugeValidity(lat::Lattice{T}, A::Function ; _kwargs::Dict{Symbol, <:Any} = Dict{Symbol, Any}(), accuracy::Int64 = 5) :: Bool where{T}

        latticePrimitives   =   lat.size .* lat.uc.primitives
        AShifts     =   A.(latticePrimitives ; _kwargs...)

        return prod(isinteger.(round.(AShifts / (2 * pi), digits = accuracy)))
    end


@doc """
```julia
ModifyGauge!(lat::Lattice{T}, A::Function ; n::Int64 = 100, _kwargs::Dict = Dict()) where{T}
```
Modifies the lattice hoppings using the fixed gauge vector potential `A` on the lattice `lat` by multiplying the bond matrices with the corresponding phase factors.
"""
    function ModifyGauge!(lat::Lattice{T}, A::Function ; n::Int64 = 100, _kwargs::Dict{Symbol, <:Any} = Dict{Symbol, Any}()) where{T}

        positions   =   getindex.(Ref(lat.positions), collect(1:lat.length))
        r1s         =   repeat(positions, 1, size(lat.bondSites, 2))
        r2s         =   getindex.(Ref(lat.positions), lat.bondSites) .+ [sum(shift .* lat.size .* lat.uc.primitives) for shift in lat.bondShifts]

        phases      =   exp.( im .* GetBondPhase.(Ref(A), r1s, r2s ; n = n, _kwargs =_kwargs))
        lat.bondMats=   lat.bondMats .* phases

    end


@doc """
```julia
GetStringPhase(Monopoles::Vector{Vector{Float64}}, r1::Vector{Float64}, r2::Vector{Float64} ;  flux::Float64 = Float64(pi)) --> ComplexF64
```
Get the Dirac string between the given `Monopoles` phase cutting the bond on the lattice between sites `r1` and `r2`.
"""
    function GetStringPhase(Monopoles::Vector{Vector{Float64}}, r1::Vector{Float64}, r2::Vector{Float64} ;  flux::Float64 = Float64(pi)) :: ComplexF64

        t = SegmentIntersection(Monopoles, [r1, r2])
        r = SegmentIntersection([r1, r2], Monopoles)

        if t>=0.0 && t<1.0 && r>=0.0 && r<1.0
            return exp(im * flux)
        else
            return 1.0
        end
    end


@doc """
```julia
InsertMonopolePair!(lat::Lattice{T}, Monopoles::Vector{Vector{Float64}} ; flux::Float64 = Float64(pi) ) 
```
Inserts a pair of monopole <--> Anti-Monopole on the lattice `lat` at the given positions `Monopoles` with the given flux `flux` through a Dirac
"""
    function InsertMonopolePair!(lat::Lattice{T}, Monopoles::Vector{Vector{Float64}} ; flux::Float64 = Float64(pi) ) where{T}

        positions   =   getindex.(getindex.(Ref(lat.positions), collect(1:lat.length)) , Ref(1)) 
        r1s         =   repeat(positions, 1, size(lat.bondSites, 2))
        r2s         =   getindex.(getindex.(Ref(lat.positions), lat.bondSites) , Ref(1))

        dists       =   norm.(r2s .- r1s)
        check       =   isapprox.(dists, lat.bondDists, atol = 1e-6, rtol = 1e-6)

        phases      =   GetStringPhase.(Ref(Monopoles), r1s, r2s ; flux = flux)
        phases      =   phases .^ check
        lat.bondMats=   lat.bondMats .* phases

        return findall(==(false), isapprox.(phases, Ref(1.0)))

    end








































end