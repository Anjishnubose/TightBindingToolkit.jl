module Flux

    export GetBondPhase, CheckGaugeValidity, ModifyGauge!, GetStringPhase, InsertMonopolePair!

    using LinearAlgebra

    using ..TightBindingToolkit.Useful: SegmentIntersection
    using ..TightBindingToolkit.LatticeStruct: Lattice, FillLattice!


    function GetBondPhase(A::Function, r1::Vector{Float64}, r2::Vector{Float64} ; n::Int64 = 100, _kwargs::Dict = Dict()) :: ComplexF64

        ts  =   collect(range(0.0, 1.0, length = n + 1))
        rs  =   Ref(r1) .+ (Ref(r2 - r1) .* ts)
        As  =   A.(rs ; _kwargs...)

        dr  =   (r2 - r1) ./ n
        ProjectedAs =   dot.(As, Ref(dr))

        return sum(ProjectedAs) - 1/2 * (ProjectedAs[begin] + ProjectedAs[end])
    end


    function CheckGaugeValidity(lat::Lattice{T}, A::Function ; _kwargs::Dict = Dict(), accuracy::Int64 = 5) :: Bool where{T}

        latticePrimitives   =   lat.size .* lat.uc.primitives
        AShifts     =   A.(latticePrimitives ; _kwargs...)

        return prod(isinteger.(round.(AShifts / (2 * pi), digits = accuracy)))
    end


    function ModifyGauge!(lat::Lattice{T}, A::Function ; n::Int64 = 100, _kwargs::Dict = Dict()) where{T}

        positions   =   getindex.(getindex.(Ref(lat.positions), collect(1:lat.length)) , Ref(1)) 
        r1s         =   repeat(positions, 1, size(lat.BondSites, 2))
        r2s         =   getindex.(getindex.(Ref(lat.positions), lat.BondSites) , Ref(1))

        phases      =   exp.( im .* GetBondPhase.(Ref(A), r1s, r2s ; n = n, _kwargs=_kwargs))
        lat.BondMats=   lat.BondMats .* phases

    end


    function GetStringPhase(Monopoles::Vector{Vector{Float64}}, r1::Vector{Float64}, r2::Vector{Float64} ;  flux::Float64 = Float64(pi)) :: ComplexF64

        t = SegmentIntersection(Monopoles, [r1, r2])
        r = SegmentIntersection([r1, r2], Monopoles)

        if t>=0.0 && t<1.0 && r>=0.0 && r<1.0
            return exp(im * flux)
        else
            return 1.0
        end
    end


    function InsertMonopolePair!(lat::Lattice{T}, Monopoles::Vector{Vector{Float64}} ; flux::Float64 = Float64(pi) ) where{T}

        positions   =   getindex.(getindex.(Ref(lat.positions), collect(1:lat.length)) , Ref(1)) 
        r1s         =   repeat(positions, 1, size(lat.BondSites, 2))
        r2s         =   getindex.(getindex.(Ref(lat.positions), lat.BondSites) , Ref(1))

        dists       =   norm.(r2s .- r1s)
        check       =   isapprox.(dists, lat.BondDists, atol = 1e-6, rtol = 1e-6)

        phases      =   GetStringPhase.(Ref(Monopoles), r1s, r2s ; flux = flux)
        phases      =   phases .^ check
        lat.BondMats=   lat.BondMats .* phases

        return findall(==(false), isapprox.(phases, Ref(1.0)))

    end








































end