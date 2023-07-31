module DesignLattice

    export CreateLattice, ModifyLattice!, ScaleLatticeBonds!, ScaleLattice!, RemoveLatticeBonds!

    using LinearAlgebra, Logging

    using ..TightBindingToolkit.UCell: Bond, UnitCell
    using ..TightBindingToolkit.Parameters: Param, CreateUnitCell!, ModifyUnitCell!
    using ..TightBindingToolkit.LatticeStruct: Lattice, FillLattice!


@doc """
```julia
CreateLattice(uc::UnitCell{T}, param::Param{T}, size::Vector{Int64} , index::Int64=length(param.value) ; null_dist::Float64 = -1.0, null_label::String = "-") --> Lattice{T} where {T}
CreateLattice(uc::UnitCell{T}, params::Vector{Param{T}}, size::Vector{Int64} , indices::Vector{Int64}=length.(getproperty.(params, :value)) ; null_dist::Float64 = -1.0, null_label::String = "-") :: Lattice{T} where {T}
```
Creates a lattice using `Param` objetcs given in `params`.

"""
    function CreateLattice(uc::UnitCell{T}, param::Param{T}, size::Vector{Int64} ; index::Int64=length(param.value), null_dist::Float64 = -1.0, null_label::String = "-") :: Lattice{T} where {T}

        CreateUnitCell!(uc, param, index)
        lattice     =   Lattice(uc, size ; null_dist = null_dist, null_label = null_label)
        FillLattice!(lattice)

        return lattice
    end

    function CreateLattice(uc::UnitCell{T}, params::Vector{Param{T}}, size::Vector{Int64} ; null_dist::Float64 = -1.0, null_label::String = "-") :: Lattice{T} where {T}

        CreateUnitCell!(uc, params)
        lattice     =   Lattice(uc, size ; null_dist = null_dist, null_label = null_label)
        FillLattice!(lattice)

        return lattice
    end
    
    
@doc """
```julia
ModifyLattice!(lattice::Lattice{T}, param::Param{T}) where {T}
ModifyLattice!(lattice::Lattice{T}, params::Vector{Param{T}}) where {T}
```
Modifies a lattice when the `Param` objects given in `params` are modified.

"""
    function ModifyLattice!(lattice::Lattice{T}, param::Param{T}) where {T}

        ModifyUnitCell!(lattice.uc, param)
        FillBonds!(lattice)
    end

    function ModifyLattice!(lattice::Lattice{T}, params::Vector{Param{T}} ) where {T}

        ModifyUnitCell!(lattice.uc, params)
        FillBonds!(lattice)
    end


@doc """
```julia
ScaleLatticeBonds!(lattice::Lattice{T} , label::String , scaling::Float64) where {T}
```
Scales a lattice bond with the given `label` and by the given `scaling`.

"""
    function ScaleLatticeBonds!(lattice::Lattice{T} , label::String , scaling::Float64) where {T}
        @assert !isnan(scaling) && !isinf(scaling) "scaling is NaN or Inf."
        lattice.BondMats[findall( ==(label) , lattice.BondLabels )] .= scaling .* lattice.BondMats[findall( ==(label) , lattice.BondLabels )]
    
    end


@doc """
```julia
ScaleLattice!(lattice::Lattice{T}, param::Param{T}) where {T}
ScaleLattice!(lattice::Lattice{T}, params::Vector{Param{T}}) where {T}
```
Scales a lattice bond assuming that the `Param` objects got their strengths modified.

"""
    function ScaleLattice!(lattice::Lattice{T}, param::Param{T}) where {T}
        scaling     =   param.value[end] / param.value[end - 1]
        ScaleLatticeBonds!(lattice, param.label, scaling )

    end

    function ScaleLattice!(lattice::Lattice{T}, params::Vector{Param{T}}) where {T}

        ScaleLattice!.(Ref(lattice), params)
    end


@doc """
```julia
RemoveLatticeBonds!(lattice::Lattice{T} , label::String ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}
```
Removes a lattice bond with the given `label`.

"""
    function RemoveLatticeBonds!(lattice::Lattice{T} , label::String ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}

        indices     =   findall( ==(label) , lattice.BondLabels )

        lattice.BondMats[indices]   .=  0.0 .* lattice.BondMats[indices]
        lattice.BondLabels[indices] .=  Ref(null_label)
        lattice.BondDists[indices]  .=  Ref(null_dist)    
    end





end