module Parameters
    export Param, AddAnisotropicBond!, AddIsotropicBonds!, CreateUnitCell!, ModifyUnitCell!, GetParams

    using ..TightBindingToolkit.Useful: GetAllOffsets
    using ..TightBindingToolkit.UCell: Bond, UnitCell, IsSameBond
    using ..TightBindingToolkit.DesignUCell: RemoveBonds!
    import ..TightBindingToolkit.DesignUCell: AddAnisotropicBond!, AddIsotropicBonds!

    using LinearAlgebra, Logging

@doc """
`Param{T}` is a data type representing a general tight-binding parameter, which can span multiple bonds.

# Attributes
- `value        ::  Vector{ Float64 }`: the strength of the parameter (or even the full history of it if its changed).
- `unitBonds    ::  Vector{ Bond{T} }`: All the bonds this parameter lives on. These bonds are supposed to have "unit" strength, and ultimately get scaled by the `value` when making the `UnitCell`.
- `label        ::  String`: some string label to mark the parameter.  
- `dist         ::  Float64`: the distance of the bonds the parameter lives on.

Initialize this structure using 
```julia
Param( value::Float64 )
Param( value::Float64 , rank::Int64 )
```

"""
    mutable struct Param{T}  
        value       ::  Vector{Float64}
        unitBonds   ::  Vector{Bond{T}} 
        label       ::  String
        dist        ::  Float64
    
        function Param( value::Float64 ) 
            @warn "Rank not passed to Param object. Choosing default value of 2!"
            return new{2}( [value] , Bond{2}[], "", -1.0 )
        end
    
        function Param( value::Float64, rank::Int64 ) 
            return new{rank}( [value] , Bond{rank}[], "", -1.0 )
        end    
    end 


@doc """
```julia
AddAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
AddAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
```
Add a bond with the given attributes to `param`.
If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.

"""
    function AddAnisotropicBond!( param::Param{T}, uc::UnitCell{T} , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String ) where {T}
	
        @assert uc.localDim == 1
        dims 	=	repeat([uc.localDim], T)
    
        if base <= length(uc.basis) && target <= length(uc.basis)
            if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
                push!( param.unitBonds , Bond( base , target , offset , ComplexF64.(reshape([mat], dims...)) , dist, label ) )

                if param.label==""
                    param.label     =   label
                    param.dist      =   dist
                else
                    @assert param.label == label && param.dist == dist "Inconsistent label or distance given"
                end
            else
                @warn "Issue with bond between $(base) and $(target) at distance $(dist)"
            end
        else
            @warn "One or both of those basis sites, [$(base), $(target)] have not been added to the UnitCell object."
        end
    end
    
    function AddAnisotropicBond!( param::Param{T}, uc::UnitCell{T} , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Array{<:Number, T} , dist::Float64, label::String ) where {T}

        dims 	=	repeat([uc.localDim], T)
        @assert size(mat) == Tuple(dims) "Interaction matrix has the wrong dimension!"
    
        if base <= length(uc.basis) && target <= length(uc.basis)
            if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
                push!( param.unitBonds , Bond( base , target , offset , ComplexF64.(mat) , dist, label ) )

                if param.label==""
                    param.label     =   label
                    param.dist      =   dist
                else
                    @assert param.label == label && param.dist == dist "Inconsistent label or distance given"
                end
            else
                @warn "Issue with bond between $(base) and $(target) at distance $(dist)"
            end
        else
            @warn "One or both of those basis sites, [$(base), $(target)] have not been added to the UnitCell object."
        end
    end


@doc """
```julia
AddIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))
AddIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )
```
Add a set of "isotropic" bonds, which are the same for each pair of sites at the given distance. 
If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.
The input `checkOffsetRange` must be adjusted depending on the input distance. 
The optional input `subs` is meant for isotropic bonds when only a subset of sublattices are involved.

"""
    function AddIsotropicBonds!( param::Param{T}, uc::UnitCell{T} , dist::Float64 , mat::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis))) where {T}

        @assert uc.localDim == 1
        dims 	=	repeat([uc.localDim], T)
        offsets =	GetAllOffsets(checkOffsetRange, length(uc.primitives))    
    
        for i in subs
            for j in subs
                for offset in offsets

                    if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
                        proposal 	=	Bond(i, j, offset, ComplexF64.(reshape([mat], dims...)), dist, label)

                        if sum(IsSameBond.( Ref(proposal) , param.unitBonds ))==0
                            push!( param.unitBonds , proposal )
    
                            if param.label==""
                                param.label     =   label
                                param.dist      =   dist
                            else
                                @assert param.label == label && param.dist == dist "Inconsistent label or distance given"
                            end
                        end
                    end
                end
    
            end
        end
    end
    
    function AddIsotropicBonds!( param::Param{T}, uc::UnitCell{T} , dist::Float64 , mat::Array{<:Number, T} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) ) where {T}
    
        dims 	=	repeat([uc.localDim], T)
        @assert size(mat) == Tuple(dims) "Interaction matrix has the wrong dimension!"
        offsets 		=	GetAllOffsets(checkOffsetRange, length(uc.primitives))    
    
        for i in subs
            for j in subs
                for offset in offsets

                    if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
                        proposal 	=	Bond(i, j, offset, ComplexF64.(mat), dist, label)

                        if sum(IsSameBond.( Ref(proposal) , param.unitBonds ))==0
                            push!( param.unitBonds , proposal )
    
                            if param.label==""
                                param.label     =   label
                                param.dist      =   dist
                            else
                                @assert param.label == label && param.dist == dist "Inconsistent label or distance given"
                            end
                        end
                    end
                end
    
            end
        end
    end


@doc """
```julia
CreateUnitCell!(uc::UnitCell, param::Param , index::Int64=length(param.value))
CreateUnitCell!(uc::UnitCell, params::Vector{Param}, indices::Vector{Int64}=length.(getproperty.(params, :value)))
```
Add bonds corrsponding to a `param` to `UnitCell`, scaled with the `param.value[index]`. Also includes the broadcasted call.

"""
    function CreateUnitCell!(uc::UnitCell{T}, param::Param{T} , index::Int64=length(param.value)) where {T}
        bonds   =   deepcopy(param.unitBonds)
        map(x -> x.mat = param.value[index] * x.mat, bonds)
        append!(uc.bonds, bonds)
        ##### TODO: Check this works properly!!!

    end

    function CreateUnitCell!(uc::UnitCell{T}, params::Vector{Param{T}}, indices::Vector{Int64}=length.(getproperty.(params, :value))) where {T}
        @assert length(indices)==length(params) "Inconsistent input of params"
        CreateUnitCell!.(Ref(uc), params, indices)
    end


@doc """
```julia
ModifyUnitCell!(uc::UnitCell, param::Param)
ModifyUnitCell!(uc::UnitCell, params::Vector{Param})
```
Modify all bonds in `UnitCell` corresponding to given `param`, taking the latest value in `param.value`. 

"""
    function ModifyUnitCell!(uc::UnitCell{T}, param::Param{T}) where {T}

        RemoveBonds!(uc, param.label)
        CreateUnitCell!(uc, param)
    end

    function ModifyUnitCell!(uc::UnitCell{T}, params::Vector{Param{T}}) where {T}
        
        ModifyUnitCell!.(Ref(uc), params)
    end


@doc """
```julia
GetParams(uc::UnitCell) --> Vector{Param}
```
For legacy purposes. 
If you have a `UnitCell` built using the old technique of adding bonds directly, you can get a vector of `Param` using this function, corresponding to each unique bond type already present in `UnitCell`.

"""
    function GetParams(uc::UnitCell{T}) :: Vector{Param{T}} where {T}
        
        params  =   Param{T}[]
        for bond in uc.bonds
            labels      =   getproperty.(params, :label)

            if bond.label in labels
                index   =   findfirst(==(bond.label), labels)
                AddAnisotropicBond!(params[index], uc, bond.base, bond.target, bond.offset, bond.mat, bond.dist, bond.label)
            else
                newParam    =   Param(1.0, T)
                push!(params, newParam)
                AddAnisotropicBond!(params[end], uc, bond.base, bond.target, bond.offset, bond.mat, bond.dist, bond.label)
            end
        end

        return params

    end


end
