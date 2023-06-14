module Parameters

    using ..TightBindingToolkit.UCell:Bond, UnitCell, getAllOffsets, isSameBond, RemoveBonds!
    import ..TightBindingToolkit.UCell:addAnisotropicBond!, addIsotropicBonds!

    export Param, addAnisotropicBond!, addIsotropicBonds!, CreateUnitCell!, ModifyUnitCell!, GetParams

    using LinearAlgebra

    @doc """
	`Param` is a data type representing a general tight-binding parameter, which can span multiple bonds.
	
	# Attributes
	- `value        ::  Vector{ Float64 }`: the strength of the parameter (or even the full history of it if its changed).
	- `unitBonds    ::  Vector{ Bond }`: All the bonds this parameter lives on. These bonds are supposed to have "unit" strength, and ultimately get scaled by the `value` when making the `UnitCell`.
	- `label        ::  String`: some string label to mark the parameter.  
	- `dist         ::  Float64`: the distance of the bonds the parameter lives on.
	
	Initialize this structure using 
	```julia
	Param( value::Float64 )
	```
	
    """
    mutable struct Param  
        value       ::  Vector{Float64}
        unitBonds   ::  Vector{Bond} 
        label       ::  String
        dist        ::  Float64

        Param( value::Float64 ) = new{}( [value] , Bond[], "", -1.0 )
    end 


    @doc """
	```julia
	addAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
	addAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
	```
	Add a bond with the given attributes to `param`.
	If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.

	"""
    function addAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
        
        @assert uc.localDim == 1

        if base <= length(uc.basis) && target <= length(uc.basis)
            if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
                push!( param.unitBonds , Bond( base , target , offset , [mat;;] , dist, label ) )
                if param.label==""
                    param.label     =   label
                    param.dist      =   dist
                else
                    @assert param.label == label && param.dist == dist
                end
            else
                println( "Issue with bond between " , base , " and " , target , " at distance " , dist )
            end
        else
            println("One or both of those basis sites have not been added to the UnitCell object.")
        end
    end

    function addAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
        
        @assert size(mat) == (uc.localDim, uc.localDim) "Intertaction matrix has the wrong dimension!"

        if base <= length(uc.basis) && target <= length(uc.basis)
            if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
                push!( param.unitBonds , Bond( base , target , offset , mat , dist, label ) )
                if param.label==""
                    param.label     =   label
                    param.dist      =   dist
                else
                    @assert param.label == label && param.dist == dist
                end
            else
                println( "Issue with bond between " , base , " and " , target , " at distance " , dist )
            end
        else
            println("One or both of those basis sites have not been added to the UnitCell object.")
        end
    end


    @doc """
	```julia
	addIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))
	addIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )
	```
	Add a set of "isotropic" bonds, which are the same for each pair of sites at the given distance. 
	If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.
	The input `checkOffsetRange` must be adjusted depending on the input distance. 
	The optional input `subs` is meant for isotropic bonds when only a subset of sublattices are involved.

    """
    function addIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))

        @assert uc.localDim == 1
        offsets 		=	getAllOffsets(checkOffsetRange, length(uc.primitives))    

        for i in subs
            for j in subs
                for offset in offsets
                    if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
                        proposal 	=	Bond(i, j, offset, [mats;;], dist, label)
                        if sum(isSameBond.( Ref(proposal) , param.unitBonds ))==0
                            push!( param.unitBonds , proposal )

                            if param.label==""
                                param.label     =   label
                                param.dist      =   dist
                            else
                                @assert param.label == label && param.dist == dist
                            end
                        end
                    end
                end

            end
        end
    end

    function addIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )

        @assert size(mats) == (uc.localDim, uc.localDim) "Intertaction matrix has the wrong dimension!"
        offsets 		=	getAllOffsets(checkOffsetRange, length(uc.primitives))    

        for i in subs
            for j in subs
                for offset in offsets
                    if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
                        proposal 	=	Bond(i, j, offset, mats, dist, label)
                        if sum(isSameBond.( Ref(proposal) , param.unitBonds ))==0
                            push!( param.unitBonds , proposal )

                            if param.label==""
                                param.label     =   label
                                param.dist      =   dist
                            else
                                @assert param.label == label && param.dist == dist
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
    function CreateUnitCell!(uc::UnitCell, param::Param , index::Int64=length(param.value))
        append!(uc.bonds, deepcopy(param.unitBonds))
        map(x -> x.mat = param.value[index] * x.mat, uc.bonds[findall(==(param.label), getproperty.(uc.bonds, :label))])

    end

    function CreateUnitCell!(uc::UnitCell, params::Vector{Param}, indices::Vector{Int64}=length.(getproperty.(params, :value)))
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
    function ModifyUnitCell!(uc::UnitCell, param::Param)

        RemoveBonds!(uc, param.label)
        CreateUnitCell!(uc, param)
    end

    function ModifyUnitCell!(uc::UnitCell, params::Vector{Param})
        
        ModifyUnitCell!.(Ref(uc), params)
    end


    @doc """
	```julia
    	GetParams(uc::UnitCell) --> Vector{Param}
    	```
	For legacy purposes. 
    	If you have a `UnitCell` built using the old technique of adding bonds directly, you can get a vector of `Param` using this function, corresponding to each unique bond type already present in `UnitCell`.

    """
    function GetParams(uc::UnitCell) :: Vector{Param}
        
        params  =   Param[]
        for bond in uc.bonds
            if bond.label in getproperty.(params, :label)
                index   =   findfirst(==(bond.label), getproperty.(params, :label))
                addAnisotropicBond!(params[index], uc, bond.base, bond.target, bond.offset, bond.mat, bond.dist, bond.label)
            else
                newParam    =   Param(1.0)
                push!(params, newParam)
                addAnisotropicBond!(params[end], uc, bond.base, bond.target, bond.offset, bond.mat, bond.dist, bond.label)
            end
        end

        return params

    end


end
