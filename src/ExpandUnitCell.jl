module ExpandUCell
    export ChangePrimitives!, ExpandUnitCell, ExpandBonds!

    using LinearAlgebra

    using ..TightBindingToolkit.Useful: Meshgrid
    using ..TightBindingToolkit.UCell: Bond, UnitCell, GetRealSpacePositions, AddBasisSite!, IsSameBond, FlipBond
    using ..TightBindingToolkit.DesignUCell: AddAnisotropicBond!


@doc """
```julia
ChangePrimitives!(uc::UnitCell{T}, newPrimitives::Vector{Vector{Float64}})
```
Changes the pirmitive vectors of the given `UnitCell` assuming the sublattices stay the same. Changes the bonds accordingly.

"""
    function ChangePrimitives!(uc::UnitCell{T}, newPrimitives::Vector{Vector{Float64}} ; OffsetRange::Int64 = 2, accuracy::Int64 = 6) where {T}

        oldPrimitives 	=	deepcopy(uc.primitives)
        uc.primitives 	=	newPrimitives

        DistanceDict 	=	GetRealSpacePositions(uc ; OffsetRange = OffsetRange, accuracy = accuracy)
    
        newBonds 		=	Bond{T}[]

        for (sub, basis) in enumerate(uc.basis)

            bases 		=	getproperty.(uc.bonds, :base)
            bonds 	    =	uc.bonds[findall(==(sub), bases)]

            for bond in bonds
                position 	    =	uc.basis[bond.target] + sum(bond.offset .* oldPrimitives)
                target, offset 	=	DistanceDict[round.(position, digits = accuracy)]

                proposal        =   Bond(bond.base, target, offset, bond.mat, bond.dist, bond.label)
                if isempty(newBonds)
                    push!( newBonds , proposal )
                else
                    if sum(IsSameBond.( Ref(proposal) , newBonds ))==0
                        push!( newBonds , proposal )
                    end
                end
            end

            targets 		=	getproperty.(uc.bonds, :target)
            bonds 	        =	FlipBond.(uc.bonds[findall(==(sub), targets)])

            for bond in bonds
                position 	    =	uc.basis[bond.target] + sum(bond.offset .* oldPrimitives)
                target, offset 	=	DistanceDict[round.(position, digits = accuracy)]

                proposal    =   Bond(bond.base, target, offset, bond.mat, bond.dist, bond.label)
                if isempty(newBonds)
                    push!( newBonds , proposal )
                else
                    if sum(IsSameBond.( Ref(proposal) , newBonds ))==0
                        push!( newBonds , proposal )
                    end
                end
            end


        end
        @assert length(newBonds) == length(uc.bonds)
        uc.bonds 		=	newBonds
    end


@doc """
```julia
ExpandUnitCell(ucOG::UnitCell{T}, scaling::Vector{Int64} ; OffsetRange::Int64 = 2)
```
Returns a UnitCell which is an integer multiple of the given UnitCell (the primitives are scaled by the given scalings along each direction). The bonds are properly redefined amongst the new sublattices and new primitive vectors.

"""
    function ExpandUnitCell(ucOG::UnitCell{T}, scaling::Vector{Int64} ; OffsetRange::Int64 = 2, accuracy::Int64 = 6) :: UnitCell{T} where {T}

        asNew 		=	scaling .* ucOG.primitives
        ucNew 		=	UnitCell(asNew, ucOG.localDim, T)
        NewOffsets 	=	Meshgrid(scaling .- 1 ; starts = zeros(Int64, length(scaling)))

        for offset in NewOffsets
            for (iBasis, basis) in enumerate(ucOG.basis)
                b 	=	basis + sum(offset .* ucOG.primitives)
                AddBasisSite!(ucNew, b, ucOG.fields[iBasis], ucOG.OnSiteMats)
            end
        end

        DistanceDictNew 	=	GetRealSpacePositions(ucNew ; OffsetRange = OffsetRange, accuracy = accuracy)

        bases 		=	getproperty.(ucOG.bonds, :base)
        for sub in 1:length(ucNew.basis)
            oldSub 	=	 ((sub - 1) % length(ucOG.basis)) + 1
            bonds 	=	ucOG.bonds[findall(==(oldSub), bases)]

            for bond in bonds
                targetPosition 	=	ucNew.basis[sub] + ((ucOG.basis[bond.target] - ucOG.basis[bond.base] ) + sum(bond.offset .* ucOG.primitives))

                target, offset 	=	DistanceDictNew[round.(targetPosition, digits = accuracy)]
                AddAnisotropicBond!(ucNew, sub, target, offset, bond.mat, bond.dist, bond.label)
            end
        end

        return ucNew

    end


@doc """
```julia
ExpandBonds!(ucOG::UnitCell{T}, ucNew::UnitCell{T} ; OffsetRange::Int64 = 2)
```
An in-place version of `ExpandUnitCell` which works with the new primitive vectors given in `ucNew`.

"""
    function ExpandBonds!(ucOG::UnitCell{T}, ucNew::UnitCell{T} ; OffsetRange::Int64 = 2, accuracy::Int64 = 6) where {T}

        DistanceDictNew 	=	GetRealSpacePositions(ucNew ; OffsetRange = OffsetRange, accuracy = accuracy)

        bases 		=	getproperty.(ucOG.bonds, :base)
        
        for sub in 1:length(ucNew.basis)
            oldSub 	=	 ((sub - 1) % length(ucOG.basis)) + 1

            bonds 	=	ucOG.bonds[findall(==(oldSub), bases)]
            for bond in bonds
                targetPosition 	=	ucNew.basis[sub] + ((ucOG.basis[bond.target] - ucOG.basis[bond.base] ) + sum(bond.offset .* ucOG.primitives))

                target, offset 	=	DistanceDictNew[round.(targetPosition, digits = accuracy)]
                AddAnisotropicBond!(ucNew, sub, target, offset, bond.mat, bond.dist, bond.label)
            end
        end

    end



end