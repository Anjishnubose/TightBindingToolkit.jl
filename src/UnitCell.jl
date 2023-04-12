using LinearAlgebra
using MKL

Id3 = Matrix(1.0I,3,3)
############## Pauli matrices ############################################
Pauli        =   zeros(ComplexF64, (4, 2, 2))
Pauli[1,:,:] =   [0.0+0.0*im 1.0+0.0*im; 1.0+0.0*im 0.0+0.0*im]
Pauli[2,:,:] =   [0.0+0.0*im 0.0-1.0*im; 0.0+1.0*im 0.0+0.0*im]
Pauli[3,:,:] =   [1.0+0.0*im 0.0+0.0*im; 0.0+0.0*im -1.0+0.0*im]
Pauli[4,:,:] =   [1.0+0.0*im 0.0+0.0*im; 0.0+0.0*im 1.0+0.0*im]

PauliVec     =   [ Pauli[1,:,:] , Pauli[2,:,:] , Pauli[3,:,:] ]
SpinVec      =   0.5 .* [ Pauli[1,:,:] , Pauli[2,:,:] , Pauli[3,:,:] ]

"""
Bond stucture contains the following attributes to classify a bond
	base	:	sublattice of the initial site
	target	:	sublattice of the final site
	offset 	:	the distance of the unit cells in which these sublattices belong to, in units of the lattice basis vectors
	mat 	:	the matrtix describing this bond --> can be hopping for partons, or spin-exchange for spins.
	dist 	:	the distance b/w the two sites = length of bond
	label 	:	some string label to mark the bond "type"  	
"""
mutable struct Bond
    base	::  Int64
    target  ::  Int64
	offset 	::	Vector{Int64}
	mat 	::	Matrix{}
	dist 	::	Float64
	label 	::	String
end

"""
Function to check if two bond objects are describing the same physical bond, just inverted! 
"""
function isSameBond( Bond1::Bond , Bond2::Bond )
	if Bond1.base==Bond2.target && Bond1.target == Bond2.base && Bond1.offset == -Bond2.offset
		@assert Bond1.mat ≈ Bond2.mat' "Matrices are inconsistent b/w flipped bonds"
		return true
	else
		return false
	end
end

mutable struct UnitCell
    primitives  :: Vector{ Vector{ Float64 } } # Lattice basis vectors
    basis       :: Vector{ Vector{ Float64 } } # Sublattice positions
    
    """
    The 'bonds' attribute contains ( base, target, offset, matrix , distance, label )
    """
    bonds       :: Vector{Bond}
    fields      :: Vector{ Vector{Float64}} # Weiss fields for spinons, Zeeman field for Spins
    localDim    :: Int # Local Hilbert space dimension ( 3 for classical spins, 2 for partons )
    BC          :: Vector{ Int } # Boundary condition

    UnitCell( a1 , localDim ) = new{}( [ a1 ] , [] , [] , [] , localDim , [] )
    UnitCell( a1 , a2 , localDim ) = new{}( [ a1 , a2 ] , [] , [] , [] , localDim , [] )
    UnitCell( a1 , a2 , a3 , localDim ) = new{}( [ a1 , a2 , a3 ] , [] , [] , [] , localDim , [] )

end 

"""
Function to calculate distance b/w sites (0, base), and (offset, target)
"""
function getDistance(uc::UnitCell, base, target, offset)
	return norm( sum( offset.*uc.primitives ) + (uc.basis[target] - uc.basis[base] ) )
end

function addBasisSite!( uc::UnitCell , position::Vector{Float64} )
    push!( uc.basis , position )
    push!( uc.fields , zeros(Float64 , 3) )
end

function addBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )
    push!( uc.basis , position )
    push!( uc.fields , field )
end

"""
Use the function below for Anisotropic Bonds.
"""
function addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{} , dist::Float64, label::String )
	@assert size(mat) == (uc.localDim, uc.localDim) "Intertaction matrix has the wrong dimension!"
    if base <= length(uc.basis) && target <= length(uc.basis)
		if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
        	push!( uc.bonds , Bond( base , target , offset , mat , dist, label ) )
		else
            println( "Issue with bond between " , base , " and " , target , " at distance " , dist )
		end
    else
        println("One or both of those basis sites have not been added to the UnitCell object.")
    end
end

function getAllOffsets(OffsetRange::Int64, dim::Int64)
	if dim==1
		offsets 	=	collect(OffsetRange:-OffsetRange:-1)
	elseif  dim==2
		offsets 	=	reshape([[i, j] for i in OffsetRange:-OffsetRange:-1, j in OffsetRange:-OffsetRange:-1], (2*OffsetRange+1)^2)
	elseif  dim==3
		offsets 	=	reshape([[i, j, k] for i in OffsetRange:-OffsetRange:-1, j in OffsetRange:-OffsetRange:-1, k in OffsetRange:-OffsetRange:-1], (2*OffsetRange+1)^3)
	else
		println("Does not work for dimensions = ", dim)
	end
	return offsets
end


"""
Use the function below for Isotropic Bonds.
"""
function addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Matrix{} , label::String, checkOffsetRange::Int64 )

	offsets 		=	getAllOffsets(checkOffsetRange, length(uc.primitives))    

    for i in 1:length(uc.basis)
        for j in 1:length(uc.basis)
            for offset in offsets
                if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
                    proposal 	=	Bond(i, j, offset, mats, dist, label)
					if sum(isSameBond.( Ref(proposal) , uc.bonds ))==0
						push!( uc.bonds , proposal )
					end
                end
            end

        end
    end
end

"""
Function to modify specific bond matrices given a bond distance or a bond label
"""
function ModifyBonds!(uc::UnitCell, dist::Float64, newMat::Matrix{})
	distances 	=	getfield.(uc.bonds, :dist)
	map(x -> x.mat = newMat, uc.bonds[findall(≈(dist), distances)])
end

function ModifyBonds!(uc::UnitCell, label::String, newMat::Matrix{})
	labels 	=	getfield.(uc.bonds, :label)
	map(x -> x.mat = newMat, uc.bonds[findall(==(label), labels)])
end

function RemoveBonds!(uc::UnitCell, label::String )
	labels 	=	getfield.(uc.bonds, :label)
	deleteat!( uc.bonds , findall(==(label), labels) )
	# map(x -> x.mat = newMat, uc.bonds[findall(==(label), labels)])
end

function RemoveBonds!(uc::UnitCell, dist::Float64 )
	labels 	=	getfield.(uc.bonds, :dist)
	deleteat!( uc.bonds , findall(≈(dist), labels) )
	# map(x -> x.mat = newMat, uc.bonds[findall(==(label), labels)])
end

"""
Function to modify site fields
"""
function ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})
	uc.fields[site] 	=	newField
end

function ModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})
	uc.fields 	=	newField
end



