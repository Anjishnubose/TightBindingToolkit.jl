module UCell
	export Bond , isSameBond , UnitCell , isSameUnitCell , getDistance , addBasisSite! , addAnisotropicBond! , addIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!

	using LinearAlgebra

	@doc """
	`Bond{T<:Number}` is a data type representing a general bond on a lattice.

	# Attributes
	- `base::Int64`: sub-lattice of the initial site on the bond.
	- `base::Int64`: sub-lattice of the final site on the bond.
	- `offset::Vector{Int64}`: the difference of the unit cells in which these sublattices belong to, in units of the lattice basis vectors.
	- `mat::Matrix{T}`: the matrtix describing this bond --> can be hopping for partons, or spin-exchange for spins.
	- `dist::Float64`: the distance b/w the two sites = length of bond.
	- `label::String`: some string label to mark the bond type.

	"""
	mutable struct Bond{T<:Number}
		base	::  Int64
		target  ::  Int64
		offset 	::	Vector{Int64}
		mat 	::	Matrix{T}
		dist 	::	Float64
		label 	::	String
	end

	@doc """
	Function to check if two bond objects are describing the same physical bond, just inverted! 
	"""
	function isSameBond( Bond1::Bond , Bond2::Bond ) :: Bool
		if Bond1.base==Bond2.target && Bond1.target == Bond2.base && Bond1.offset == -Bond2.offset && Bond1.label==Bond2.label
			return true
		else
			return false
		end
	end


	@doc """
	`UnitCell` is a data type representing a general unit cell of a lattice.

	# Attributes
	- `primitives  :: Vector{ Vector{ Float64 } }`: primitive bases vectors of the lattice.
	- `basis       :: Vector{ Vector{ Float64 } }`: positions of the sub-lattices.
	- `bonds       :: Vector{Bond}`: the set of all bonds defining a lattice.
	- `fields      :: Vector{ Vector{Float64}}` : the fields oneach basis site.
	- `localDim    :: Int64`: Local Hilbert space dimension ( e.g. 3 for classical spins, 2 for spin-1/2 electrons ).
	- `BC          :: Vector{ ComplexF64 }`: boundary conditions, in the form of [e^{ιθ_i}]. e.g. θ=0 for PBC, θ=π for APBC, and so on.

	Initialize this structure using 
	```julia
	UnitCell( as::Vector{Vector{Float64}} , localDim::Int64)
	```
	"""
	mutable struct UnitCell
		primitives  :: Vector{ Vector{ Float64 } } # Lattice basis vectors
		basis       :: Vector{ Vector{ Float64 } } # Sublattice positions
		
		"""
		The 'bonds' attribute contains ( base, target, offset, matrix , distance, label )
		"""
		bonds       :: Vector{Bond}
		fields      :: Vector{ Vector{Float64}} # Zeeman fields
		localDim    :: Int64 # Local Hilbert space dimension ( 3 for classical spins, 2 for partons )
		BC          :: Vector{ ComplexF64 } # Boundary condition
		
		UnitCell( as::Vector{Vector{Float64}} , localDim::Int64) = new{}( as , Vector{Float64}[] , Bond[] , Vector{Float64}[] , localDim , ones(ComplexF64, length(as)) )

	end 

	@doc """
	```julia
	isSameUnitCell(uc_1::UnitCell, uc_2::UnitCell) --> Bool
	```
	Function to check if two unit cell live on the same underlying lattice or not, i.e. have the same `primitives`, same `sublattices`, and same `localDim`.

	"""
	function isSameUnitCell(uc_1::UnitCell, uc_2::UnitCell) :: Bool
		
		check 	=	true
		for attribute in [:primitives, :basis, :localDim]
			prop1 	=	getproperty(uc_1, attribute)
			prop2 	=	getproperty(uc_2, attribute)

			if length(prop1)!=length(prop2)
				check 	=	false
				return check
			else
				check 	=	check && isapprox(prop1, prop2, atol=1e-3, rtol=1e-3)
			end
		end

		return check
	end


	@doc """
	```julia
	getDistance(uc::UnitCell, base::Int64, target::Int64, offset::Vector{Int64}) --> Float64
	```
	get the distance between site at position (0, `base`) and (R, `target`), where R = `offset`, when written in units of the unit cell primitive vectors.

	"""
	function getDistance(uc::UnitCell, base::Int64, target::Int64, offset::Vector{Int64}) :: Float64
		return norm( sum( offset.*uc.primitives ) + (uc.basis[target] - uc.basis[base] ) )
	end

	@doc """
	```julia
	addBasisSite!( uc::UnitCell , position::Vector{Float64} )
	addBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )
	```
	Add a sublattice to the `UnitCell`  at the given real-space position, with an on-site `field`.

	"""
	function addBasisSite!( uc::UnitCell , position::Vector{Float64} )
		push!( uc.basis , position )
		push!( uc.fields , zeros(Float64 , 4) )
	end

	function addBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )
		push!( uc.basis , position )
		push!( uc.fields , field )
	end


	@doc """
	```julia
	getAllOffsets(OffsetRange::Int64, dim::Int64) --> Vector{Vector{Int64}}
	```
	Given a range, returns the set of all possible Bond offsets such that each element of the offset vector lies in [-`OffsetRange`, `OffsetRange`].

	"""
	function getAllOffsets(OffsetRange::Int64, dim::Int64) :: Vector{Vector{Int64}}
		if dim==1
			offsets 	=	collect([i] for i in OffsetRange:-1:-OffsetRange)
		elseif  dim==2
			offsets 	=	reshape([[i, j] for i in OffsetRange:-1:-OffsetRange, j in OffsetRange:-1:-OffsetRange], (2*OffsetRange+1)^2)
		elseif  dim==3
			offsets 	=	reshape([[i, j, k] for i in OffsetRange:-1:-OffsetRange, j in OffsetRange:-1:-OffsetRange, k in OffsetRange:-1:-OffsetRange], (2*OffsetRange+1)^3)
		else
			println("Does not work for dimensions = ", dim)
		end
		return offsets
	end


	@doc """
	```julia
	addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
	addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
	```
	Add a bond with the given attributes to `UnitCell`.
	If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.

	"""
	function addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
		
		@assert uc.localDim == 1

		if base <= length(uc.basis) && target <= length(uc.basis)
			if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
				push!( uc.bonds , Bond( base , target , offset , [mat;;] , dist, label ) )
			else
				println( "Issue with bond between " , base , " and " , target , " at distance " , dist )
			end
		else
			println("One or both of those basis sites have not been added to the UnitCell object.")
		end
	end

	function addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
		
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


	@doc raw"""
	```julia
	addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))
	addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )
	```
	Add a set of "isotropic" bonds, which are the same for each pair of sites at the given distance. 
	If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.
	The input `checkOffsetRange` must be adjusted depending on the input distance. 
	The optional input `subs` is meant for isotropic bonds when only a subset of sublattices are involved.

	"""
	function addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))

		@assert uc.localDim == 1
		offsets 		=	getAllOffsets(checkOffsetRange, length(uc.primitives))    

		for i in subs
			for j in subs
				for offset in offsets
					if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
						proposal 	=	Bond(i, j, offset, [mats;;], dist, label)
						if sum(isSameBond.( Ref(proposal) , uc.bonds ))==0
							push!( uc.bonds , proposal )
						end
					end
				end

			end
		end
	end

	function addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )

		@assert size(mats) == (uc.localDim, uc.localDim) "Intertaction matrix has the wrong dimension!"
		offsets 		=	getAllOffsets(checkOffsetRange, length(uc.primitives))    

		for i in subs
			for j in subs
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


	@doc """
	```julia
	ModifyBonds!(uc::UnitCell, dist::Float64, newMat::Matrix{<:Number})
	ModifyBonds!(uc::UnitCell, label::String, newMat::Matrix{<:Number})
	```
	Modify an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`, to the given bond matrix.

	"""
	function ModifyBonds!(uc::UnitCell, dist::Float64, newMat::Matrix{<:Number})
		distances 	=	getfield.(uc.bonds, :dist)
		@assert size(newMat)==size(uc.bonds[1].mat)  "New entry is incompatible with old bond matrix"
		map(x -> x.mat = newMat, uc.bonds[findall(≈(dist), distances)])
	end

	function ModifyBonds!(uc::UnitCell, label::String, newMat::Matrix{<:Number})
		labels 	=	getfield.(uc.bonds, :label)
		@assert size(newMat)==size(uc.bonds[1].mat)  "New entry is incompatible with old bond matrix"
		map(x -> x.mat = newMat, uc.bonds[findall(==(label), labels)])
	end


	@doc """
	```julia
	ScaleBonds!(uc::UnitCell, dist::Float64, scale::Number)
	ScaleBonds!(uc::UnitCell, label::String, scale::Number)
	```
	Scale the matrix of an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`, by the given scaling factor.

	"""
	function ScaleBonds!(uc::UnitCell, dist::Float64, scale::Number)
		distances 	=	getfield.(uc.bonds, :dist)
		map(x -> x.mat = scale * x.mat, uc.bonds[findall(≈(dist), distances)])
	end

	function ScaleBonds!(uc::UnitCell, label::String, scale::Number)
		labels 	=	getfield.(uc.bonds, :label)
		map(x -> x.mat = scale * x.mat, uc.bonds[findall(==(label), labels)])
	end


	@doc """
	```julia
	RemoveBonds!(uc::UnitCell, dist::Float64)
	ScaleBonds!(uc::UnitCell, label::String)
	```
	Remove an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`.

	"""
	function RemoveBonds!(uc::UnitCell, label::String )
		labels 	=	getfield.(uc.bonds, :label)
		deleteat!( uc.bonds , findall(==(label), labels) )
	end

	function RemoveBonds!(uc::UnitCell, dist::Float64 )
		labels 	=	getfield.(uc.bonds, :dist)
		deleteat!( uc.bonds , findall(≈(dist), labels) )
	end


	@doc """
	```julia
	ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})
	ModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})
	```
	Modify the on-site fields in the `UnitCell`, either one at a time, or all of them.

	"""
	function ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})
		uc.fields[site] 	=	newField
	end

	function ModifyFields!(uc::UnitCell, newField::Vector{Float64}, dim::Int64)
		@assert length(newField) == length(uc.basis)
		setindex!.(uc.fields, newField, Ref(dim))
	end

	function ModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})
		uc.fields 	=	newField
	end


	@doc """
	```julia
	ModifyIsotropicFields!(uc::UnitCell, newField::Vector{Float64})
	ModifyIsotropicFields!(uc::UnitCell, newField::Float64, dim::Int64)
	```
	Modify the on site field uniformly, on all sublattices. The optional argument `dim` is if you want to only modify one of the 4 elements of on-site fields (3 Zeeman and 1 chemical potential).
	
	"""
	function ModifyIsotropicFields!(uc::UnitCell, newField::Vector{Float64})
		uc.fields 	= 	repeat(newField, length(uc.basis))				
	end
	
	function ModifyIsotropicFields!(uc::UnitCell, newField::Float64, dim::Int64)
		map(x -> x[dim] = newField, uc.fields )			
	end


end
