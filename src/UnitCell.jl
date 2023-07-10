module UCell
	export Bond , BondRank, IsSameBond , UnitCell , IsSameUnitCell , AddBasisSite! , GetDistance , GetRealSpacePositions, IsSameUnitCell

	using LinearAlgebra, Logging
	using ..TightBindingToolkit.Useful: GetAllOffsets

	@doc """
	`Bond{T}` is a data type representing a general bond on a lattice.

	# Attributes
	- `base::Int64`: sub-lattice of the initial site on the bond.
	- `target::Int64`: sub-lattice of the final site on the bond.
	- `offset::Vector{Int64}`: the difference of the unit cells in which these sublattices belong to, in units of the lattice basis vectors.
	- `mat::Array{ComplexF64, T}`: an array describing this bond --> can be hopping for partons, or spin-exchange for spins. The rank of the array, `T`, is determined by the specific  model being looked at. Eg. rank = 2 for a free parton Hamiltonian.
	- `dist::Float64`: the distance b/w the two sites = length of bond.
	- `label::String`: some string label to mark the bond type.

	"""
	mutable struct Bond{T}
		base	::  Int64
		target  ::  Int64
		offset 	::	Vector{Int64}
		mat 	::	Array{ComplexF64, T}
		dist 	::	Float64
		label 	::	String
	end


	@doc """
	Function to return rank of a bond or a collection of bonds.
	
	"""
	function BondRank(bond::Bond{T}) :: Int64 where{T}
		return T
	end
	
	function BondRank(bonds::Array{Bond{T}}) :: Int64 where{T}
		return T
	end


	@doc """
	Function to check if two bond objects are describing the same physical bond, just inverted! 
	"""
	function IsSameBond( Bond1::Bond , Bond2::Bond ) :: Bool
		return Bond1.base==Bond2.target && Bond1.target == Bond2.base && Bond1.offset == -Bond2.offset && Bond1.label==Bond2.label && BondRank(Bond1) == BondRank(Bond2)
	end


	@doc """
	`UnitCell{T}` is a data type representing a general unit cell of a lattice.

	# Attributes
	- `primitives  :: Vector{ Vector{ Float64 } }`: primitive bases vectors of the lattice.
	- `basis       :: Vector{ Vector{ Float64 } }`: positions of the sub-lattices.
	- `bonds       :: Vector{Bond{T}}`: the set of all bonds defining a lattice.
	- `fields      :: Vector{ Vector{Float64}}` : the fields oneach basis site.
	- `localDim    :: Int64`: Local Hilbert space dimension ( e.g. 3 for classical spins, 2 for spin-1/2 electrons ).
	- `BC          :: Vector{ ComplexF64 }`: boundary conditions, in the form of [e^{ιθ_i}]. e.g. θ=0 for PBC, θ=π for APBC, and so on.

	Initialize this structure using 
	```julia
	UnitCell( as::Vector{Vector{Float64}} , localDim::Int64)
	UnitCell( as::Vector{Vector{Float64}} , localDim::Int64, rank::Int64)
	```
	"""
	mutable struct UnitCell{T}
		primitives  :: Vector{ Vector{ Float64 } } # Lattice basis vectors
		basis       :: Vector{ Vector{ Float64 } } # Sublattice positions
		
		"""
		The 'bonds' attribute contains ( base, target, offset, matrix , distance, label )
		"""
		bonds       :: Vector{Bond{T}}
		fields      :: Vector{ Vector{Float64}} # Zeeman fields
		localDim    :: Int64 # Local Hilbert space dimension ( 3 for classical spins, 2 for partons )
		BC          :: Vector{ ComplexF64 } # Boundary condition
	
		function UnitCell( as::Vector{Vector{Float64}} , localDim::Int64)
			@warn "Bond rank not passed when constructing UnitCell. Choosing default value of 2."
			return new{2}( as , Vector{Float64}[] , Bond{2}[] , Vector{Float64}[] , localDim , ones(ComplexF64, length(as)) )
		end
	
		function UnitCell( as::Vector{Vector{Float64}} , localDim::Int64, rank::Int64)
			return new{rank}( as , Vector{Float64}[] , Bond{rank}[] , Vector{Float64}[] , localDim , ones(ComplexF64, length(as)) )
		end
	
	end 


	@doc """
	```julia
	AddBasisSite!( uc::UnitCell , position::Vector{Float64} )
	AddBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )
	```
	Add a sublattice to the `UnitCell`  at the given real-space position, with an on-site `field`.

	"""
	function AddBasisSite!( uc::UnitCell , position::Vector{Float64} )
		push!( uc.basis , position )
		d 	=	(uc.localDim^2)
		@warn "No On-Site field passed to basis site. Choosing default value of $(zeros(Float64, d)) given the local Hilbert space of UnitCell."
		push!( uc.fields , zeros(Float64 , d) )
	end

	function AddBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )
		push!( uc.basis , position )
		push!( uc.fields , field )
	end


	@doc """
	```julia
	GetDistance(uc::UnitCell, base::Int64, target::Int64, offset::Vector{Int64}) --> Float64
	```
	get the distance between site at position (0, `base`) and (R, `target`), where R = `offset`, when written in units of the unit cell primitive vectors.

	"""
	function GetDistance(uc::UnitCell, base::Int64, target::Int64, offset::Vector{Int64}) :: Float64
		return norm( sum( offset.*uc.primitives ) + (uc.basis[target] - uc.basis[base] ) )
	end


	@doc """
	```julia
	GetRealSpacePositions(uc::UnitCell{T} ; OffsetRange::Int64 = 2) --> Dict
	```
	Returns a dictionary whose keys are vectors in the cartesian coordinates (rounded off to `accuracy` digits), with values giving the corresponding sublattice and Unit Cell coordinate of a lattice site at that position.

	"""
	function GetRealSpacePositions(uc::UnitCell{T} ; OffsetRange::Int64 = 2, accuracy::Int64 = 6) :: Dict{Vector{Float64}, Tuple{Int64, Vector{Int64}}} where {T}

		offsets 	=	GetAllOffsets(OffsetRange, length(uc.primitives))
		DistanceDict  	=	Dict{Vector{Float64}, Tuple{Int64, Vector{Int64}}}()
	
		for offset in offsets
			for sub in 1:length(uc.basis)
	
				position 	=	uc.basis[sub] + sum(offset .* uc.primitives)
				DistanceDict[round.(position, digits = accuracy)] 	=	(sub, offset)
	
			end
		end
	
		return DistanceDict
	end


	@doc """
	```julia
	IsSameUnitCell(uc_1::UnitCell, uc_2::UnitCell) --> Bool
	```
	Function to check if two unit cell live on the same underlying lattice or not, i.e. have the same `primitives`, same `sublattices`, and same `localDim`.

	"""
	function IsSameUnitCell(uc_1::UnitCell, uc_2::UnitCell) :: Bool
		
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


end
