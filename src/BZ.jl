using LinearAlgebra
include("UnitCell.jl")

"""
Function to get reciprocal Lattice vectors
"""
function getRLVs( uc::UnitCell )
    if length(uc.primitives) == 2
        a1 = vcat( uc.primitives[1] , 0.0 )
        a2 = vcat( uc.primitives[2] , 0.0 )
        a3 = [ 0.0 , 0.0 , 1.0 ]
        V   = dot( a1 , cross( a2 , a3 ) )
        b1  = 2 * pi / V * ( cross( a2 , a3 ) )
        b2  = 2 * pi / V * ( cross( a3 , a1 ) )
        return [ b1[1:2] , b2[1:2] ]

    elseif length(uc.primitives) == 3
        V   = dot( uc.primitives[1] , cross( uc.primitives[2] , uc.primitives[3] ) )
        b1  = 2 * pi / V * ( cross( uc.primitives[2] , uc.primitives[3] ) )
        b2  = 2 * pi / V * ( cross( uc.primitives[3] , uc.primitives[1] ) )
        b3  = 2 * pi / V * ( cross( uc.primitives[1] , uc.primitives[2] ) )
        return [ b1 , b2 , b3 ]
    else
        print("This code only works for 2D and 3D lattices. ")
    end
end

mutable struct BZ
    basis	    ::  Vector{Vector{Float64}}
    gridSize    ::  Int64
    k1s         ::  Matrix{Float64}
    k2s         ::  Matrix{Float64}
	ks   	    ::	Matrix{Vector{Float64}}
    kExps       ::  Dict

    BZ(gridSize)    =   new{}([], gridSize, Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0),  Array{Vector{Float64}}(undef, 0, 0), Dict())
end

"""
Function to get the Brillouin Zone of size kSize ---> returns a matrix of kVecs = [[kx, ky] for kx, ky in the brillouin Zone]
If you want a matrix of kVecs, just reshapa(BZ(uc, N), N, N) to get a matrix where each row is a fixed kx and each column is a fixed ky.
"""
function fillBZ!(bz::BZ, uc::UnitCell, offsetRange::Int64=1)

    bz.basis=   getRLVs(uc)
    @assert length(bz.basis)==2 "Sorry, the code is only written for d=2 right now!"

    k1Inds  =   repeat(1:bz.gridSize, outer=bz.gridSize)
    k2Inds  =   repeat(1:bz.gridSize, inner=bz.gridSize)
    bz.k1s  =   reshape((2 .* k1Inds .- (bz.gridSize + 1)) ./ (2*bz.gridSize), bz.gridSize, bz.gridSize)
    bz.k2s  =   reshape((2 .* k2Inds .- (bz.gridSize + 1)) ./ (2*bz.gridSize), bz.gridSize, bz.gridSize)

    bz.ks   =   (bz.k1s .* Ref(bz.basis[1])) .+ (bz.k2s .* Ref(bz.basis[2]))
    bz.kExps=   kExp(uc, bz, offsetRange)

end

"""
Function to return a dictionary containing all possible exp(I k.δ), where δ = n1 * a1 + n2 * a2, offset = [n1, n2], and a1, a2 are the primitives
"""
function kExp(uc::UnitCell, bz::BZ, offsetRange::Int64)

    offsets     =   getAllOffsets(offsetRange, length(uc.primitives))
    kDict       =   Dict()
    for offset in offsets
        kDict[string(offset)]   =   exp.( im .* dot.(bz.ks, Ref(sum(offset .* uc.primitives))))
    end

    return kDict
end



