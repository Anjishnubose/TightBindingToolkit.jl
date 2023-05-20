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

###### Add a method to sample Qs outside the first BZ
mutable struct BZ
    basis	        ::  Vector{Vector{Float64}}
    gridSize        ::  Int64
    k1s             ::  Matrix{Float64}
    k2s             ::  Matrix{Float64}
	ks   	        ::	Matrix{Vector{Float64}}
    HighSymPoints   ::  Dict
    shift           ::  Vector{Float64}

    BZ(gridSize)    =   new{}([], gridSize, Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0),  Array{Vector{Float64}}(undef, 0, 0), Dict(), [])
end


function Monkhorst(ind::Int64, N::Int64) :: Float64
    return (2 * ind - (N + 1)) / (2 * N)
end

function VecAngle(a::Vector{Float64}, b::Vector{Float64})
    return acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

"""
Function to get the Brillouin Zone of size kSize ---> returns a matrix of kVecs = [[kx, ky] for kx, ky in the brillouin Zone]
If you want a matrix of kVecs, just reshapa(BZ(uc, N), N, N) to get a matrix where each row is a fixed kx and each column is a fixed ky.
"""
function fillBZ!(bz::BZ, uc::UnitCell, offsetRange::Int64=1 ; shift::Vector{Float64}=zeros(Float64, length(uc.primitives)))

    bz.basis    =   getRLVs(uc)
    @assert length(bz.basis)==2 "Sorry, the code is only written for d=2 right now!"
    bz.shift    =   shift

    k1Inds  =   repeat(1:bz.gridSize, outer=bz.gridSize)
    k2Inds  =   repeat(1:bz.gridSize, inner=bz.gridSize)
    ##### Monkhorst grid : Choose kSize = 6 * n + 3 to cover high symmetry points, when shift = [0.0, 0.0]
    bz.k1s  =   reshape(Monkhorst.(k1Inds, Ref(bz.gridSize)), bz.gridSize, bz.gridSize)
    bz.k2s  =   reshape(Monkhorst.(k2Inds, Ref(bz.gridSize)), bz.gridSize, bz.gridSize)

    bz.ks   =   (bz.k1s .* Ref(bz.basis[1])) .+ (bz.k2s .* Ref(bz.basis[2])) .+ Ref(shift)

    bz.HighSymPoints["Gamma"]   =   [0.0, 0.0]
    bz.HighSymPoints["M1"]      =  @. 0.5 * bz.basis[1] + 0.0 * bz.basis[2]
    bz.HighSymPoints["M2"]      =  @. 0.0 * bz.basis[1] + 0.5 * bz.basis[2]
    bz.HighSymPoints["M3"]      =  @. 0.5 * bz.basis[1] + 0.5 * bz.basis[2]

    if isapprox(VecAngle(bz.basis[1], bz.basis[2]), 2*pi/3, atol=1e-4, rtol=1e-4)
        bz.HighSymPoints["K1"]      =   @. (2/3) * bz.basis[1] + (1/3) * bz.basis[2]
        bz.HighSymPoints["K2"]      =   @. (1/3) * bz.basis[1] + (2/3) * bz.basis[2]
    elseif isapprox(VecAngle(bz.basis[1], bz.basis[2]), pi/3, atol=1e-4, rtol=1e-4)
        bz.HighSymPoints["K1"]      =   @. (1/3) * bz.basis[1] + (1/3) * bz.basis[2]
        bz.HighSymPoints["K2"]      =   @. (2/3) * bz.basis[1] + (2/3) * bz.basis[2]
    end

end


function ReduceQ(Q::Vector{Float64}, bz::BZ) :: Vector{Float64}
    U           =   reduce(hcat, bz.basis) ##### basis transformation from x, y -> b1, b2
    Q_reduced   =   inv(U) * (Q .- bz.shift)   ##### Q in the basis of b1 and b2
    Q_reduced   =   Q_reduced .- round.(Q_reduced)  ##### Q in the basis of b1 and b2, and shifted back to the first BZ
    Q_reduced   =   sum(Q_reduced .* bz.basis) .+ bz.shift ##### Q back in the x ,y basis, but shifted to the first BZ
    return Q_reduced
end


"""
function which takes in a Q and returns index [Q1, Q2] (if they exist) s,t Q = Q1 * b1 + Q2 * b2, where b1 and b2 are reciprocal lattice vectors
"""
function GetQIndex(Q::Vector{Float64}, bz::BZ ; nearest::Bool = false) :: Vector{Int64}

    Q_reduced   =   ReduceQ(Q, bz)

    if !nearest
        inds    =   findfirst(isapprox(Q_reduced, rtol=1e-5, atol=1e-5), bz.ks)
        @assert !isnothing(inds)  Q_reduced, "Given momentum does not exist on the BZ grid" 
    else
        inds    =   findmin(norm.(bz.ks .- Ref(Q_reduced)))[2]
    end
    return [inds[1], inds[2]]
end


"""
function to get a path b/w two given points on a 2-d integer grid
"""
function GetIndexPath(start::Vector{Int64}, ending::Vector{Int64} ; exclusive::Bool=true) :: Vector{Vector{Int64}}
    @assert length(start)==length(ending)==2 "The function can only return a path in 2-d"
    diff    =   ending - start

    if diff[1]!=0   ##### Slope is not zero
        xs  =   collect(start[1] : sign(ending[1] - start[1]) : ending[1] - sign(ending[1] - start[1]) * exclusive)
        ys  =   round.(Int64 , start[2] .+ (diff[2] / diff[1]) .* (xs .- start[1]))
    else
        ys  =   collect(start[2] : sign(ending[2] - start[2]) : ending[2] - sign(ending[2] - start[2]) * exclusive)
        xs  =   repeat([start[1]], length(ys))
    end

    path    =   hcat(xs, ys)
    path    =   Vector{eltype(path)}[eachrow(path)...]

    return path
end

"""
function to return a path in the brillouin zone connecting to given monentas
"""
function getBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) :: Vector{Vector{Float64}}

    start_ind   =   GetQIndex(start, bz ; nearest=nearest)
    end_ind     =   GetQIndex(ending, bz ; nearest=nearest)
    path_ind    =   GetIndexPath(start_ind, end_ind ; exclusive = exclusive)

    return getindex.(Ref(bz.ks), first.(path_ind), last.(path_ind))
end

"""
function to return a path connecting points A -> B -> C -> D -> ... when given a vector of momentum points [A, B, C, D, ...]
"""
function CombinedBZPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) :: Vector{Vector{Float64}}

    if closed
        starts  =   copy(points)
        endings =   circshift(starts, -1)
    else
        starts  =   copy(points[1:end-1])
        endings =   copy(points[2:end])
    end

    paths   =   getBZPath.(Ref(bz), starts, endings ; nearest = nearest , exclusive = true)
    paths   =   reduce(vcat, paths)

    return paths
end