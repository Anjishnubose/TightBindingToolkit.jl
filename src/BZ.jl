using LinearAlgebra
include("UnitCell.jl")


@doc """
```julia
getRLVs( uc::UnitCell ) --> Vector{Vector{Float64}}
```
Returns the reciprocal lattice vectors corresponding to the given Unit Cell.

"""
function getRLVs( uc::UnitCell ) :: Vector{Vector{Float64}}
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

@doc """
`BZ` is a data type representing a discretized Brillouin Zone in momentum space.

# Attributes
 - `basis           :: Vector{ Vector{ Float64 } }`: reciprocal lattice vectors of the Brillouin Zone.
 - `gridSize        :: Int64`: The number of points along each dimension of the grid.
 - `k1s             :: Matrix{Float64}`: the [`Monkhorst`](@ref) grid corresponding to k along b1.
 - `k2s             :: Matrix{Float64}`: the [`Monkhorst`](@ref) grid corresponding to k along b2.
 - `ks   	        :: Matrix{Vector{Float64}}`: The grid of momentum points [kx, ky].
 - `HighSymPoints   :: Dict`: A dictionary containing the HIgh-Symmetry points Γ, K(2), and M(3).
 - `shift           ::  Vector{Float64}` : how shifted the grid is from its centre point at the Γ point.

Initialize this structure using 
```julia
BZ(gridSize::Int64)
```
"""
mutable struct BZ
    basis	        ::  Vector{Vector{Float64}}
    gridSize        ::  Int64
    k1s             ::  Matrix{Float64}
    k2s             ::  Matrix{Float64}
	ks   	        ::	Matrix{Vector{Float64}}
    HighSymPoints   ::  Dict
    shift           ::  Vector{Float64}

    BZ(gridSize::Int64)    =   new{}([], gridSize, Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0),  Array{Vector{Float64}}(undef, 0, 0), Dict(), [])
end


@doc raw"""
```julia
Monkhorst(ind::Int64, N::Int64) --> Float64
```
The Monkhorst grid is defined as follows
`` [\frac{2i - (N+1)}{2N}, i∈[1, N]] ``

"""
function Monkhorst(ind::Int64, N::Int64) :: Float64
    return (2 * ind - (N + 1)) / (2 * N)
end


function VecAngle(a::Vector{Float64}, b::Vector{Float64})
    return acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end


@doc """
```julia
fillBZ!(bz::BZ, uc::UnitCell, offsetRange::Int64=1 ; shift::Vector{Float64}=zeros(Float64, length(uc.primitives)))
```
Fills the `BZ` object with the relevant attributes, after it has been initialized as `BZ(gridsize=N)`.

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


@doc """
```julia
ReduceQ(Q::Vector{Float64}, bz::BZ) --> Vector{Float64}
```
Reduces a given momentum back to the range covered by the discretized Brillouin Zone.

"""
function ReduceQ(Q::Vector{Float64}, bz::BZ) :: Vector{Float64}
    U           =   reduce(hcat, bz.basis) ##### basis transformation from x, y -> b1, b2
    Q_reduced   =   inv(U) * (Q .- bz.shift)   ##### Q in the basis of b1 and b2
    Q_reduced   =   Q_reduced .- round.(Q_reduced)  ##### Q in the basis of b1 and b2, and shifted back to the first BZ
    Q_reduced   =   sum(Q_reduced .* bz.basis) .+ bz.shift ##### Q back in the x ,y basis, but shifted to the first BZ
    return Q_reduced
end


@doc """
```julia
GetQIndex(Q::Vector{Float64}, bz::BZ ; nearest::Bool = false) --> Vector{Int64}
```
Returns the index in the discretized `BZ` of the momentum point corresponding to the fiven momentum `Q`. 
If the input `nearest` is set to `true`, will return the index of the momentum point on the grid closest to `Q`, if `Q` does not exist on the grid. 

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


@doc """
```julia
GetIndexPath(start::Vector{Int64}, ending::Vector{Int64} ; exclusive::Bool=true) --> Vector{Vector{Int64}}
```
Returns a path in index-space of the discretized `BZ` which joins the two sets of indices `start` and `ending`. 
If the input `exclusive` is set to `true`, the returned path will NOT contain the `ending` point itself.

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


@doc """
```julia
getBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) --> Vector{Vector{Float64}}
```
Returns the actual path in momentum-space of the discretized `BZ` which joins the two momentums `start` and `ending`. 
The optional input `nearest` is the same as in [`GetQIndex`](@ref), and `exclusive` is the same as in [`GetIndexPath`](@ref).

"""
function getBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) :: Vector{Vector{Float64}}

    start_ind   =   GetQIndex(start, bz ; nearest=nearest)
    end_ind     =   GetQIndex(ending, bz ; nearest=nearest)
    path_ind    =   GetIndexPath(start_ind, end_ind ; exclusive = exclusive)

    return getindex.(Ref(bz.ks), first.(path_ind), last.(path_ind))
end


@doc """
```julia
CombinedBZPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) --> Vector{Vector{Float64}}
```
Returns a path in momentum-space of the discretized `BZ` which joins the given momentum points present in `points` as point[1] --> point[2] --> ... --> point[end] --> point[1].
The optional input `nearest` is the same as in [`GetQIndex`](@ref), and `closed` determines if the path is a clsoed loop or not.

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