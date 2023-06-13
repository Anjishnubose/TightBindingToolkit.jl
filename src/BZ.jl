module BZone
    export getRLVs , BZ , Monkhorst , fillBZ! , ReduceQ , GetQIndex , GetIndexPath , getBZPath , CombinedBZPath
 
    using LinearAlgebra
    using ..TightBindingToolkit.UCell: UnitCell


    @doc """
    ```julia
    getRLVs( uc::UnitCell ) --> Vector{Vector{Float64}}
    ```
    Returns the reciprocal lattice vectors corresponding to the given Unit Cell.

    """
    function getRLVs( uc::UnitCell ) :: Vector{Vector{Float64}}
        if length(uc.primitives) == 1
            return [2 *pi ./ uc.primitives[1]]

        elseif length(uc.primitives) == 2
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
    - `kInds           :: Vector{Matrix{Float64}}`: the [`Monkhorst`](@ref) grid corresponding to k along bs.
    - `ks   	       :: Matrix{Vector{Float64}}`: The grid of momentum points [kx, ky].
    - `HighSymPoints   :: Dict`: A dictionary containing the HIgh-Symmetry points Γ, K(2), and M(3).
    - `shift           :: Vector{Int64}` : how shifted the grid is from its centre point at the Γ point, in units of `1/gridSize`.

    Initialize this structure using 
    ```julia
    BZ(gridSize::Int64)
    ```
    """
    mutable struct BZ
        basis	        ::  Vector{Vector{Float64}}
        gridSize        ::  Int64
        kInds           ::  Vector{Matrix{Float64}}
        ks   	        ::	Matrix{Vector{Float64}}
        HighSymPoints   ::  Dict
        shift           ::  Vector{Int64}

        BZ(gridSize::Int64)    =   new{}(Vector{Float64}[], gridSize, Matrix{Float64}[], Array{Vector{Float64}}(undef, 0, 0), Dict(), Int64[])
    end


    @doc raw"""
    ```julia
    Monkhorst(ind::Int64, N::Int64) --> Float64
    Monkhorst(ind::Int64, N::Int64, shift::Int64, BC::Float64) --> Float64
    ```
    The usual Monkhorst grid is defined as follows
    `` [\frac{2i - (N+1)}{2N}, i∈[1, N]] ``
    The modified Monkhorst grid takes into account the desired boundary condition `BC`, and an integer `shift` (if required, to change the starting point), and shifts the momentum grid accordingly.

    """
    function Monkhorst(ind::Int64, N::Int64) :: Float64
        return (2 * ind - (N + 1)) / (2 * N)
    end

    function Monkhorst(ind::Int64, N::Int64, shift::Int64, BC::Float64) :: Float64
        if N%2==0.0
            return (1 / N) * (ind + shift - (N / 2) + (BC / (2 * pi)))
        else
            return (1 / N) * (ind + shift - ((N + 1) / 2) + (BC / (2 * pi)))
        end
    end


    function VecAngle(a::Vector{Float64}, b::Vector{Float64})
        return acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
    end


    function meshgrid(grid::Vector{Int64})
        inds    =   []
        for dim in 1:length(grid)
            push!(inds, 1:grid[dim])
        end

        inds    =   collect(Base.product(inds...))
        return inds
    end


    @doc """
    ```julia
    fillBZ!(bz::BZ, uc::UnitCell, offsetRange::Int64=1 ; shift::Vector{Float64}=zeros(Float64, length(uc.primitives)))
    ```
    Fills the `BZ` object with the relevant attributes, after it has been initialized as `BZ(gridsize=N)`.

    """
    function fillBZ!(bz::BZ, uc::UnitCell; shift::Vector{Int64}=zeros(Int64, length(uc.primitives)))

        bz.basis    =   getRLVs(uc)
        bz.shift    =   shift
        dims 		=	length(uc.primitives)

        @assert length(bz.basis)==2 "Sorry, the code is only written for d=2 right now!"

        inds    =   meshgrid(repeat([bz.gridSize], dims))

        for dim in 1:dims
            index   =   getindex.(inds, dim)
            push!(bz.kInds, Monkhorst.(index, Ref(bz.gridSize), Ref(shift[dim]), Ref(angle(uc.BC[dim]))))
        end

        bz.ks   =   (bz.kInds[1] .* Ref(bz.basis[1])) 
        for dim in 2:dims
            bz.ks   .+=     (bz.kInds[dim] .* Ref(bz.basis[dim]))
        end 

        bz.HighSymPoints["Gamma"]   =   zeros(Float64, dims)

        if dims==2

            bz.HighSymPoints["M1"]      =  @. 0.5 * bz.basis[1] + 0.0 * bz.basis[2]
            bz.HighSymPoints["M2"]      =  @. 0.0 * bz.basis[1] + 0.5 * bz.basis[2]
            bz.HighSymPoints["M3"]      =  @. 0.5 * bz.basis[1] + 0.5 * bz.basis[2]
            bz.HighSymPoints["M4"]      =  @. -0.5 * bz.basis[1] + 0.0 * bz.basis[2] ###Added this to test

            if isapprox(VecAngle(bz.basis[1], bz.basis[2]), 2*pi/3, atol=1e-4, rtol=1e-4)
                bz.HighSymPoints["K1"]      =   @. (2/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2"]      =   @. (1/3) * bz.basis[1] + (2/3) * bz.basis[2]
            elseif isapprox(VecAngle(bz.basis[1], bz.basis[2]), pi/3, atol=1e-4, rtol=1e-4)
                bz.HighSymPoints["K1"]      =   @. (1/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2"]      =   @. (2/3) * bz.basis[1] + (2/3) * bz.basis[2]
            end
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
        Q_reduced   =   inv(U) * (Q)  ##### Q in the basis of b1 and b2
        Q_reduced   =   Q_reduced .- round.(Q_reduced .- (bz.shift ./ bz.gridSize))  ##### Q in the basis of b1 and b2, and shifted back to the first BZ
        Q_reduced   =   sum((Q_reduced) .* bz.basis)   ##### Q back in the x ,y basis, but shifted to the first BZ
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
        return collect(Tuple(inds))
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

end