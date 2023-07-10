module Chern
    export FindLinks , FieldStrength , ChernNumber, CheckValidity

    using ..TightBindingToolkit.Hams:Hamiltonian, IsBandGapped

    using LinearAlgebra


    @doc """
    ```julia
    FindLinks(Ham::Hamiltonian, subset::Vector{Int64}) --> Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
    ```
    Function to get the linking matrices on each neighbouring point in the `BZ`.
    On a bond connecting k_i and k_j, the linking matrix U is defined such that U[m, n] = <v^m[k_i]|v^n[k_j]> where states[k_j[1], k_j[2], :, m] 
    v^m[k_j], the mth eigenstate at momentum k_j.

    """
    function FindLinks(Ham::Hamiltonian, subset::Vector{Int64})::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
        shifted_1   =   circshift(Ham.states, [-1, 0])
        shifted_2   =   circshift(Ham.states, [0, -1])

        Link_1        =   det.(selectdim.(adjoint.(Ham.states), 1, Ref(subset)) .* selectdim.(shifted_1, 2, Ref(subset)))
        Link_2        =   det.(selectdim.(adjoint.(Ham.states), 1, Ref(subset)) .* selectdim.(shifted_2, 2, Ref(subset)))
        ##### selectdim(x, 1, v) = x[v, :] and similarly selectdim(x, 2, v) = x[:, v]
        ##### selectdim.(M, 1, Ref(subset)) ---> [[M[subset, :] for all k points]] 
        return (Link_1, Link_2)
    end


    @doc """
    ```julia
    FieldStrength(Links::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}) --> Matrix{ComplexF64}
    ```
    Function to calculate the product of the links over each plaquette on the BZ grid. This is the generalized Bery curvature for multiple degenerate bands.

    """
    function FieldStrength(Links::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}})::Matrix{ComplexF64}

        Fields   =   Links[1] .* circshift(Links[2], [-1, 0]) .* circshift(conj.(Links[1]), [0, -1]) .* conj.(Links[2])
        """
        (k+a2)*<--(Links[1])†[k+a2]---*
            |                       ^
            |                       |
        (Links[2])†[k]         Links[2][k+a1]
            |                       |
            |                       |
        (k)*-->Links[1][k]-------->*(k+a1)
        """
        return Fields
    end


    function CheckValidity(Ham::Hamiltonian, subset::Vector{Int64})
        bandGapped  =   IsBandGapped(Ham)
        bandGapped[begin, begin]    =   true
        bandGapped[end, end]        =   true

        occupied_bands  =   collect(extrema(subset))    #### Just need to check that the lowest and highest bands given are gapped or not
        @assert subset == collect(UnitRange(occupied_bands...)) "Cannot skip bands in between!"

        valence_bands   =   clamp.(occupied_bands + [-1, 1], Ref(1), Ref(length(Ham.bands[begin])))   ##### one below and one above the lowest and highest band given
        check           =   prod(getindex.(Ref(bandGapped), occupied_bands, valence_bands))

        @assert check "Given subset of bands have band touchings / degeneracies with other bands. Chern number is not well defined! "

    end


    @doc """
    ```julia
    ChernNumber(Ham::Hamiltonian, subset::Vector{Int64}) --> Float64
    ```
    Function to get Chern numbers given a `Hamiltonian` and a `subset` of bands

    """
    function ChernNumber(Ham::Hamiltonian, subset::Vector{Int64} ; check_validity::Bool = false)::Float64

        if check_validity
            CheckValidity(Ham, subset)
        end

        Links   =   FindLinks(Ham, subset)
        Field   =   FieldStrength(Links)
        chern   =   (1/(2*pi)) * sum(angle.(Field))
        return chern
    end

end