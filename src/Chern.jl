module Chern
    export FindLinks , FieldStrength , ChernNumber, CheckValidity, PartialChernNumber, FilledChernNumber, OccupiedChernNumber, KuboChern

    using ..TightBindingToolkit.Hams:Hamiltonian, IsBandGapped, GetVelocity!
    using ..TightBindingToolkit.Useful: DistFunction
    using ..TightBindingToolkit.BZone:BZ

    using LinearAlgebra


@doc """
```julia
FindLinks(Ham::Hamiltonian, subset::Vector{Int64}) --> Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
```
Function to get the linking matrices on each neighbouring point in the `BZ`.
On a bond connecting k_i and k_j, the linking matrix U is defined such that U[m, n] = <v^m[k_i]|v^n[k_j]> where states[k_j[1], k_j[2]][:, m] = v^m[k_j], the mth eigenstate at momentum k_j.

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


@doc """
```julia
PartialChernNumber(Ham::Hamiltonian, band::Int64, mu::Float64) --> Float64
Function to get the Chern number of a partially filled band given a `Hamiltonian`, a `band` index, and a chemical potential `mu`.

"""
    function PartialChernNumber(Ham::Hamiltonian, band::Int64, mu::Float64, T::Float64)::Float64

        @assert band in UnitRange(1, length(Ham.bands[begin])) "Band index out of range!"

        Links   =   FindLinks(Ham, [band])
        Field   =   FieldStrength(Links)

        energies=   getindex.(Ham.bands, Ref(band))
        filled  =   DistFunction(energies; mu = mu, T = T, stat = -1)
        chern   =   (1/(2*pi)) * sum((angle.(Field)) .* filled )

        return chern
    end


@doc """
```julia
FilledChernNumber(Ham::Hamiltonian, mu::Float64) --> Float64
Function to get the Chern number of bands filled upto a given chemical potential `mu`.

"""
    function FilledChernNumber(Ham::Hamiltonian, mu::Float64, T::Float64)::Float64

        # filled_bands    =   searchsortedfirst.(Ham.bands, Ref(mu)) .- 1

        # if findmax(filled_bands)[1] == 0
        #     @warn "Chemical potential is below the lowest band. Chern number is not well defined!"
        #     return 0.0

        # else
        return sum(PartialChernNumber.(Ref(Ham), collect(1:length(Ham.bands[begin])), Ref(mu), Ref(T)))
        # end

    end

    function OccupiedChernNumber(Ham::Hamiltonian, mu::Float64, T::Float64)
        chern = 0.0

        for band in 1:length(Ham.bands[begin])
            link = FindLinks(Ham, [band])
            field = FieldStrength(link)
            occupation = DistFunction(getindex.(Ham.bands, band); mu = mu, T = T)
            chern += (1/(2*pi)) * sum((angle.(field)) .* occupation)
        end

        return chern
    end

    function KuboChern(Ham::Hamiltonian, bz::BZ, mu::Float64)

        Vx = conj.(permutedims.(Ham.states)) .* Ham.velocity[1] .* Ham.states
        Vy = conj.(permutedims.(Ham.states)) .* Ham.velocity[2] .* Ham.states

        chern = 0.0 + im*0.0
        for k in eachindex(Ham.bands)
            Es = Ham.bands[k]
            vx = Vx[k]
            vy = Vy[k]

            ind = searchsortedfirst(Es, mu)
            if ind == 1 || ind == length(Es)
                continue
            else
                for i in 1:ind-1
                    for j in ind:length(Es)
                        chern += (vx[i, j] * vy[j, i] - vx[j, i] * vy[i, j]) / ((Es[j] - Es[i])^2)
                    end
                end
            end

        end

        b1 = [bz.basis[1];0.0]
        b2 = [bz.basis[2];0.0]
        bzUnitArea = cross(b1, b2)[3]/(4*pi^2)

        return imag(chern)*bzUnitArea*2*pi/length(Ham.bands)

    end


end
