include("UnitCell.jl")
include("BZ.jl")
include("Hamiltonian.jl")

using LinearAlgebra

"""
Function to get the linking matrices on each neighbouring point in the brillouin zone grid.
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

"""
Function to calculate the product of the links over each plaquette on the BZ grid
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

"""
Function to get Chern numbers given a Hamiltonian and a subset of bands
"""
function ChernNumber(Ham::Hamiltonian, subset::Vector{Int64})::Float64
    Links   =   FindLinks(Ham, subset)
    Field   =   FieldStrength(Links)
    chern   =   (1/(2*pi)) * sum(angle.(Field))
    return chern
end