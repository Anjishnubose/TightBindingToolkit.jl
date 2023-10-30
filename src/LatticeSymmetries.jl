module LatticeSymmetries
    export Translations, Degeneracies, FindQuantumNumbers

    using LinearAlgebra, StatsBase

    using ..TightBindingToolkit.UCell: UnitCell
    using ..TightBindingToolkit.LatticeStruct: Lattice
    using ..TightBindingToolkit.LatHam: LatticeHamiltonian
    using ..TightBindingToolkit.LatModel: LatticeModel


    function Translations(lattice::Lattice{T} ; primitive::Int64, with_local::Bool = false) :: Matrix{ComplexF64} where {T}

        Ta     =   zeros(ComplexF64, lattice.length, lattice.length) 
        id     =   Matrix(I, length(lattice.uc.primitives), length(lattice.uc.primitives))

        for j in 1:lattice.length

            translated  =   mod.(lattice.sites[j][end] + id[primitive, :], lattice.size)
            Ta[lattice.sites((lattice.sites[j][begin], translated)), j]     =   1.0
        end

        if with_local
            Ta  =   kron(Ta, Matrix(I, lattice.uc.localDim, lattice.uc.localDim))
        end

        return Ta
    end


    function Degeneracies(H::LatticeHamiltonian ; tol::Int64 = 6) :: Vector{UnitRange{Int64}}

        counts      =   countmap(round.(H.bands, digits = tol))

        if haskey(counts, 0.0) && haskey(counts, -0.0)
            counts[0.0] =   counts[0.0] + counts[-0.0]
            delete!(counts, -0.0)
        end

        energies    =   sort(collect(keys(counts)))
        counts      =   getindex.(Ref(counts), energies)

        endings     =   cumsum(counts)
        beginnings  =   [1; endings[1:end-1] .+ 1]
        
        degeneracies=   UnitRange.(beginnings, endings)

        return degeneracies
    end


    function FindQuantumNumbers(M::LatticeModel, SymMatrix::Matrix{ComplexF64} ; tol::Int64 = 6, till_band::Int64) :: Vector{Vector{ComplexF64}}

        degens  =   Degeneracies(M.Ham ; tol = tol)

        nBlocks =   searchsortedfirst(cumsum(length.(collect.(degens))), till_band)
        wavefunctions   =   getindex.(Ref(M.Ham.states), Ref(1:size(M.Ham.states, 1)), degens[1:nBlocks])

        SymMatrixBlocks =   adjoint.(wavefunctions) .* Ref(SymMatrix) .* wavefunctions
        quantumNumbers  =   eigvals.(SymMatrixBlocks)

        return quantumNumbers
    end


end