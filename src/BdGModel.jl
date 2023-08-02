module BdG
    export BdGModel, FindFilling, GetMu!, GetFilling!, GetGk!, GetGr!, SolveModel!, GetGap!, FreeEnergy

    using ..TightBindingToolkit.Useful: DistFunction, DeriDistFunction, BinarySearch, FFTArrayofMatrix
    using ..TightBindingToolkit.UCell: UnitCell
    using ..TightBindingToolkit.DesignUCell: ModifyIsotropicFields!
    using ..TightBindingToolkit.BZone:BZ, MomentumPhaseFFT
    using ..TightBindingToolkit.Hams:Hamiltonian, ModifyHamiltonianField!
    
    using LinearAlgebra, Tullio, TensorCast, Logging

    import ..TightBindingToolkit.TBModel:FindFilling, GetMu!, GetFilling!, GetGk!, GetGr!, SolveModel!, GetGap!, FreeEnergy

@doc """
    `BdGModel` is a data type representing a general Tight Binding system with pairing.

    # Attributes
    - `uc_hop  ::  UnitCell`: the Unit cell of the lattice with the hoppings.
    - `uc_pair ::  UnitCell`: the Unit cell of the lattice with the pairings.
    - `bz      ::  BZ`: The discretized Brillouin Zone.
    - `Ham     ::  Hamiltonian`: the Hamiltonian at all momentum-points.
    - `T       ::  Float64`: the temperature of the system.
    - `filling ::  Float64`: The filling of the system.
    - `mu      ::  Float64`: The chemical potential of the system.
    - `stat    ::  Int64` : ±1 for bosons and fermions.
    - `gap     ::  Float64` : the energy gap of excitations at the given filling.
    - `Gk      ::  Array{Matrix{ComplexF64}}` : An Array (corresponding to the grid of k-points in `BZ`) of Greens functions.
    - `Fk      ::  Array{Matrix{ComplexF64}}` : An Array (corresponding to the grid of k-points in `BZ`) of anomalous Greens functions.

    Initialize this structure using 
    ```julia
    BdGModel(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)
    ```
    You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
"""
    mutable struct BdGModel
        uc_hop  ::  UnitCell
        uc_pair ::  UnitCell
        bz      ::  BZ
        Ham     ::  Hamiltonian
        """
        Thermodynamic properties
        """
        T       ::  Float64         ##### Temperature
        filling ::  Float64         ##### Filling fraction
        mu      ::  Float64         ##### Chemical potential
        gap     ::  Float64
        stat    ::  Int64           ##### +1 for Boson, -1 for Fermions
        """
        Correlations
        """
        Gk      ::  Array{Matrix{ComplexF64}}
        Fk      ::  Array{Matrix{ComplexF64}}
        Gr      ::  Array{Matrix{ComplexF64}}
        Fr      ::  Array{Matrix{ComplexF64}}
        
        BdGModel(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1) = new{}(uc_hop, uc_pair, bz, Ham, T, filling, mu, -999.0, stat ,Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...),Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc_hop.primitives))...))
        ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
    end


@doc """
```julia
FindFilling(M::BdGModel) --> Float64
FindFilling(mu::Float64, M::BdGModel)
```
Find filling at given temperature and chemical potential that takes BdGModel object as argument instead of Hamiltonian, since for a BdG case, the filling depends on the wavefunctions also.
Because of this, if you want to calculate the filling at a different chemical potential, you have to modify the Hamiltonian, the UnitCells, rediagonalize, and recalculate eveyrthing.

"""
    function FindFilling(M::BdGModel) :: Float64
        @assert M.Ham.is_BdG==true "Use other format for pure hopping Hamiltonian"

        N       =   length(M.Ham.bands[begin])
        Eks     =   reshape(getindex.(M.Ham.bands, Ref(1:N÷2)), prod(M.bz.gridSize))
        U11s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(1:N÷2)), prod(M.bz.gridSize))
        U21s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(1:N÷2)), prod(M.bz.gridSize))

        nFs     =   DistFunction.(Eks; T=M.T, mu=0.0, stat=M.stat)

        @tullio filling := ((conj(U11s[k1][i, j]) * U11s[k1][i, j] * nFs[k1][j]) 
                          + (conj(U21s[k1][i, j]) * U21s[k1][i, j] * (1 - nFs[k1][j])))

        return real(filling) / (prod(M.bz.gridSize) * M.uc_hop.localDim * length(M.uc_hop.basis))
    end

    function FindFilling(mu::Float64, M::BdGModel) :: Float64
        ModifyHamiltonianField!(M.Ham, M.uc_hop, repeat([mu], outer=length(M.uc_hop.basis)))
        return FindFilling(M)
    end


"""
```julia
GetFilling!(M::BdGModel)
```
Function to get filling given a chemical potential!

"""
    function GetFilling!(M::BdGModel)
        M.filling   =   FindFilling(M.mu, M)
    end


"""
```julia
GetMu!(M::BdGModel ;  tol::Float64=0.001)
```
Function to get the correct chemical potential given a filling.
"""
    function GetMu!(M::BdGModel ;  mu_tol::Float64 = 0.001, filling_tol::Float64 = 1e-6)
        M.mu    =   BinarySearch(M.filling, M.Ham.bandwidth, FindFilling, (M,) ; x_tol = mu_tol, target_tol = filling_tol)
        @info "Found chemical potential μ = $(M.mu) for given filling = $(M.filling)."
    
    end


"""
```julia
GetGk!(M::BdGModel)
```
Finding the Greens functions, and anomalous greens functions in momentum space at some chemical potential.
"""
    function GetGk!(M::BdGModel)

        N       =   length(M.Ham.bands[begin])
        Eks     =   reshape(getindex.(M.Ham.bands, Ref(1:N÷2)) , prod(M.bz.gridSize))   ##### Only the negative energies from the bdG spectrum

        U11s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(1:N÷2)) , prod(M.bz.gridSize))          ##### The 4 different quadrants in the Unitary relating the BdG quasiparticles and the nambu basis
        U21s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(1:N÷2)) , prod(M.bz.gridSize)) 
        U12s    =   reshape(getindex.(M.Ham.states, Ref(1:N÷2), Ref(N÷2 + 1: N)) , prod(M.bz.gridSize)) 
        U22s    =   reshape(getindex.(M.Ham.states, Ref(N÷2 + 1: N), Ref(N÷2 + 1: N)) , prod(M.bz.gridSize))

        nFs     =   DistFunction.(Eks; T=M.T, mu=0.0, stat=M.stat)

        @reduce Gk[k1][i, j] |= sum(l) ((conj(U11s[k1][i, l]) * U11s[k1][j, l] * nFs[k1][l]) 
                                      + (conj(U12s[k1][i, l]) * U12s[k1][j, l] * (1 - nFs[k1][l])))

        @reduce Fk[k1][i, j] |= sum(l) ((conj(U11s[k1][i, l]) * U21s[k1][j, l] * nFs[k1][l]) 
                                      + (conj(U12s[k1][i, l]) * U22s[k1][j, l] * (1 - nFs[k1][l])))
        M.Gk    =   reshape(Gk, M.bz.gridSize...)
        M.Fk    =   reshape(Fk, M.bz.gridSize...)

    end


@doc """
```julia
GetGr!(M::BdGModel)
```
Finding the equal-time Greens functions and anomalous Greens function in real space of a `BdGModel`.

"""
    function GetGr!(M::BdGModel)
        
        phaseShift   =  MomentumPhaseFFT(M.bz, M.uc_hop)

        Gr           =  FFTArrayofMatrix(M.Gk)  
        M.Gr         =  Gr .* phaseShift   

        Fr           =  FFTArrayofMatrix(M.Fk)  
        M.Fr         =  Fr .* phaseShift   
    end


@doc """
```julia
SolveModel!(M::BdGModel)
```
one-step function to find all the attributes in  BdGModel after it has been initialized.
"""
    function SolveModel!(M::BdGModel ; get_correlations::Bool = true, get_gap::Bool = false, verbose::Bool = true, mu_tol::Float64 = 1e-3, filling_tol::Float64 = 1e-6)
        @assert M.Ham.is_BdG==true "Use other format for pure hopping Hamiltonian"

        if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
            GetFilling!(M)
        else
            GetMu!(M ; mu_tol = mu_tol, filling_tol = filling_tol)
        end

        if get_gap

            energies    =   sort(reduce(vcat, M.Ham.bands))
            M.gap       =   energies[min(floor(Int64, length(energies)*0.5) + 1, length(energies))] - energies[floor(Int64, length(energies)*0.5)]
        end

        if get_correlations
            GetGk!(M)
            GetGr!(M)
        end
        
        if verbose
            @info "System Filled!"
        end
    end


@doc """
```julia
GetGap!(M::BdGModel)
```
Calculate the BdG gap of the system.
"""
    function GetGap!(M::BdGModel)
        energies    =   sort(reduce(vcat, M.Ham.bands))
        M.gap       =   energies[min(floor(Int64, length(energies)*0.5) + 1, length(energies))] - energies[floor(Int64, length(energies)*0.5)]

    end


@doc """
```julia
FreeEnergy(M::BdGModel; F0::Float64 = 0.0) --> Float64
```
Calculate the free energy of the given `BdGModel`.

"""
    function FreeEnergy(M::BdGModel ; F0::Float64 = 0.0) :: Float64

        Es      =   reduce(vcat, M.Ham.bands)
        F       =   log.(1 .+ exp.(-(Es .- M.mu) / (M.T)))
        F       =   -M.T * (sum(F) / length(F))

        return F - F0
    end


end