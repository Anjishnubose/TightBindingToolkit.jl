module TBModel

    export Model , FindFilling , GetMu! , GetFilling! , GetCount , GetGk! , GetGr!, SolveModel!

    using ..TightBindingToolkit.Useful: Meshgrid, DistFunction, DeriDistFunction, FFTArrayofMatrix, BinarySearch
    using ..TightBindingToolkit.UCell:UnitCell
    using ..TightBindingToolkit.BZone:BZ, MomentumPhaseFFT
    using ..TightBindingToolkit.Hams:Hamiltonian

    using LinearAlgebra, Logging


    @doc """
    `Model` is a data type representing a general Tight Binding system.

    # Attributes
    - `uc      ::  UnitCell`: the Unit cell of the lattice.
    - `bz      ::  BZ`: The discretized Brillouin Zone.
    - `Ham     ::  Hamiltonian`: the Hamiltonian at all momentum-points.
    - `T       ::  Float64`: the temperature of the system.
    - `filling ::  Float64`: The filling of the system.
    - `mu      ::  Float64`: The chemical potential of the system.
    - `stat    ::  Int64` : ±1 for bosons and fermions.
    - `gap     ::  Float64` : the energy gap of excitations at the given filling.
    - `Gk      ::  Array{Matrix{ComplexF64}}` : An array (corresponding to the grid of k-points in `BZ`) of Greens functions.

    Initialize this structure using 
    ```julia
    Model(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)
    ```
    You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
    """
    mutable struct Model
        uc      ::  UnitCell
        bz      ::  BZ
        Ham     ::  Hamiltonian
        T       ::  Float64         ##### Temperature
        filling ::  Float64         ##### Filling fraction
        mu      ::  Float64         ##### Chemical potential
        stat    ::  Int64           ##### +1 for bosons, -1 for fermions
        gap     ::  Float64
        Gk      ::  Array{Matrix{ComplexF64}}
        Gr      ::  Array{Matrix{ComplexF64}}
        
        Model(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1) = new{}(uc, bz, Ham, T, filling, mu, stat, 0.0, Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc.primitives))...))
        ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
    end


    


    @doc """
    ```julia
    FindFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1) --> Float64
    ```
    Find filling at given `T`=temperature and `mu`=chemical potential, for a given `bands`.

    """
    function FindFilling( mu::Float64, bands::Vector{Float64}, T::Float64, stat::Int64=-1) :: Float64
        return sum(DistFunction(bands; T=T, mu=mu, stat=stat)) / length(bands)
    end

    

    @doc """
    ```julia
    GetMu!(M::Model; tol::Float64=0.001)
    ```
    Function to get chemical potential for a given `Model`, within a tolerance.

    """
    function GetMu!(M::Model ;  mu_tol::Float64 = 0.001, filling_tol::Float64 = 1e-6)
    
        if M.T≈0.0
            energies  =   sort(reduce(vcat, M.Ham.bands))
            M.mu      =   (energies[floor(Int64, length(energies)*M.filling)] + energies[floor(Int64, length(energies)*M.filling)+1])/2
        else
            M.mu      =   BinarySearch(M.filling, M.Ham.bandwidth, FindFilling, (reduce(vcat, M.Ham.bands), M.T, M.stat) ; x_tol=mu_tol, target_tol = filling_tol)
            @info "Found chemical potential μ = $(M.mu) for given filling = $(M.filling)."
        end
    end


    @doc """
    ```julia
    GetFilling!(M::Model)
    ```
    Find filling for a given `Model`.

    """
    function GetFilling!(M::Model)
        if M.T≈0.0
            energies    =   sort(reduce(vcat, M.Ham.bands))
            M.filling   =   searchsortedlast(energies, M.mu) / length(energies)
        else
            M.filling   =   FindFilling( M.mu, reduce(vcat, M.Ham.bands), M.T, M.stat)
        end
    end


    @doc """
    ```julia
    GetCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) --> Matrix{Float64}
    ```
    Function to return a diagonal matrix whose entries are M[i, i] = θ(-(E^i(k)-μ)) ----> 1 if the energy is below the chemical potential, otherwise 0.

    """
    function GetCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) :: Matrix{Float64}
        return diagm(DistFunction(Es; T=T, mu=mu, stat=stat))
    end


    @doc """
    ```julia
    GetGk!(M::Model)
    ```
    Finding the equal-time Greens functions in momentum space of a `Model`.

    """
    function GetGk!(M::Model)
        quasiCount 	=	GetCount.(M.Ham.bands, Ref(M.mu), Ref(M.T), Ref(M.stat))   ##### Matrix with 1s and 0s along the diagonal. The count of the quasiparticles at each k point determined by the bands and the chemical potential
        M.Gk        =   transpose.(M.Ham.states .* quasiCount .* adjoint.(M.Ham.states))
    end


    @doc """
    ```julia
    GetGr!(M::Model)
    ```
    Finding the equal-time Greens functions in real space of a `Model`.

    """
    function GetGr!(M::Model)
        
        Gr           =  FFTArrayofMatrix(M.Gk)  
        phaseShift   =  MomentumPhaseFFT(M.bz, M.uc)
        M.Gr         =  Gr .* phaseShift   
    end


    @doc """
    ```julia
    SolveModel!(M::Model)
    ```
    one-step function to find all the attributes in Model after it has been initialized.

    """
    function SolveModel!(M::Model ; get_correlations::Bool = true, verbose::Bool = true, mu_tol::Float64 = 1e-3, filling_tol::Float64 = 1e-6)
        @assert M.Ham.is_BdG==false "BdG Hamiltonians should be solved using a BdGModel"

        if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
            GetFilling!(M)
        else
            GetMu!(M ; mu_tol = mu_tol, filling_tol = filling_tol)
        end

        energies  =   sort(reduce(vcat, M.Ham.bands))
        M.gap     =   energies[floor(Int64, length(energies)*M.filling) + 1] - energies[floor(Int64, length(energies)*M.filling)]
    
        if get_correlations
            GetGk!(M)
            GetGr!(M)
        end
        
        if verbose
            @info "System Filled!"
        end
    end

end