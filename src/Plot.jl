module PlotTB
    export plot_UnitCell! ,plot_band_contour!, plot_band_structure!, plot_FS!

    using LinearAlgebra, LaTeXStrings, Plots

    using ..TightBindingToolkit.UCell:UnitCell, getAllOffsets
    using ..TightBindingToolkit.BZone:BZ, GetQIndex, getBZPath, CombinedBZPath, ReduceQ
    using ..TightBindingToolkit.Hams:Hamiltonian
    using ..TightBindingToolkit.TBModel:Model
    using ..TightBindingToolkit.BdG:BdGModel


    @doc """
    ```julia
    plot_UnitCell!(uc::UnitCell ; range::Int64 = 1, cmp::Symbol = :matter, plot_conjugate::Bool=true) --> Plots.plot()
    ```
    Function to plot the `UnitCell`. 
    `range` determines the range of UnitCells plotted in real-space. default is ±1.
    `cmp` determines the colorScheme of chosen to differentiate b.w different hopping types.
    `plot_conjugate` switches on whether the conjugate of the bonds in `UnitCell` are also plotted or not.

    """
    function plot_UnitCell!(uc::UnitCell ; range::Int64 = 1, cmp::Symbol = :matter, plot_conjugate::Bool=true)
        
        dim         =   length(uc.primitives)
        @assert dim == 2 "Unit Cell plotting only works for 2d right now!"
        offsets     =   getAllOffsets(range, dim)

        p       =   plot(aspect_ratio=:equal, grid=false)
        ##### Plotting sites
        for offset in offsets
            shift   =   sum(offset .* uc.primitives)
            sites   =   uc.basis .+ Ref(shift)
            ##### Lattice sites
            scatter!(Tuple.(sites), label = "", markercolor=:orange, markersize=16)

            if offset == zeros(Int64, dim)
                for (i, site) in enumerate(sites)
                    annotate!(site..., Plots.text(L"\mathbf{%$i}", :bottom, :left, 16, :black))
                end
                ##### Lattice primitive vectors
                for primitive in uc.primitives
                    base    =   uc.basis[begin]
                    target  =   uc.basis[begin] .+ primitive
                    plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=true , color=:green , linewidth=6 , label = "", linealpha=0.75)
                end

            end

        end

        ##### Different labels, colors, and thickness of each bond types.
        labels      =   unique(getproperty.(uc.bonds, :label))
        cmp         =   cgrad(cmp, max(2 , length(labels)), categorical=true, rev=true)
        thickness   =   collect(4:-2/(max(2 , length(labels)) -1):2)   
        counts      =   zeros(Int64, length(labels))
        
        ##### Plotting bonds
        for bond in uc.bonds
            base    =   uc.basis[bond.base]
            target  =   uc.basis[bond.target] .+ sum(bond.offset .* uc.primitives)
            target_conj     =   uc.basis[bond.target] .- sum(bond.offset .* uc.primitives)

            index   =   findfirst(==(bond.label), labels)
            counts[index]   +=  1

            if counts[index] == 1
                plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=true , color=cmp[index] , linewidth=thickness[index] , label = L"%$(bond.label)", linealpha=0.75)
            else
                plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=true , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=0.75)
            end
            ##### Plotting the conjugate of each bond
            if plot_conjugate
                plot!([target_conj[begin] , base[begin]] , [target_conj[end], base[end]] , arrow=true , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=0.75)
            end
        end

        xlabel!(L"x")
        ylabel!(L"y")
        title!("Unit Cell")

        return p

    end


    @doc """
    ```julia
    plot_band_contour!(Ham::Hamiltonian , bz::BZ , band_index::Int64) --> Plots.plot()
    ```
    Function to draw equal energy contours of the bands in `Hamiltonian`, specifically for the band with the given `band_index`.

    """
    function plot_band_contour!(Ham::Hamiltonian , bz::BZ , band_index::Int64 ; cmp::Symbol = :turbo)
        pyplot()
        @assert band_index <= length(Ham.bands[begin]) "Given band does not exist in the given Hamiltonian!"
        @assert length(size(Ham.bands)) == 2 "Contour plots only work for 2d Hamiltonians"

        xgrid = getindex.(bz.ks , 1)
        ygrid = getindex.(bz.ks , 2)

        plt   = contourf(xgrid , ygrid , getindex.(Ham.bands , band_index), levels=20, color=cmp, size=(600,600))
        xlabel!(L"k_x")
        ylabel!(L"k_y")
        title!(L"Band :  %$band_index")
        return plt
    end 


    @doc """
    ```julia
    plot_band_structure!(M<:Union{Model, BdGModel}, path::Vector{Vector{Float64}},  band_index::Vector{Int64} = collect(1:length(M.Ham.bands[begin])) ; labels::Vector{} = repeat([""], length(path)), closed::Bool=true, nearest::Bool=true) --> Plots.plot()
    ```
    Function to plot band structures of the `Model` along a `path` in the BZ determined by the given critical points.
    Can take in multiple bands into account ∈ `band_index`.
    `labels` are the Plot labels of the critical points.

    """
    function plot_band_structure!(M::T, path::Vector{Vector{Float64}},  band_index::Vector{Int64} = collect(1:length(M.Ham.bands[begin])) ; labels::Vector{} = repeat([""], length(path)), closed::Bool=true, nearest::Bool=true) where {T<:Union{Model, BdGModel}}
        
        bzpath     = CombinedBZPath(M.bz, path ; nearest=nearest, closed = closed)
        path_index = GetQIndex.(bzpath, Ref(M.bz))

        bands_from_index = getindex.(Ref(M.Ham.bands), CartesianIndex.(Tuple.(path_index)))

        label_indices    = getindex.(findmin.([norm.(Ref(ReduceQ(x,M.bz)).-bzpath) for x in path]),2)

        plt = plot(grid=false)
        for j in band_index
            plot!(getindex.(bands_from_index , j), labels=L"Band : %$j", lw = 2.0)
        end
        if typeof(M)==BdGModel
            hline!([0.0], linestyle = :dash, label="") ##### Effective chemical potential in BdG
        else
            hline!([M.mu], linestyle = :dash, label="") 
        end

        xticks!(label_indices, labels)
        xlabel!("Path", guidefontsize = 9)
        ylabel!("Energy", guidefontsize = 9)
        title!("Band Structure along path", titlefontsize = 12)

        return plt
    end


    @doc """
    ```julia
    plot_FS!(Ham::Hamiltonian , bz::BZ , Efermi::Vector{Float64} , band_index::Vector{Int64})--> Plots.plot()
    ```
    Function to draw the fermi surface at `Efermi` for the given `Hamiltonian` on the given `BZ`. 

    """
    function plot_FS!(Ham::Hamiltonian , bz::BZ , Efermi::Vector{Float64} , band_index::Vector{Int64} ; cmp::Symbol = :turbo)
        pyplot()
        @assert length(size(Ham.bands)) == 2 "Fermi surface plots only work for 2d Hamiltonians"
        offsets     =   getAllOffsets(1, 2)

        plt = Plots.plot(aspect_ratio=:equal)

        for j in band_index
            energies    =   getindex.(Ham.bands , j)

            for offset in offsets
                shift   =   sum(offset .* bz.basis)
                ks      =   bz.ks .+ Ref(shift) 

                x       =   getindex.(ks, 1)
                y       =   getindex.(ks, 2)
                contour!(x, y, energies, levels=Efermi, color=cmp, size=(600,600), cbar = false, lw=2)

                for key in keys(bz.HighSymPoints)
                    p = Tuple(shift .+ bz.HighSymPoints[key])
                    scatter!(p, label = "", markercolor=:orange)
                    annotate!(p..., Plots.text(L"%$key", :bottom, :left, 8))
                end

                if offset == [0, 0]
                    boundary_1  =   vcat(x[begin, :], x[:, end], x[end, :], x[:, begin])
                    boundary_2  =   vcat(y[begin, :], y[:, end], y[end, :], y[:, begin])
                    scatter!(boundary_1, boundary_2, label = "", alpha=0.5, markercolor="orange", markersize=0.5)
                end

            end
        end

        for (i, b) in enumerate(bz.basis)
            Plots.plot!([0 , getindex(b , 1)] , [0 , getindex(b , 2)] , arrow=true , color=:black , lw=1 , label = L"b_{%$i}", linestyle=:dash)
        end

        xlabel!(L"k_x", guidefontsize = 9)
        ylabel!(L"k_y", guidefontsize = 9)
        Ef  =   round.(Efermi, digits=2)
        title!("Fermi Surface at " * L"E_f = %$Ef", titlefontsize = 12)

        max_range   =   maximum(norm.(bz.basis))
        xlims!((-max_range , max_range))
        ylims!((-max_range , max_range))

        return plt
    end

end