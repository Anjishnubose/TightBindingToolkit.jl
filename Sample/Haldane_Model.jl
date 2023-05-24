using TightBindingToolkit, Plots

"""
This script sets up the famous Haldane Model on a honeycomb lattice, and finds it chern numbers.
"""

"""
The lattice primitive vector for a honeycomb lattice : a1 and a2.
"""
a1 = [ 1/2 , sqrt(3)/2 ]
a2 = [-1/2 , sqrt(3)/2 ]
UC = UnitCell( [a1 , a2] , 1) ##### localDim = 1 cause its a spinless Hamiltonian

"""
The two sub-lattices of a honeycomb lattice.
"""
b1 = [ 0.0 , 0.0 ]
b2 = [ 0.0 , 1/sqrt(3) ]
addBasisSite!( UC , b1 )
addBasisSite!( UC , b2 )

"""
Adding structure to the lattice now, through the bond objects.
"""
################## On-site ###################
const M   =   1.0
addAnisotropicBond!( UC , 1 , 1 , [ 0 , 0 ]  ,  M , 0.0 , "M")
addAnisotropicBond!( UC , 2 , 2 , [ 0 , 0 ]  , -M , 0.0 , "M")
################# Nearest Neighbour hopping #################
const t1  =   -1.0
const NNdistance  =   1/sqrt(3)
addIsotropicBonds!(UC , NNdistance, t1 , "t1")
################# 2nd Nearest Neighbour #################
const t2  =   (-1.0) * im
const secondNNdistance = 1.0 
addAnisotropicBond!( UC , 1 , 1 , [ 1 , 0 ]  , t2 , secondNNdistance , "t2")
addAnisotropicBond!( UC , 1 , 1 , [ 0 , -1 ] , t2 , secondNNdistance , "t2")
addAnisotropicBond!( UC , 1 , 1 , [ -1 , 1 ] , t2 , secondNNdistance , "t2")
addAnisotropicBond!( UC , 2 , 2 , [ 0 , 1 ]  , t2 , secondNNdistance , "t2")
addAnisotropicBond!( UC , 2 , 2 , [ -1 , 0 ] , t2 , secondNNdistance , "t2")
addAnisotropicBond!( UC , 2 , 2 , [ 1 , -1 ] , t2 , secondNNdistance , "t2")

"""
Making the Brillouin Zone with kSize * kSize discrete points.
"""
const kSize   =   6 * 20 + 3   ##### a Monkhorst grid of size N = 6n+3 covers all the High-symemtry points.
bz            =   BZ(kSize)
fillBZ!(bz, UC)

"""
Spanning over different t2/M values
"""
t2s     =   range(-1.0, 1.0, 50) * M  
gaps    =   Float64[]
mus     =   Float64[]
Cherns  =   Float64[]
for (i, t2Val) in enumerate(t2s)

    if i>1
        ScaleBonds!(UC, "t2", t2Val / t2s[i-1])
    else
        ScaleBonds!(UC, "t2", t2Val / (-1))
    end
    ##### Since t2 bonds are anisotropic, the simpler function ModifyBonds! wont work.
    """
    Constructing the Hamiltonian in k-space, and diagonalizing it.
    """
    H       =   Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)

    """
    Filling the model with fermions at half-filling, with the default temperature T=0.0.
    """
    HaldaneModel       =   Model(UC, bz, H ; filling=0.5)
    SolveModel!(HaldaneModel)
    push!(mus, HaldaneModel.mu)
    push!(gaps, HaldaneModel.gap)

    """
    Then the chern numbers of the first band in the Haldane Model becomes : 
    """
    chern   =   ChernNumber(H, [1])
    push!(Cherns, chern)
end

# println("The chemical potential at half filling, μ = ", HaldaneModel.mu)
# println("The energy gap at half filling, gap = ", HaldaneModel.gap)

# println("The chern number of the lowest band is ",chern)
gapPlot     = plot(t2s, gaps, labels="gap", lw=3)
chernPlot   = plot(t2s, Cherns, labels="chern", lw=3)





