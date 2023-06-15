using Plots
using TightBindingToolkit

"""
This script sets up the famous Haldane Model on a honeycomb lattice, and finds it chern numbers.
"""

"""
The lattice primitive vector for a honeycomb lattice : a1 and a2.
"""
a1 = [ 1/2 , sqrt(3)/2 ]
a2 = [-1/2 , sqrt(3)/2 ]
UC = UnitCell( [a1 , a2] , 2) ##### localDim = 2 cause we are Constructing a spinful version of the Haldane Model.

"""
The two sub-lattices of a honeycomb lattice.
"""
b1 = [ 0.0 , 0.0 ]
b2 = [ 0.0 , 1/sqrt(3) ]
addBasisSite!( UC , b1 )
addBasisSite!( UC , b2 )

SpinVec     =   SpinMats(1//2)

"""
Adding structure to the lattice now, through the bond objects.
"""
################## On-site ###################
const M     =   1.0
MParam      =   Param(M)   
addAnisotropicBond!( MParam, UC , 1 , 1 , [ 0 , 0 ]  ,  SpinVec[4], 0.0 , "M")
addAnisotropicBond!( MParam, UC , 2 , 2 , [ 0 , 0 ]  , -SpinVec[4], 0.0 , "M")
################# Nearest Neighbour hopping #################
const t1    =   -1.0
t1Param     =   Param(t1)   
const NNdistance  =   1/sqrt(3)
addIsotropicBonds!(t1Param, UC , NNdistance, SpinVec[4] , "t1")
################# 2nd Nearest Neighbour #################
const t2    =   (-1.0) 
t2Param     =    Param(t2)  
const secondNNdistance = 1.0 
addAnisotropicBond!(t2Param, UC , 1 , 1 , [ 1 , 0 ]  , im * SpinVec[4], secondNNdistance , "t2")
addAnisotropicBond!(t2Param, UC , 1 , 1 , [ 0 , -1 ] , im * SpinVec[4], secondNNdistance , "t2")
addAnisotropicBond!(t2Param, UC , 1 , 1 , [ -1 , 1 ] , im * SpinVec[4], secondNNdistance , "t2")
addAnisotropicBond!(t2Param, UC , 2 , 2 , [ 0 , 1 ]  , im * SpinVec[4], secondNNdistance , "t2")
addAnisotropicBond!(t2Param, UC , 2 , 2 , [ -1 , 0 ] , im * SpinVec[4], secondNNdistance , "t2")
addAnisotropicBond!(t2Param, UC , 2 , 2 , [ 1 , -1 ] , im * SpinVec[4], secondNNdistance , "t2")

CreateUnitCell!(UC, [MParam, t1Param, t2Param])
"""
Making the Brillouin Zone with kSize * kSize discrete points.
"""
const kSize   =   6 * 20 + 3   ##### a Monkhorst grid of size N = 6n+3 covers all the High-symemtry points.
bz            =   BZ(kSize)    ##### Make the BZ explicitly in 2d
fillBZ!(bz, UC)

"""
Spanning over different t2/M values
"""
t2s     =   range(-1.0, 1.0, 11) * M  
gaps    =   Float64[]
mus     =   Float64[]
Cherns  =   Float64[]

for t2Val in t2s

    push!(t2Param.value, t2Val)
    ModifyUnitCell!(UC, [t2Param])

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
    Then the chern numbers of the lowest two bands (1 and 2) in the Haldane Model becomes : 
    """
    chern   =   ChernNumber(H, [1, 2])
    push!(Cherns, chern)
end

chernPlot   = plot(t2s, Cherns, labels="chern", lw=3, markershape = :circle )





