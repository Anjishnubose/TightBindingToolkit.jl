using Plots, LaTeXStrings, TightBindingToolkit

"""
This script sets up the famous Kitaev chain in 1d.
"""

"""
The lattice primitive vector for a honeycomb lattice : a1 and a2.
"""
a1 = [ 1.0 ]
UC_hop  = UnitCell( [a1] , 1) ##### localDim = 1 cause its a spinless Hamiltonian
UC_pair = UnitCell( [a1] , 1) ##### localDim = 1 cause its a spinless Hamiltonian

"""
The single sub-lattice of the chain
"""
b1 = [ 0.0 ]
addBasisSite!( UC_hop  , b1 )
addBasisSite!( UC_pair , b1 )

"""
Adding structure to the lattice now, through the bond objects.
"""
################# Nearest Neighbour hopping #################
const t1    =   1.0
t1Param     =   Param(t1)   
const NNdistance  =   1.0
addIsotropicBonds!(t1Param, UC_hop , NNdistance, 1.0 , "t1")
################# Nearest Neighbour pairing #################
const p1    =   0.5
p1Param     =   Param(p1)   
addIsotropicBonds!(p1Param, UC_pair , NNdistance, 1.0 , "p1")

CreateUnitCell!(UC_hop , [t1Param])
CreateUnitCell!(UC_pair, [p1Param])

##### BZ with size kSize
const n       =   40
const kSize   =   6 * n + 3
bz            =   BZ(kSize, 1)
fillBZ!(bz, UC_hop)

##### Some model parameters
const T     =   0.001
const mu    =   0.0
ModifyIsotropicFields!(UC_hop, mu, 4)
H       =   Hamiltonian(UC_hop , UC_pair , bz)
DiagonalizeHamiltonian!(H)

KitaevModel   =   BdGModel(UC_hop , UC_pair, bz, H; T=T, mu=mu)
SolveModel!(KitaevModel)

kPoints = [bz.HighSymPoints["G"], bz.HighSymPoints["M1"]]

plot_band_structure!(KitaevModel, kPoints ; labels = [L"$\Gamma$", L"$M_1$"])



