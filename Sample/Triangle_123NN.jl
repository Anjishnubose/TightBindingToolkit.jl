using Plots
using TightBindingToolkit

"""
This script sets up a simple triangular lattice with first, second, and third neighbour hoppings.
It then calculates the magnetic susceptibility at zero energy, along a path in the brillouin zone through the Γ, K, and M points.
"""

"""
The lattice primitive vector for a triangular lattice : a1 and a2.
"""
a1 = [ 1/2 , sqrt(3)/2 ]
a2 = [-1/2 , sqrt(3)/2 ]
UC = UnitCell( [a1 , a2] , 2) ##### localDim=2 since we are working with spin-1/2 particles now

"""
Unit cell has only one sub-lattice.
"""
b1 = [ 0.0 , 0.0 ]
addBasisSite!( UC , b1 )

"""
Adding structure to the lattice now, through the bond objects.
"""
SpinVec     =   SpinMats(1//2) ##### Working with spin-1/2
################ Nearest Neighbour #################
const t1    =   -1.0
t1Param     =   Param(t1)   
const NNdistance  =   1.0
addIsotropicBonds!(t1Param, UC , NNdistance, SpinVec[4] , "t1")
############### 2nd Nearest Neighbour #################
const t2    =   0.25
t2Param     =   Param(t2)  
const secondNNdistance  =   sqrt(3)
addIsotropicBonds!(t2Param, UC , secondNNdistance, SpinVec[4] , "t2" ; checkOffsetRange = 2)
############### 3rd Nearest Neighbour #################
const t3    =   -0.25
t3Param     =   Param(t3)   
const thirdNNdistance  =   2.0
addIsotropicBonds!(t3Param, UC , thirdNNdistance, SpinVec[4] , "t3" ; checkOffsetRange = 3)

CreateUnitCell!(UC, [t1Param, t2Param, t3Param])

"""
Making the Brillouin Zone with kSize * kSize discrete points.
"""
const n       =   50
const kSize   =   6 * n + 3 ##### a Monkhorst grid of size N = 6n+3 covers all the High-symemtry points.
bz            =   BZ(kSize)
fillBZ!(bz, UC)

"""
A path in the discretized BZ which goes through Γ -> K -> M -> Γ.
"""
path = CombinedBZPath(bz, [bz.HighSymPoints["Gamma"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; nearest=true)
##### The keyword nearest=true since the Monkhorst grid can never actually have the M points in it, so you have to find the point on the grid closest to them.

"""
Some model parameters : temperature and chemical potential
"""
const T     =   0.001
const mu    =   0.0   
H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

"""
Filling the model at this chemical potential with spin-1/2 fermions.
"""
TriangleModel   =   Model(UC, bz, H; T=T, mu=mu)
SolveModel!(TriangleModel)
println("The filling at μ = 0 is ", TriangleModel.filling)

Omegas  =   [0.0]
const spread  = 1e-3  
dos = DOS(Omegas, H)    ##### the single particle density of states

suscep  =   susceptibility(TriangleModel, path, Omegas ; eta=spread)
FillChis!(suscep)
##### Plotting the real part of the magnetic susceptibility at 0 energy, and momentum points along the path.
plot(real.(suscep.chis["zz"][1, :]), labels="χ_zz(Ω=0, Q)", lw=2.0, markershape=:circle, markercolor="red")



