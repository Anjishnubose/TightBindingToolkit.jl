using Plots, LaTeXStrings, TightBindingToolkit

a1  = [ 1/2 , sqrt(3)/2 ]
a2  = [-1/2 , sqrt(3)/2 ]

b1  = [ 0.0 , 0.0 ]
b2  = [ 0.0 , 1.0/sqrt(3) ]

const firstNNdistance  = 1.0/sqrt(3)
const secondNNdistance = 1.0
const thirdNNdistance  = 2/sqrt(3)

#### Parton Lattice 
UC = UnitCell( [a1 , a2] , 2 , 2)

AddBasisSite!( UC , b1 )
AddBasisSite!( UC , b2 )

#################### Bonds #########################################
SpinVec     =   SpinMats(1//2)
######## NN #################
const t1  =   -1.0
AddIsotropicBonds!(UC , firstNNdistance, t1 * SpinVec[4] , "t1")
######## 2NN #################
const t2  =   -0.00 * im
AddAnisotropicBond!( UC , 1 , 1 , [ 1 , 0 ]  , t2 * SpinVec[3] , secondNNdistance , "t2")
AddAnisotropicBond!( UC , 1 , 1 , [ 0 , -1 ] , t2 * SpinVec[3] , secondNNdistance , "t2")
AddAnisotropicBond!( UC , 1 , 1 , [ -1 , 1 ] , t2 * SpinVec[3] , secondNNdistance , "t2")
AddAnisotropicBond!( UC , 2 , 2 , [ 0 , 1 ]  , t2 * SpinVec[3] , secondNNdistance , "t2")
AddAnisotropicBond!( UC , 2 , 2 , [ -1 , 0 ] , t2 * SpinVec[3] , secondNNdistance , "t2")
AddAnisotropicBond!( UC , 2 , 2 , [ 1 , -1 ] , t2 * SpinVec[3] , secondNNdistance , "t2")
######## 3NN #################
const t3  =   -0.00
AddIsotropicBonds!(UC , thirdNNdistance, t3 * SpinVec[3] , "t3", checkOffsetRange = 3)
# p       =   plot_UnitCell!(UC ; range=1)

##### BZ with size kSize
const n       =   100
const kSize   =   6 * n 
bz            =   BZ(kSize)
FillBZ!(bz, UC)

H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

Omegas  =   collect(range(H.bandwidth..., 501))
dos     =   DOS(Omegas, H)




