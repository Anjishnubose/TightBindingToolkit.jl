using LinearAlgebra, Plots, LaTeXStrings

include("../src/TightBindingToolkit.jl")
using .TightBindingToolkit

a1   =  [ 3/2, sqrt(3)/2]
a2   =  [-3/2, sqrt(3)/2]

b1   =  [ 0.0 , 0.0 ]
b2   =  [ 0.0 , 1.0/sqrt(3) ]
b3   =  [-1.0 , 1.0/sqrt(3) ]
b4   =  [-1/2 , sqrt(3)/2   ]
b5   =  [ 1.0 , 1.0/sqrt(3) ]
b6   =  [ 1/2 , sqrt(3)/2   ]

firstNNdistance     = 1.0/sqrt(3)
secondNNdistance    = 1.0
thirdNNdistance     = 2/sqrt(3)

UC      =   UnitCell([a1, a2], 2, 2)

AddBasisSite!.(Ref(UC), [b1, b2, b3, b4, b5, b6])

SpinVec     =   SpinMats(1//2)

const tIntra  =     1.0
const tInter  =     -1.0

tIntraParam     =   Param(tIntra, 2)
AddAnisotropicBond!(tIntraParam, UC, 3, 1, [ 0,  1], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 1, 5, [-1,  0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 4, 2, [ 0,  0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 2, 6, [ 0,  0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 5, 4, [ 0, -1], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 6, 3, [ 1,  0], SpinVec[4], firstNNdistance, "Intra hopping")


tInterParam     =   Param(tInter, 2)
AddAnisotropicBond!(tInterParam, UC, 2, 1, [ 0,  0], SpinVec[4], firstNNdistance, "Inter hopping")
AddAnisotropicBond!(tInterParam, UC, 3, 4, [ 0,  0], SpinVec[4], firstNNdistance, "Inter hopping")
AddAnisotropicBond!(tInterParam, UC, 5, 6, [ 0,  0], SpinVec[4], firstNNdistance, "Inter hopping")

CreateUnitCell!(UC, [tIntraParam, tInterParam])

Plot_UnitCell!(UC ; bond_cmp=:viridis, site_size=10.0, plot_conjugate=false, plot_lattice=true)

kSize   =   6*50 + 3
bz      =   BZ([kSize, kSize])
FillBZ!(bz, UC)

H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

VortexModel     =   Model(UC, bz, H ; T=0.001, mu=0.0)
SolveModel!(VortexModel ; get_gap = true)

# Plot_FS!(H, bz, collect(0.05:0.05:0.5), [7]; cbar=true)

Plot_Band_Structure!(VortexModel , [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; labels = [L"\Gamma", L"K_1", L"M_2"] )