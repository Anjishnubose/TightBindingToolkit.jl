using LinearAlgebra, Plots, LaTeXStrings

include("../src/TightBindingToolkit.jl")
using .TightBindingToolkit

a1   =  [ 3/2, sqrt(3)/2]
a2   =  [-3/2, sqrt(3)/2]

b1   =  [ 0.0 , 0.0 ]
b2   =  [ 1/2 , -1/(2 * sqrt(3)) ]
b3   =  [ 1/2 , -sqrt(3)/(2) ]
b4   =  [ 0.0 , -2/sqrt(3)   ]
b5   =  [-1/2 , -sqrt(3)/(2) ]
b6   =  [-1/2 , -1/(2 * sqrt(3))]

firstNNdistance     = 1.0/sqrt(3)
secondNNdistance    = 1.0
thirdNNdistance     = 2/sqrt(3)

UC      =   UnitCell([a1, a2], 2, 2)
UC.BC   =   zeros(ComplexF64, 2)

AddBasisSite!.(Ref(UC), [b1, b2, b3, b4, b5, b6])

SpinVec     =   SpinMats(1//2)

const tIntra  =     1.0
const tInter  =     -10.0

tIntraParam     =   Param(tIntra, 2)
AddAnisotropicBond!(tIntraParam, UC, 1, 2, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 2, 3, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 3, 4, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 4, 5, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 5, 6, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")
AddAnisotropicBond!(tIntraParam, UC, 6, 1, [ 0, 0], SpinVec[4], firstNNdistance, "Intra hopping")


tInterParam     =   Param(tInter, 2)
AddAnisotropicBond!(tInterParam, UC, 1, 4, [ 1,  1], SpinVec[4], firstNNdistance, "Inter hopping")
AddAnisotropicBond!(tInterParam, UC, 3, 6, [ 0, -1], SpinVec[4], firstNNdistance, "Inter hopping")
AddAnisotropicBond!(tInterParam, UC, 5, 2, [-1,  0], SpinVec[4], firstNNdistance, "Inter hopping")

CreateUnitCell!(UC, [tIntraParam, tInterParam])


###################### k-Space #################################

UCPlot = Plot_UnitCell!(UC ; bond_cmp=:viridis, site_size=10.0, plot_conjugate=false, plot_lattice=true)

# kSize   =   6*50 + 3
# bz      =   BZ([kSize, kSize])
# FillBZ!(bz, UC)

# H       =   Hamiltonian(UC, bz)
# DiagonalizeHamiltonian!(H)

# VortexModel     =   Model(UC, bz, H ; T=0.001, mu=0.0)
# SolveModel!(VortexModel ; get_gap = true)

# Plot_FS!(H, bz, collect(0.05:0.05:0.5), [7]; cbar=true)
# Plot_Band_Structure!(VortexModel , [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; labels = [L"\Gamma", L"K_1", L"M_2"] )

###################### Real Space ######################################

lattice     =   Lattice(UC, [6, 6])
FillLattice!(lattice)

H = LatticeHamiltonian(lattice)
DiagonalizeHamiltonian!(H)

plot(H.bands, marker = :circle, label = "Energies")
Plot_Lattice!(lattice ; bond_cmp=:viridis, site_size=10.0)

# M   =   LatticeModel(lattice, H)

######################## Symmetry ##########################

# Ta1 = Translations(lattice ; primitive  = 1, with_local = true)
# Ta2 = Translations(lattice ; primitive  = 2, with_local = true)

# m1 = GetPhase.(FindQuantumNumbers(M, Ta1 ; till_band = length(H.bands)÷2))
# m2 = GetPhase.(FindQuantumNumbers(M, Ta2 ; till_band = length(H.bands)÷2))

# mGround = mod.([sum(sum.(m1)), sum(sum.(m2))], Ref(1.0))