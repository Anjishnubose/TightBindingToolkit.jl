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
UC.BC   =   [0.0, 0.0]

AddBasisSite!.(Ref(UC), [b1, b2, b3, b4, b5, b6])

SpinVec     =   SpinMats(1//2)

const tIntra  =     1.0
const tInter  =     -10.0
const tInter  =     -1.5

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

kSize   =   6*50 + 3
bz      =   BZ([kSize, kSize])
FillBZ!(bz, UC)

H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

VortexModel     =   Model(UC, bz, H ; T=0.001, mu=0.0)
SolveModel!(VortexModel ; get_gap = true)

R3      =   zeros(Float64, 6, 6)
R3[3, 1]=   1.0
R3[4, 2]=   1.0
R3[5, 3]=   1.0
R3[6, 4]=   1.0
R3[1, 5]=   1.0
R3[2, 6]=   1.0
R3      =   kron(R3, Matrix(I, 2, 2))

R2      =   zeros(Float64, 6, 6)
R2[4, 1]=   1.0
R2[5, 2]=   1.0
R2[6, 3]=   1.0
R2[1, 4]=   1.0
R2[2, 5]=   1.0
R2[3, 6]=   1.0
R2      =   kron(R2, Matrix(I, 2, 2))

Gammaind=   GetQIndex(bz.HighSymPoints["G"], bz)
K1ind   =   GetQIndex(bz.HighSymPoints["K1"], bz)
K2ind   =   GetQIndex(bz.HighSymPoints["K2"], bz)

OrbitalsGamma   =   H.states[Gammaind...]
OrbitalsK1  =   H.states[K1ind...]
OrbitalsK2  =   H.states[K2ind...]

REigGamma   =   eigen((OrbitalsGamma' * R3 * OrbitalsGamma)[1:6, 1:6]).values
REigK1      =   eigen((OrbitalsK1' * R3 * OrbitalsK1)[1:6, 1:6]).values
REigK2      =   eigen((OrbitalsK2' * R3 * OrbitalsK2)[1:6, 1:6]).values



# Plot_FS!(H, bz, collect(0.05:0.05:0.5), [7]; cbar=true)
# Plot_Band_Structure!(VortexModel , [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; labels = [L"\Gamma", L"K_1", L"M_2"] )

###################### Real Space ######################################

# lattice     =   Lattice(UC, [6, 6])
# FillLattice!(lattice)

# # Plot_Lattice!(lattice ; bond_cmp=:viridis, site_size=8.0)

# H = LatticeHamiltonian(lattice)
# DiagonalizeHamiltonian!(H)

# plot(H.bands, marker = :circle, label = "Energies")
Plot_Lattice!(lattice ; bond_cmp=:viridis, site_size=10.0)

# M   =   LatticeModel(lattice, H)

# ######################## Symmetry ##########################

# Ta1 = Translations(lattice ; primitive  = 1, with_local = true)
# Ta2 = Translations(lattice ; primitive  = 2, with_local = true)
# Ta1 = Translations(lattice ; primitive  = 1, with_local = true)
# Ta2 = Translations(lattice ; primitive  = 2, with_local = true)

# m1 = GetPhase.(FindQuantumNumbers(M, Ta1 ; till_band = length(H.bands)÷2))
# m2 = GetPhase.(FindQuantumNumbers(M, Ta2 ; till_band = length(H.bands)÷2))
# m1 = GetPhase.(FindQuantumNumbers(M, Ta1 ; till_band = length(H.bands)÷2))
# m2 = GetPhase.(FindQuantumNumbers(M, Ta2 ; till_band = length(H.bands)÷2))

# mGround = mod.([sum(sum.(m1)), sum(sum.(m2))], Ref(1.0))
# mGround = mod.([sum(sum.(m1)), sum(sum.(m2))], Ref(1.0))