using Plots, LinearAlgebra, LaTeXStrings#, TightBindingToolkit
include("../src/TightBindingToolkit.jl")
using .TightBindingToolkit

##### Triangle lattice
const a1  =   [ 0.5, sqrt(3)/2]
const a2  =   [-0.5, sqrt(3)/2]

const b1  =   [0.0, 0.0]

HoppingUC       =   UnitCell([a1, a2], 1, 2)
PairingUC       =   UnitCell([a1, a2], 1, 2)

AddBasisSite!(HoppingUC, b1)
AddBasisSite!(PairingUC, b1)

const t     =   1.0
##### HoppingParams
const t1    =   t
t1Param     =   Param(-t1, 2)
AddIsotropicBonds!(t1Param, HoppingUC, 1.0, 1.0, "s Hopping")

const Δ1    =   1.0 * t
f1Param     =   Param(Δ1, 2)
AddAnisotropicBond!(f1Param, PairingUC, 1, 1, [ 1,  0], -1.0 , 1.0, "f_x^3-3xy^2 Hopping")
AddAnisotropicBond!(f1Param, PairingUC, 1, 1, [ 0,  1],  1.0 , 1.0, "f_x^3-3xy^2 Hopping")
AddAnisotropicBond!(f1Param, PairingUC, 1, 1, [ 1, -1],  1.0 , 1.0, "f_x^3-3xy^2 Hopping")

const Δ2    =   1.0 * t
f2Param     =   Param(-Δ2 / (3 * sqrt(3)), 2)
AddAnisotropicBond!(f2Param, PairingUC, 1, 1, [ 1,  1],   im, sqrt(3), "f_y^3-3yx^2 Hopping")
AddAnisotropicBond!(f2Param, PairingUC, 1, 1, [ 2, -1],  -im, sqrt(3), "f_y^3-3yx^2 Hopping")
AddAnisotropicBond!(f2Param, PairingUC, 1, 1, [-1,  2],  -im, sqrt(3), "f_y^3-3yx^2 Hopping")

HoppingParams   =   [t1Param]
CreateUnitCell!(HoppingUC, HoppingParams)

PairingParams   =   [f1Param, f2Param]
CreateUnitCell!(PairingUC, PairingParams)

const n       =   20
const kSize   =   6 * n + 3  
bz            =   BZ([kSize, kSize])
FillBZ!(bz, HoppingUC)

##### Thermodynamic parameters
const T         =   0.001
const filling   =   0.5
const stat      =   -1

mus             =   collect(LinRange(-8.0, 4.0, 61))
Cherns          =   Float64[]
gaps            =   Float64[]
fillings        =   Float64[]

H       =   Hamiltonian(HoppingUC, PairingUC, bz)
DiagonalizeHamiltonian!(H)

M       =   BdGModel(HoppingUC, PairingUC, bz, H ; mu = 0.0, T = T)
SolveModel!(M)

# Plot_Band_Structure!(M ,[bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; labels=[L"$\Gamma$", L"K_1", L"M_2"], nearest = true)

for mu in mus

    """
    Filling the model with fermions at half-filling, with the default temperature T=0.0.
    """
    M       =   BdGModel(HoppingUC, PairingUC, bz, H ; mu = mu, T = T)
    SolveModel!(M ;  get_gap = true)

    push!(gaps, M.gap)
    push!(fillings, M.filling)

    """
    Then the chern numbers of the lowest two bands (1 and 2) in the Haldane Model becomes : 
    """
    chern   =   ChernNumber(H, [1])
    push!(Cherns, chern)
end
