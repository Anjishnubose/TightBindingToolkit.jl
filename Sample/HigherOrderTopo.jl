using Plots, LinearAlgebra, LaTeXStrings, TightBindingToolkit

##### Triangle lattice
const a1  =   [ 0.5, sqrt(3)/2]
const a2  =   [-0.5, sqrt(3)/2]

const b1  =   [0.0, 0.0]

HoppingUC       =   UnitCell([a1, a2], 2, 2)

AddBasisSite!(HoppingUC, b1)

SpinVec     =   SpinMats(1//2)
const t     =   1.0
##### HoppingParams
const t1    =   t
t1Param     =   Param(-t1, 2)
AddIsotropicBonds!(t1Param, HoppingUC, 1.0, SpinVec[3], "s Hopping")

const Δ1    =   1.5 * t
f1Param     =   Param(Δ1, 2)
AddAnisotropicBond!(f1Param, HoppingUC, 1, 1, [ 1,  0],  -im * SpinVec[2], 1.0, "f_x^3-3xy^2 Hopping")
AddAnisotropicBond!(f1Param, HoppingUC, 1, 1, [ 0,  1],   im * SpinVec[2], 1.0, "f_x^3-3xy^2 Hopping")
AddAnisotropicBond!(f1Param, HoppingUC, 1, 1, [ 1, -1],   im * SpinVec[2], 1.0, "f_x^3-3xy^2 Hopping")

const Δ2    =   1.5 * t
f2Param     =   Param(-Δ2 / (3 * sqrt(3)), 2)
AddAnisotropicBond!(f2Param, HoppingUC, 1, 1, [ 1,  1],   im * SpinVec[1], sqrt(3), "f_y^3-3yx^2 Hopping")
AddAnisotropicBond!(f2Param, HoppingUC, 1, 1, [ 2, -1],  -im * SpinVec[1], sqrt(3), "f_y^3-3yx^2 Hopping")
AddAnisotropicBond!(f2Param, HoppingUC, 1, 1, [-1,  2],  -im * SpinVec[1], sqrt(3), "f_y^3-3yx^2 Hopping")

const δ     =   0.5 * t
dParam      =   Param((3 * t - δ), 2)
AddIsotropicBonds!(dParam, HoppingUC, 0.0, SpinVec[3], "Oribital gap")

HoppingParams   =   [t1Param, f1Param, f2Param, dParam]
CreateUnitCell!(HoppingUC, HoppingParams)


const n       =   20
const kSize   =   6 * n + 3  
bz            =   BZ([kSize, kSize])
FillBZ!(bz, HoppingUC)

##### Thermodynamic parameters
const T         =   0.001
const filling   =   0.5
const stat      =   -1

deltas          =   collect(LinRange(-0.5, 0.5, 51))
Cherns          =   Float64[]

# H       =   Hamiltonian(HoppingUC, bz)
# DiagonalizeHamiltonian!(H)

# M       =   Model(HoppingUC, bz, H ; filling=0.5)
# SolveModel!(M)

# Plot_Band_Structure!(M ,[bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]] ; labels=[L"$\Gamma$", L"K_1", L"M_2"], nearest = true)

for delta in deltas

    push!(dParam.value, (3 * t - delta))
    ModifyUnitCell!(HoppingUC, [dParam])

    """
    Constructing the Hamiltonian in k-space, and diagonalizing it.
    """
    H       =   Hamiltonian(HoppingUC, bz)
    DiagonalizeHamiltonian!(H)

    """
    Filling the model with fermions at half-filling, with the default temperature T=0.0.
    """
    M       =   Model(HoppingUC, bz, H ; filling=0.5)
    SolveModel!(M)

    """
    Then the chern numbers of the lowest two bands (1 and 2) in the Haldane Model becomes : 
    """
    chern   =   ChernNumber(H, [1])
    push!(Cherns, chern)
end
