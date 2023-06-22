module TightBindingToolkit

# Write your package code here.
include("SpinMats.jl")
using .SpinMatrices
export SpinMats

include("UnitCell.jl")
using .UCell
export Bond , isSameBond , UnitCell , isSameUnitCell , getDistance , addBasisSite! , addAnisotropicBond! , addIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!, Lookup

include("Params.jl")
using .Parameters
export Param, addAnisotropicBond!, addIsotropicBonds!, CreateUnitCell!, ModifyUnitCell!, GetParams

include("BZ.jl")
using .BZone
export getRLVs , BZ , fillBZ! , Monkhorst, ReduceQ , GetQIndex , GetIndexPath , CombinedIndexPath , getBZPath , CombinedBZPath

include("Hamiltonian.jl")
using .Hams
export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , DiagonalizeHamiltonian! , DOS, ModifyHamiltonianField!

include("Model.jl")
using .TBModel
export Model , dist , distDer , findFilling , getMu! , getFilling! , getCount , getGk! , SolveModel! , BinarySearch

include("BdGModel.jl")
using .BdG
export BdGModel, findFilling, getMu!, getFilling!, getGk!, SolveModel!

include("Chern.jl")
using .Chern
export FindLinks , FieldStrength , ChernNumber

include("Susceptibility.jl")
using .suscep
export susceptibility , FindChi , FillChis!

include("Plot.jl")
using .PlotTB
export plot_UnitCell! , plot_band_contour!, plot_band_structure!, plot_FS!

end
