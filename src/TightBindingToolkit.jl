module TightBindingToolkit

# Write your package code here.
include("SpinMats.jl")
using .SpinMatrices
export SpinMats

include("UnitCell.jl")
using .UCell
export Bond , isSameBond , UnitCell , getDistance , addBasisSite! , addAnisotropicBond! , addIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!

include("Params.jl")
using .Parameters
export Param, addAnisotropicBond!, addIsotropicBonds!, CreateUnitCell!, ModifyUnitCell!, GetParams

include("BZ.jl")
using .BZone
export getRLVs , BZ , fillBZ! , Monkhorst, ReduceQ , GetQIndex , GetIndexPath , getBZPath , CombinedBZPath

include("Hamiltonian.jl")
using .Hams
export Hamiltonian , FillHamiltonian , FullHamiltonian , DiagonalizeHamiltonian! , DOS

include("Model.jl")
using .TBModel
export Model , dist , distDer , findFilling , getMu! , getFilling! , getCount , getGk! , SolveModel!

include("Chern.jl")
using .Chern
export FindLinks , FieldStrength , ChernNumber

include("Susceptibility.jl")
using .suscep
export susceptibility , FindChi , FillChis!

end
