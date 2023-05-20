module TightBindingToolkit

# Write your package code here.
include("SpinMats.jl")
export SpinMats

include("UnitCell.jl")
export Bond , isSameBond , UnitCell , getDistance , addBasisSite! , addAnisotropicBond! , addIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!

include("BZ.jl")
export getRLVs , BZ , fillBZ! , ReduceQ , GetQIndex , GetIndexPath , getBZPath , CombinedBZPath

include("Hamiltonian.jl")
export Hamiltonian , FillHamiltonian , DiagonalizeHamiltonian! , DOS

include("Model.jl")
export Model , dist , distDer , findFilling , getMu! , getFilling! , getCount , getGk! , SolveModel!

include("Chern.jl")
export FindLinks , FieldStrength , ChernNumber

include("Susceptibility.jl")
export susceptibility , FindChi , FillChis!

end
