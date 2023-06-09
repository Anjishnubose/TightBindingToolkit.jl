module TightBindingToolkit

####### Write your package code here.

##### Module to get SU(2) representation for any half-integer / integer spin.
include("SpinMats.jl")
using .SpinMatrices
export SpinMats

include("Useful.jl")
using .Useful
export GetAllOffsets, VecAngle, Meshgrid, BinarySearch, DistFunction, DeriDistFunction , GetIndexPath, FFTArrayofMatrix

##### Module for defining Unit Cell structure and some funcitons to manipulate it.
include("UnitCell.jl")
using .UCell
export Bond , BondRank, IsSameBond , FlipBond , UnitCell , IsSameUnitCell , AddBasisSite! , GetDistance , GetRealSpacePositions, IsSameUnitCell

include("DesignUnitCell.jl")
using .DesignUCell
export AddAnisotropicBond! , AddIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!, Lookup

include("ExpandUnitCell.jl")
using .ExpandUCell
export ChangePrimitives!, ExpandUnitCell, ExpandBonds!

##### Module to define Param structure which represents a general hopping parameter, and can be used to construct a UnitCell.
include("Params.jl")
using .Parameters
export Param, AddAnisotropicBond!, AddIsotropicBonds!, CreateUnitCell!, ModifyUnitCell!, GetParams

##### Module to define Brillouin Zone structure and some functions to manipulate it.
include("BZ.jl")
using .BZone
export GetRLVs , BZ , Monkhorst , FillBZ! , ReduceQ , GetQIndex , CombinedIndexPath , GetBZPath , CombinedBZPath , MomentumPhaseFFT

##### module to define a Hamiltonian structure and solving for some of its properties.
include("Hamiltonian.jl")
using .Hams
export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , DiagonalizeHamiltonian! , DOS, ModifyHamiltonianField!, IsBandGapped

##### Module to define a Tight-Binding Model structure which takes into account thermodynamical parameters such as temperature and filling etc.
include("Model.jl")
using .TBModel
export Model , FindFilling , GetMu! , GetFilling! , GetCount , GetGk! , GetGr!, SolveModel!

##### Module to define the equivalent but for bdG systems with pairing.
include("BdGModel.jl")
using .BdG
export BdGModel, FindFilling , GetMu! , GetFilling! , GetGk! , GetGr!, SolveModel!

##### Module to calculate generalized Chern numbers
include("Chern.jl")
using .Chern
export FindLinks , FieldStrength , ChernNumber, CheckValidity

##### Module to calculate generalized magnetic susceptibility
include("Susceptibility.jl")
using .suscep
export Susceptibility , FindChi , FillChis!

##### Module containing some plotting functions
include("Plot.jl")
using .PlotTB
export Plot_UnitCell! , Plot_Band_Contour!, Plot_Band_Structure!, Plot_FS!

end
