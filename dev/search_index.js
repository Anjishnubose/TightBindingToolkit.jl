var documenterSearchIndex = {"docs":
[{"location":"Model/#TightBindingToolkit.Model","page":"Tight Binding Model","title":"TightBindingToolkit.Model","text":"","category":"section"},{"location":"Model/","page":"Tight Binding Model","title":"Tight Binding Model","text":"Modules = [TightBindingToolkit, TightBindingToolkit.TBModel]\nPrivate = false\nPages   = [\"Model.jl\"]","category":"page"},{"location":"Model/#TightBindingToolkit.TBModel.Model","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.Model","text":"Model is a data type representing a general Tight Binding system.\n\nAttributes\n\nuc      ::  UnitCell: the Unit cell of the lattice.\nbz      ::  BZ: The discretized Brillouin Zone.\nHam     ::  Hamiltonian: the Hamiltonian at all momentum-points.\nT       ::  Float64: the temperature of the system.\nfilling ::  Float64: The filling of the system.\nmu      ::  Float64: The chemical potential of the system.\nstat    ::  Int64 : ±1 for bosons and fermions.\ngap     ::  Float64 : the energy gap of excitations at the given filling.\nGk      ::  Matrix{Matrix{ComplexF64}} : A matrix (corresponding to the matrix of k-points in BZ) of Greens functions.\n\nInitialize this structure using \n\nModel(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=0.0, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)\n\nYou can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.\n\n\n\n\n\n","category":"type"},{"location":"Model/#TightBindingToolkit.TBModel.SolveModel!-Tuple{Model}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.SolveModel!","text":"SolveModel!(M::Model)\n\none-step function to find all the attributes in Model after it has been initialized.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.dist-Tuple{Any}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.dist","text":"dist(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)\n\nDistribution function. stat=1 –> Bose-Einstein distribution, and stat=-1 –> Fermi distribution.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.distDer-Tuple{Any}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.distDer","text":"distDer(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)\n\nderivative of dist w.r.t the energy. stat=1 –> Bose-Einstein distribution, and stat=-1 –> Fermi distribution.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.findFilling","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.findFilling","text":"findFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1) --> Float64\n\nFind filling at given T=temperature and mu=chemical potential, for a given bands.\n\n\n\n\n\n","category":"function"},{"location":"Model/#TightBindingToolkit.TBModel.getCount-Tuple{Vector{Float64}, Float64, Float64, Int64}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.getCount","text":"getCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) --> Matrix{Float64}\n\nFunction to return a diagonal matrix whose entries are M[i, i] = θ(-(E^i(k)-μ)) ––> 1 if the energy is below the chemical potential, otherwise 0.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.getFilling!-Tuple{Model}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.getFilling!","text":"getFilling!(M::Model)\n\nFind filling for a given Model.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.getGk!-Tuple{Model}","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.getGk!","text":"getGk!(M::Model)\n\nFinding the equal-time Greens functions in momentum space of a Model.\n\n\n\n\n\n","category":"method"},{"location":"Model/#TightBindingToolkit.TBModel.getMu!","page":"Tight Binding Model","title":"TightBindingToolkit.TBModel.getMu!","text":"getMu!(M::Model, tol::Float64=0.001)\n\nFunction to get chemical potential for a given Model, within a tolerance.\n\n\n\n\n\n","category":"function"},{"location":"Params/#TightBindingToolkit.Params","page":"Parameters","title":"TightBindingToolkit.Params","text":"","category":"section"},{"location":"Params/","page":"Parameters","title":"Parameters","text":"Modules = [TightBindingToolkit, TightBindingToolkit.Parameters]\nPrivate = false\nPages = [\"Params.jl\"]","category":"page"},{"location":"Params/#TightBindingToolkit.Parameters.Param","page":"Parameters","title":"TightBindingToolkit.Parameters.Param","text":"`Param` is a data type representing a general tight-binding parameter, which can span multiple bonds.\n\n# Attributes\n- `value        ::  Vector{ Float64 }`: the strength of the parameter (or even the full history of it if its changed).\n- `unitBonds    ::  Vector{ Bond }`: All the bonds this parameter lives on. These bonds are supposed to have \"unit\" strength, and ultimately get scaled by the `value` when making the `UnitCell`.\n- `label        ::  String`: some string label to mark the parameter.  \n- `dist         ::  Float64`: the distance of the bonds the parameter lives on.\n\nInitialize this structure using \n```julia\nParam( value::Float64 )\n```\n\n\n\n\n\n","category":"type"},{"location":"Params/#TightBindingToolkit.Parameters.CreateUnitCell!","page":"Parameters","title":"TightBindingToolkit.Parameters.CreateUnitCell!","text":"CreateUnitCell!(uc::UnitCell, param::Param , index::Int64=length(param.value))\nCreateUnitCell!(uc::UnitCell, params::Vector{Param}, indices::Vector{Int64}=length.(getproperty.(params, :value)))\n\nAdd bonds corrsponding to a param to UnitCell, scaled with the param.value[index]. Also includes the broadcasted call.\n\n\n\n\n\n","category":"function"},{"location":"Params/#TightBindingToolkit.Parameters.GetParams-Tuple{UnitCell}","page":"Parameters","title":"TightBindingToolkit.Parameters.GetParams","text":"```julia\n\tGetParams(uc::UnitCell) --> Vector{Param}\n\t```\nFor legacy purposes. \n\tIf you have a `UnitCell` built using the old technique of adding bonds directly, you can get a vector of `Param` using this function, corresponding to each unique bond type already present in `UnitCell`.\n\n\n\n\n\n","category":"method"},{"location":"Params/#TightBindingToolkit.Parameters.ModifyUnitCell!-Tuple{UnitCell, Param}","page":"Parameters","title":"TightBindingToolkit.Parameters.ModifyUnitCell!","text":"ModifyUnitCell!(uc::UnitCell, param::Param)\nModifyUnitCell!(uc::UnitCell, params::Vector{Param})\n\nModify all bonds in UnitCell corresponding to given param, taking the latest value in param.value. \n\n\n\n\n\n","category":"method"},{"location":"Params/#TightBindingToolkit.UCell.addAnisotropicBond!-Tuple{Param, UnitCell, Int64, Int64, Vector{Int64}, Number, Float64, String}","page":"Parameters","title":"TightBindingToolkit.UCell.addAnisotropicBond!","text":"addAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )\naddAnisotropicBond!( param::Param, uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )\n\nAdd a bond with the given attributes to param. If given mat attribute is a number, it is converted into a 1x1 matrix when entered into the bond.\n\n\n\n\n\n","category":"method"},{"location":"Params/#TightBindingToolkit.UCell.addIsotropicBonds!-Tuple{Param, UnitCell, Float64, Number, String}","page":"Parameters","title":"TightBindingToolkit.UCell.addIsotropicBonds!","text":"```julia\naddIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))\naddIsotropicBonds!( param::Param, uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )\n```\nAdd a set of \"isotropic\" bonds, which are the same for each pair of sites at the given distance. \nIf given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.\nThe input `checkOffsetRange` must be adjusted depending on the input distance. \nThe optional input `subs` is meant for isotropic bonds when only a subset of sublattices are involved.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZ","page":"Brillouin Zone","title":"TightBindingToolkit.BZ","text":"","category":"section"},{"location":"BZ/","page":"Brillouin Zone","title":"Brillouin Zone","text":"Modules = [TightBindingToolkit, TightBindingToolkit.BZone]\nPrivate = false\nPages   = [\"BZ.jl\"]","category":"page"},{"location":"BZ/#TightBindingToolkit.BZone.BZ","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.BZ","text":"BZ is a data type representing a discretized Brillouin Zone in momentum space.\n\nAttributes\n\nbasis           :: Vector{ Vector{ Float64 } }: reciprocal lattice vectors of the Brillouin Zone.\ngridSize        :: Int64: The number of points along each dimension of the grid.\nkInds           :: Vector{Matrix{Float64}}: the Monkhorst grid corresponding to k along bs.\nks   \t       :: Matrix{Vector{Float64}}: The grid of momentum points [kx, ky].\nHighSymPoints   :: Dict: A dictionary containing the HIgh-Symmetry points Γ, K(2), and M(3).\nshift           :: Vector{Int64} : how shifted the grid is from its centre point at the Γ point, in units of 1/gridSize.\n\nInitialize this structure using \n\nBZ(gridSize::Int64)\n\n\n\n\n\n","category":"type"},{"location":"BZ/#TightBindingToolkit.BZone.CombinedBZPath-Tuple{BZ, Vector{Vector{Float64}}}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.CombinedBZPath","text":"CombinedBZPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) --> Vector{Vector{Float64}}\n\nReturns a path in momentum-space of the discretized BZ which joins the given momentum points present in points as point[1] –> point[2] –> ... –> point[end] –> point[1]. The optional input nearest is the same as in GetQIndex, and closed determines if the path is a clsoed loop or not.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.GetIndexPath-Tuple{Vector{Int64}, Vector{Int64}}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.GetIndexPath","text":"GetIndexPath(start::Vector{Int64}, ending::Vector{Int64} ; exclusive::Bool=true) --> Vector{Vector{Int64}}\n\nReturns a path in index-space of the discretized BZ which joins the two sets of indices start and ending.  If the input exclusive is set to true, the returned path will NOT contain the ending point itself.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.GetQIndex-Tuple{Vector{Float64}, BZ}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.GetQIndex","text":"GetQIndex(Q::Vector{Float64}, bz::BZ ; nearest::Bool = false) --> Vector{Int64}\n\nReturns the index in the discretized BZ of the momentum point corresponding to the fiven momentum Q.  If the input nearest is set to true, will return the index of the momentum point on the grid closest to Q, if Q does not exist on the grid. \n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.Monkhorst-Tuple{Int64, Int64}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.Monkhorst","text":"Monkhorst(ind::Int64, N::Int64) --> Float64\nMonkhorst(ind::Int64, N::Int64, shift::Int64, BC::Float64) --> Float64\n\nThe usual Monkhorst grid is defined as follows frac2i - (N+1)2N i1 N The modified Monkhorst grid takes into account the desired boundary condition BC, and an integer shift (if required, to change the starting point), and shifts the momentum grid accordingly.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.ReduceQ-Tuple{Vector{Float64}, BZ}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.ReduceQ","text":"ReduceQ(Q::Vector{Float64}, bz::BZ) --> Vector{Float64}\n\nReduces a given momentum back to the range covered by the discretized Brillouin Zone.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.fillBZ!-Tuple{BZ, UnitCell}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.fillBZ!","text":"fillBZ!(bz::BZ, uc::UnitCell, offsetRange::Int64=1 ; shift::Vector{Float64}=zeros(Float64, length(uc.primitives)))\n\nFills the BZ object with the relevant attributes, after it has been initialized as BZ(gridsize=N).\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.getBZPath-Tuple{BZ, Vector{Float64}, Vector{Float64}}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.getBZPath","text":"getBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) --> Vector{Vector{Float64}}\n\nReturns the actual path in momentum-space of the discretized BZ which joins the two momentums start and ending.  The optional input nearest is the same as in GetQIndex, and exclusive is the same as in GetIndexPath.\n\n\n\n\n\n","category":"method"},{"location":"BZ/#TightBindingToolkit.BZone.getRLVs-Tuple{UnitCell}","page":"Brillouin Zone","title":"TightBindingToolkit.BZone.getRLVs","text":"getRLVs( uc::UnitCell ) --> Vector{Vector{Float64}}\n\nReturns the reciprocal lattice vectors corresponding to the given Unit Cell.\n\n\n\n\n\n","category":"method"},{"location":"Hamiltonian/#TightBindingToolkit.Hamiltonian","page":"Hamiltonian","title":"TightBindingToolkit.Hamiltonian","text":"","category":"section"},{"location":"Hamiltonian/","page":"Hamiltonian","title":"Hamiltonian","text":"Modules = [TightBindingToolkit, TightBindingToolkit.Hams]\nPrivate = false\nPages   = [\"Hamiltonian.jl\"]","category":"page"},{"location":"Hamiltonian/#TightBindingToolkit.Hams.Hamiltonian","page":"Hamiltonian","title":"TightBindingToolkit.Hams.Hamiltonian","text":"Hamiltonian is a data type representing a general momentum-space Hamiltonian corresponding to the given UnitCell and BZ.\n\nAttributes\n\nH           :: Matrix{Matrix{ComplexF64}}: A matrix (corresponding to the matrix of k-points in BZ) of Hamiltonian matrices.\nbands       :: Matrix{Vector{Float64}}: A matrix (corresponding to the matrix of k-points in BZ) of band spectrums.\nstates      :: Matrix{Matrix{ComplexF64}}: A matrix (corresponding to the matrix of k-points in BZ) of band wavefunctions.\nbandwidth   :: Tuple{Float64, Float64} : the tuple of minimum and maximum energies in the band structure.\n\nInitialize this structure using\n\nHamiltonian(uc::UnitCell, bz::BZ)\n\n\n\n\n\n","category":"type"},{"location":"Hamiltonian/#TightBindingToolkit.Hams.DOS-Tuple{Float64, Vector{Float64}}","page":"Hamiltonian","title":"TightBindingToolkit.Hams.DOS","text":"DOS(Omega::Float64, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Float64\nDOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Vector{Float64}\n\nCalculate the Density of State correspondingto the given energies in Omegas, for the lowest bands upto till_band. The calculation is done at a finite spread of the delta-function sum. \n\n\n\n\n\n","category":"method"},{"location":"Hamiltonian/#TightBindingToolkit.Hams.DiagonalizeHamiltonian!-Tuple{Hamiltonian}","page":"Hamiltonian","title":"TightBindingToolkit.Hams.DiagonalizeHamiltonian!","text":"DiagonalizeHamiltonian!(Ham::Hamiltonian)\n\nDiagonalize the Hamiltonian at all momentum points in the BZ.\n\n\n\n\n\n","category":"method"},{"location":"Hamiltonian/#TightBindingToolkit.Hams.FillHamiltonian-Tuple{UnitCell, Vector{Float64}}","page":"Hamiltonian","title":"TightBindingToolkit.Hams.FillHamiltonian","text":"FillHamiltonian(uc::UnitCell, k::Vector{Float64} ; SpinMatrices::Vector{Matrix{ComplexF64}} = SpinMats((uc.localDim-1)//2)) :: Matrix{ComplexF64}\n\nReturns the Hamiltonian at momentum point k, corresponding to the bonds present in UnitCell.\n\n\n\n\n\n","category":"method"},{"location":"susceptibility/#TightBindingToolkit.susceptibility","page":"Magnetic susceptibility","title":"TightBindingToolkit.susceptibility","text":"","category":"section"},{"location":"susceptibility/","page":"Magnetic susceptibility","title":"Magnetic susceptibility","text":"Modules = [TightBindingToolkit, TightBindingToolkit.suscep]\nPrivate = false\nPages   = [\"Susceptibility.jl\"]","category":"page"},{"location":"susceptibility/#TightBindingToolkit.suscep.susceptibility","page":"Magnetic susceptibility","title":"TightBindingToolkit.suscep.susceptibility","text":"susceptibility is a data type representing the magnetic response, χ^ab(Q  Ω) for a general tight-binding Model.\n\nAttributes\n\nM       ::  Model: the given model.\nQs      ::  Vector{Vector{Float64}}: the set of momentum points over which χ^ab(Q  Ω) is calculated.\nOmegas  ::  Vector{Float64}: the set of energies over which χ^ab(Q  Ω) is calculated.\nSpread  ::  Float64 : the finite spread when summing over delta functions.\nchis    ::  Dict: a dictionary containing χ^ab(Q  Ω) for the different directions e.g. chis[\"xx\"] etc.\n\nInitialize this structure using \n\nsusceptibility(M::Model , Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, [], Omegas, eta, Dict())\nsusceptibility(M::Model , Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, Qs, Omegas, eta, Dict())\n\n\n\n\n\n","category":"type"},{"location":"susceptibility/#TightBindingToolkit.suscep.FillChis!-Tuple{susceptibility}","page":"Magnetic susceptibility","title":"TightBindingToolkit.suscep.FillChis!","text":"FillChis!(chi::susceptibility; fill_BZ::Bool=false, a::Int64=3, b::Int64=3)\n\nfunction to calculate susceptibility at a all given Ω=Omegas, but for all Q present in the given path, and along a fixed direction given by a and b.\n\n\n\n\n\n","category":"method"},{"location":"susceptibility/#TightBindingToolkit.suscep.FindChi-Tuple{Vector{Float64}, Float64, Model}","page":"Magnetic susceptibility","title":"TightBindingToolkit.suscep.FindChi","text":"FindChi(Q::Vector{Float64}, Omega::Float64 , M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2) --> ComplexF64\n\nfunction to calculate susceptibility at a fixed Ω=Omega, and Q, and along a fixed direction given by a and b.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingTookit.jl","page":"Introduction","title":"TightBindingTookit.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"TightBindingToolkit.jl is a Julia package meant for constructing, and obtaining useful properties of generic tight-binding models. It supports any lattice structure, with any user-defined bonds on that lattice. It also has support for any spin of the particle hopping on the lattice.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Currently supported :","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Custom Unit Cell Construction\nCorresponding Brillouin Zone Construction\nHamiltonian, given a Unit Cell and a Brillouin Zone\nDiagonalizing the Hamiltonian in momentum space to get band structures and wavefunctions\nDensity of State \nFilling the model at given chemical potential, and calculating gaps\nGetting correlation functions\nGetting Berry curvature and Chern numbers\nGetting magnetic susceptibility in any direction, at any momentum, and energy","category":"page"},{"location":"UnitCell/#TightBindingToolkit.UnitCell","page":"Unit Cell","title":"TightBindingToolkit.UnitCell","text":"","category":"section"},{"location":"UnitCell/","page":"Unit Cell","title":"Unit Cell","text":"Modules = [TightBindingToolkit, TightBindingToolkit.UCell]\nPrivate = false\nPages   = [\"UnitCell.jl\"]","category":"page"},{"location":"UnitCell/#TightBindingToolkit.UCell.Bond","page":"Unit Cell","title":"TightBindingToolkit.UCell.Bond","text":"Bond{T<:Number} is a data type representing a general bond on a lattice.\n\nAttributes\n\nbase::Int64: sub-lattice of the initial site on the bond.\nbase::Int64: sub-lattice of the final site on the bond.\noffset::Vector{Int64}: the difference of the unit cells in which these sublattices belong to, in units of the lattice basis vectors.\nmat::Matrix{T}: the matrtix describing this bond –> can be hopping for partons, or spin-exchange for spins.\ndist::Float64: the distance b/w the two sites = length of bond.\nlabel::String: some string label to mark the bond type.\n\n\n\n\n\n","category":"type"},{"location":"UnitCell/#TightBindingToolkit.UCell.UnitCell","page":"Unit Cell","title":"TightBindingToolkit.UCell.UnitCell","text":"UnitCell is a data type representing a general unit cell of a lattice.\n\nAttributes\n\nprimitives  :: Vector{ Vector{ Float64 } }: primitive bases vectors of the lattice.\nbasis       :: Vector{ Vector{ Float64 } }: positions of the sub-lattices.\nbonds       :: Vector{Bond}: the set of all bonds defining a lattice.\nfields      :: Vector{ Vector{Float64}} : the fields oneach basis site.\nlocalDim    :: Int64: Local Hilbert space dimension ( e.g. 3 for classical spins, 2 for spin-1/2 electrons ).\nBC          :: Vector{ ComplexF64 }: boundary conditions, in the form of [e^{ιθ_i}]. e.g. θ=0 for PBC, θ=π for APBC, and so on.\n\nInitialize this structure using \n\nUnitCell( as::Vector{Vector{Float64}} , localDim::Int64)\n\n\n\n\n\n","category":"type"},{"location":"UnitCell/#TightBindingToolkit.UCell.ModifyBonds!-Tuple{UnitCell, Float64, Matrix{var\"#s10\"} where var\"#s10\"<:Number}","page":"Unit Cell","title":"TightBindingToolkit.UCell.ModifyBonds!","text":"ModifyBonds!(uc::UnitCell, dist::Float64, newMat::Matrix{<:Number})\nModifyBonds!(uc::UnitCell, label::String, newMat::Matrix{<:Number})\n\nModify an existing bond in the UnitCell with the given label, or at a given distance=dist, to the given bond matrix.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.ModifyFields!-Tuple{UnitCell, Int64, Vector{Float64}}","page":"Unit Cell","title":"TightBindingToolkit.UCell.ModifyFields!","text":"ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})\nModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})\n\nModify the on-site fields in the UnitCell, either one at a time, or all of them.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.ModifyIsotropicFields!-Tuple{UnitCell, Vector{Float64}}","page":"Unit Cell","title":"TightBindingToolkit.UCell.ModifyIsotropicFields!","text":"ModifyIsotropicFields!(uc::UnitCell, newField::Vector{Float64})\nModifyIsotropicFields!(uc::UnitCell, newField::Float64, dim::Int64)\n\nModify the on site field uniformly, on all sublattices. The optional argument dim is if you want to only modify one of the 4 elements of on-site fields (3 Zeeman and 1 chemical potential).\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.RemoveBonds!-Tuple{UnitCell, String}","page":"Unit Cell","title":"TightBindingToolkit.UCell.RemoveBonds!","text":"RemoveBonds!(uc::UnitCell, dist::Float64)\nScaleBonds!(uc::UnitCell, label::String)\n\nRemove an existing bond in the UnitCell with the given label, or at a given distance=dist.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.ScaleBonds!-Tuple{UnitCell, Float64, Number}","page":"Unit Cell","title":"TightBindingToolkit.UCell.ScaleBonds!","text":"ScaleBonds!(uc::UnitCell, dist::Float64, scale::Number)\nScaleBonds!(uc::UnitCell, label::String, scale::Number)\n\nScale the matrix of an existing bond in the UnitCell with the given label, or at a given distance=dist, by the given scaling factor.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.addAnisotropicBond!-Tuple{UnitCell, Int64, Int64, Vector{Int64}, Number, Float64, String}","page":"Unit Cell","title":"TightBindingToolkit.UCell.addAnisotropicBond!","text":"addAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )\naddAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )\n\nAdd a bond with the given attributes to UnitCell. If given mat attribute is a number, it is converted into a 1x1 matrix when entered into the bond.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.addBasisSite!-Tuple{UnitCell, Vector{Float64}}","page":"Unit Cell","title":"TightBindingToolkit.UCell.addBasisSite!","text":"addBasisSite!( uc::UnitCell , position::Vector{Float64} )\naddBasisSite!( uc::UnitCell , position::Vector{Float64} , field::Vector{Float64} )\n\nAdd a sublattice to the UnitCell  at the given real-space position, with an on-site field.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.addIsotropicBonds!-Tuple{UnitCell, Float64, Number, String}","page":"Unit Cell","title":"TightBindingToolkit.UCell.addIsotropicBonds!","text":"addIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))\naddIsotropicBonds!( uc::UnitCell , dist::Float64 , mats::Matrix{<:Number} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )\n\nAdd a set of \"isotropic\" bonds, which are the same for each pair of sites at the given distance.  If given mat attribute is a number, it is converted into a 1x1 matrix when entered into the bond. The input checkOffsetRange must be adjusted depending on the input distance.  The optional input subs is meant for isotropic bonds when only a subset of sublattices are involved.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.getDistance-Tuple{UnitCell, Int64, Int64, Vector{Int64}}","page":"Unit Cell","title":"TightBindingToolkit.UCell.getDistance","text":"getDistance(uc::UnitCell, base::Int64, target::Int64, offset::Vector{Int64}) --> Float64\n\nget the distance between site at position (0, base) and (R, target), where R = offset, when written in units of the unit cell primitive vectors.\n\n\n\n\n\n","category":"method"},{"location":"UnitCell/#TightBindingToolkit.UCell.isSameBond-Tuple{Bond, Bond}","page":"Unit Cell","title":"TightBindingToolkit.UCell.isSameBond","text":"Function to check if two bond objects are describing the same physical bond, just inverted! \n\n\n\n\n\n","category":"method"},{"location":"Chern/#TightBindingToolkit.Chern","page":"Chern Numbers","title":"TightBindingToolkit.Chern","text":"","category":"section"},{"location":"Chern/","page":"Chern Numbers","title":"Chern Numbers","text":"Modules = [TightBindingToolkit, TightBindingToolkit.Chern]\nPrivate = false\nPages   = [\"Chern.jl\"]","category":"page"},{"location":"Chern/#TightBindingToolkit.Chern.ChernNumber-Tuple{Hamiltonian, Vector{Int64}}","page":"Chern Numbers","title":"TightBindingToolkit.Chern.ChernNumber","text":"ChernNumber(Ham::Hamiltonian, subset::Vector{Int64}) --> Float64\n\nFunction to get Chern numbers given a Hamiltonian and a subset of bands\n\n\n\n\n\n","category":"method"},{"location":"Chern/#TightBindingToolkit.Chern.FieldStrength-Tuple{Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}","page":"Chern Numbers","title":"TightBindingToolkit.Chern.FieldStrength","text":"FieldStrength(Links::Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}) --> Matrix{ComplexF64}\n\nFunction to calculate the product of the links over each plaquette on the BZ grid. This is the generalized Bery curvature for multiple degenerate bands.\n\n\n\n\n\n","category":"method"},{"location":"Chern/#TightBindingToolkit.Chern.FindLinks-Tuple{Hamiltonian, Vector{Int64}}","page":"Chern Numbers","title":"TightBindingToolkit.Chern.FindLinks","text":"FindLinks(Ham::Hamiltonian, subset::Vector{Int64}) --> Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}\n\nFunction to get the linking matrices on each neighbouring point in the BZ. On a bond connecting ki and kj, the linking matrix U is defined such that U[m, n] = <v^m[ki]|v^n[kj]> where states[kj[1], kj[2], :, m]  v^m[kj], the mth eigenstate at momentum kj.\n\n\n\n\n\n","category":"method"}]
}
