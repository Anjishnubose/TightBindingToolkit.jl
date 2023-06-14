using Documenter
using TightBindingToolkit

makedocs(
    build       =   "build" ,
    sitename    =   "TightBindingToolkit.jl"    ,
    modules     =   [TightBindingToolkit, TightBindingToolkit.UCell, TightBindingToolkit.Parameters, TightBindingToolkit.BZone, TightBindingToolkit.Hams, 
                            TightBindingToolkit.TBModel, TightBindingToolkit.Chern, TightBindingToolkit.suscep]   ,
    pages = [
        "Introduction"              =>  "index.md",
        "Unit Cell"                 =>  "UnitCell.md",
        "Parameters"                =>  "Params.md",
        "Brillouin Zone"            =>  "BZ.md",
        "Hamiltonian"               =>  "Hamiltonian.md",
        "Tight Binding Model"       =>  "Model.md",
        "Chern Numbers"             =>  "Chern.md",
        "Magnetic susceptibility"   =>  "susceptibility.md"  
    ]
)

deploydocs(
    repo = "github.com/Anjishnubose/TightBindingToolkit.jl.git",
    devbranch = "main"
)