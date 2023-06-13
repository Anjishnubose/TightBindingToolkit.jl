using Documenter
using TightBindingToolkit

makedocs(
    build       =   "build" ,
    sitename    =   "TightBindingToolkit.jl"    ,
    modules     =   Module[TightBindingToolkit]   ,
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