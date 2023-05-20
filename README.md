# TightBindingToolkit

[![Build Status](https://github.com/sreekar-voleti/TightBindingToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sreekar-voleti/TightBindingToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

TightBindingToolkit.jl is a Julia package meant for constructing, and obtaining useful properties of generic tight-binding models. It supports any lattice structure, with any user-defined bonds on that lattice. It also has support for any spin of the particle hopping on the lattice.

Currently supported :
* Custom Unit Cell Construction
* Corresponding Brillouin Zone Construction
* Hamiltonian, given a Unit Cell and a Brillouin Zone
* Diagonalizing the Hamiltonian in momentum space to get band structures and wavefunctions
* Density of State 
* Filling the model at given chemical potential, and calculating gaps
* Getting correlation functions
* Getting Berry curvature and Chern numbers
* Getting magnetic susceptibility in any direction, at any momentum, and energy
