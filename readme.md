 
This code is the result of a project supervised by Dennis Ogiermann and Maximilian Köhler and carried out by Esther Koch and Leonid Ryvkin. 


This projects implements an adapive grid subtype of the abstract type JuAFEM.AbstractGrid, making it possible to carry out certain Adapive Finite Element computations in the JuAFEM framework.

Currently linear quadrilateral elements in 2D and linear hexahedral elements in 3D are supported, with the possibility to refine elements isotropically into quadrants (octants) or anisotropically (halves and quarters).

# Installation

Install the JuAFEM version with support for more general grids from the Pkg REPL in julia:
```
pkg> add https://github.com/koehlerson/JuAFEM.jl.git#fe/dofhandler
```

Then install FEMRefinements.jl

```
pkg> add https://github.com/Ljon4ik4/FEMRefinements.jl.git
```
