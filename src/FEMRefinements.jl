module FEMRefinements

using JuAFEM
# using JuAFEM:AbstractCell, Cell

using LinearAlgebra
using Logging

import JuAFEM:vertices, edges, faces

export addcellset!, addfaceset!, addnodeset!,
compute_vertex_values, default_interpolation, edges, faces, getcells, getcellset, getcellsets,
getcelltype, getcoordinates!, getdim, getfaceset, getfacesets,
getnnodes, getnodeset, getnodesets, n_faces_per_cell, nnodes_per_cell,
vertices, generate_grid, getncells, getnodes, getcoordinates

export generate_grid2D, generate_grid3D, refine!, derefine!, Element, AdaptiveGrid,
getncinfo, getedgeindex, getfaceidx, getmuchmorelogicalncinfo

export facedofs, find_local_coordinate

# a few of the utility methods from JuAFEM are ugly [if not wrong]. Also filter instead of findall anod some more advanced collections operations at some places might be useful
include("elements.jl")
include("grid.jl")
include("grid_generators.jl")
include("refinements.jl")
include("utils.jl")

using ColorSchemes # do we need that?
import Makie, AbstractPlotting # do we need that?

include("visualization.jl")


# the current implementation of JuAFEM iterators does not allow for grids with variable eltype
import Base.eltype
function Base.eltype(v::Vector{Element})
    if length(v) == 0 return Element end
    return typeof(v[1])
end
export eltype

end
