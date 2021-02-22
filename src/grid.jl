mutable struct AdaptiveGrid{dim,T <: Real} <: JuAFEM.AbstractGrid{dim}
    rootcells::Vector{Element}

    cells::Vector{Element} # leafcells - don't really need it, we can just collect all leaves that are masters among the cells
    faces::Vector{Element} # leaffaces - don't really need it, we can just collect all leaves that are masters among the cells
    nodes::Vector{Node{dim,T}} # actually don't need that either

    cellsets::Dict{String,Set{Element}}
    nodesets::Dict{String,Set{Node{dim,T}}}
    facesets::Dict{String,Set{Element}} # careful about master-slave-subtleties
end

function update_leafmesh!(grid::AdaptiveGrid)
    leafcells = collect_children(grid.rootcells)
    leaffacecandidates = collect(Iterators.flatten((x -> x.boundaries).(leafcells)))
    leaffaces = filter(x -> x.ismaster, leaffacecandidates)
    nodecandidates = collect(Iterators.flatten(vertices.(leafcells)))
    nodes = unique!(x -> x, nodecandidates)
    grid.cells = leafcells
    grid.faces = leaffaces
    grid.nodes = nodes
    return grid
end

function JuAFEM.getcells(grid::AdaptiveGrid)
    return (x -> eltocell(x, grid.nodes)).(grid.cells)
end

function JuAFEM.getcells(grid::AdaptiveGrid, v::Union{Int,Vector{Int}})
    if typeof(v) <: Integer return eltocell(grid.cells[v], grid.nodes) end
    return (x -> eltocell(x, grid.nodes)).(grid.cells[v])
end

function JuAFEM.getcells(grid::AdaptiveGrid, set::String)
    leaves = collect_children(grid.cellsets[set])
    return (x -> eltocell(x, grid.nodes)).(leaves)
end

function JuAFEM.getnodes(grid::AdaptiveGrid, set::String)
    return collect(grid.nodesets[set])
end

function JuAFEM.getnodeset(grid::AdaptiveGrid, set::String)
    return Set((x -> findfirst(a -> a == x, grid.nodes)).(grid.nodesets[set]))
end

function JuAFEM.getnodesets(grid::AdaptiveGrid)
    n = Dict{String,Set{Int}}()
    for d in grid.nodesets
        n[d[1]] = getnodeset(grid, d[1])
    end
    return n
end

function JuAFEM.n_faces_per_cell(grid::AdaptiveGrid)
    return nfaces(grid.cells[1])
end

function JuAFEM.getcellset(grid::AdaptiveGrid, set::String)
    s = Set(findall(x -> ultimate_ancestor(x) ∈ grid.cellsets[set], grid.cells))
    return s
end

function JuAFEM.getcellsets(grid::AdaptiveGrid)
    c = Dict{String,Set{Int}}()
    for d in grid.cellsets
        c[d[1]] = getcellset(grid, d[1])
    end
    return c
end

function getfaceidx(grid::AdaptiveGrid, f::Element)
    i = findfirst(x -> x == f.boundaryof, grid.cells)
    j = findfirst(x -> x == f, f.boundaryof.boundaries)
    return (i, j)
end

function JuAFEM.getfaceset(grid::AdaptiveGrid, set::String) # not clear if this is correct
    leaffaces = collect_children(grid.facesets[set])
    faces = Set(get_master.(leaffaces))
    S = (x -> getfaceidx(grid, x)).(faces)
    return Set(S)
end

function JuAFEM.getfacesets(grid::AdaptiveGrid)
    f = Dict{String,Set{Tuple{Int,Int}}}()
    for d in grid.facesets
        f[d[1]] = getfaceset(grid, d[1])
    end
    return f
end

# use the ones in grid.jl of JuAFEM?
_check_setname(dict, name) = haskey(dict, name) && throw(ArgumentError("there already exists a set with the name: $name"))
_warn_emptyset(set) = length(set) == 0 && @warn("no entities added to set")

function JuAFEM.addcellset!(grid::AdaptiveGrid, name::String, cellid::Union{Set{Int},Vector{Int}})
    _check_setname(grid.cellsets,  name)
    cells = Set(cellid)
    _warn_emptyset(cells)
    grid.cellsets[name] = Set((x -> grid.cells[x]).(cells))
    grid
end

# this is probably very suboptimal for our purposes
function JuAFEM.addcellset!(grid::AdaptiveGrid, name::String, f::Function; all::Bool=true)
    _check_setname(grid.cellsets,  name)
    cells = Set{Int}()
    for (i, cell) in enumerate(getcells(grid))
        pass = all
        for node_idx in cell.nodes
            node = grid.nodes[node_idx]
            v = f(node.x)
            all ? (!v && (pass = false; break)) : (v && (pass = true; break))
        end
        pass && push!(cells, i)
    end
    _warn_emptyset(cells)
    grid.cellsets[name] = Set(grid.cells[cells])
    grid
end

function JuAFEM.addfaceset!(grid::AdaptiveGrid, name::String, faceid::Set{Tuple{Int,Int}})
    _check_setname(grid.facesets, name)
    faceset = Set(faceid)
    _warn_emptyset(faceset)
    grid.facesets[name] = (x -> (grid.cells[x[1]]).boundaries[x[2]]).(faceset)
    grid
end

function JuAFEM.addfaceset!(grid::AdaptiveGrid, name::String, f::Function; all::Bool=true)# does this work?
    _check_setname(grid.facesets, name)
    faceset = Set{Tuple{Int,Int}}()
    for (cell_idx, cell) in enumerate(getcells(grid))
        for (face_idx, face) in enumerate(faces(cell))
            pass = all
            for node in face
                v = f(node.x)
                all ? (!v && (pass = false; break)) : (v && (pass = true; break))
            end
            pass && push!(faceset, (cell_idx, face_idx))
        end
    end
    _warn_emptyset(faceset)
    grid.facesets[name] = (x -> (grid.cells[x[1]]).boundaries[x[2]]).(faceset)
    grid
end

function JuAFEM.addnodeset!(grid::AdaptiveGrid, name::String, nodeid::Union{Vector{Int},Set{Int}})
    _check_setname(grid.nodesets, name)
    grid.nodesets[name] = Set((x -> grid.nodes[x]).(Set(nodeid)))
    _warn_emptyset(grid.nodesets[name])
    grid
end

function JuAFEM.addnodeset!(grid::AdaptiveGrid, name::String, f::Function)# this can be simplified
    _check_setname(grid.nodesets, name)
    nodes = Set{Int}()
    for (i, n) in enumerate(getnodes(grid))
        f(n.x) && push!(nodes, i)
    end
    grid.nodesets[name] = Set((x -> grid.nodes[x]).(Set(nodeid)))
    _warn_emptyset(grid.nodesets[name])
    grid
end

getcelltype(grid::AdaptiveGrid) = typeof(grid.cells[1])


# these two should be obsolete

function JuAFEM.getcoordinates!(x::Vector{Vec{dim,T}}, grid::AdaptiveGrid, cell::Int) where {dim,T}
    @inbounds for i in 1:length(x)
        x[i] = grid.cells[cell].nodes[i].x
    end
end

function JuAFEM.getcoordinates(grid::AdaptiveGrid, cell::Int)
    nodeex = grid.cells[cell].nodes[1]
    dim = typeof(nodeex).parameters[1]
    T = typeof(nodeex).parameters[2]
    return (n -> n.x).(grid.cells[cell].nodes)# ::Vector{Vec{dim,T}}
end


function JuAFEM.faces(e::Element, grid::AdaptiveGrid)
    getindex = (n -> findfirst(x -> x == n, grid.nodes))
    res = (f -> getindex.(f)).(faces(e))
    return res
end


function JuAFEM.edges(e::Element, grid::AdaptiveGrid)
    getindex = (n -> findfirst(x -> x == n, grid.nodes))
    res = (b -> getindex.(b)).(edges(e))
    return res
end

function JuAFEM.vertices(e::Element, grid::AdaptiveGrid)
    getindex = (n -> findfirst(x -> x == n, grid.nodes))
    res = getindex.(vertices(e))
    return res
end


function refine!(grid::AdaptiveGrid, cells::Vector{Int},  maxdepth::Int=-1)
    cs = grid.cells[cells]
    for c in cs
        if maxdepth != -1 && getdepth(c) > maxdepth error("too deep, can not split") end
        isotropic_split!(c)
        @debug check(grid)
    end
    grid = update_leafmesh!(grid)
    return grid
end

function refine!(grid::AdaptiveGrid, cells::Vector{Int}, axes::Vector{Vector{Int}},  maxdepth::Int=-1)
    cs = grid.cells[cells]
    for i in 1:length(cs)
        if maxdepth != -1 && getdepth(cs[i]) > maxdepth error("too deep, can not split") end
        split!(cs[i], axes[i])
        @debug check(grid)
    end
    grid = update_leafmesh!(grid)
    return grid
end

function derefine!(grid::AdaptiveGrid, cells::Vector{Int})
    cs = grid.cells[cells]
    par = unique!((x -> x.parent).(cs))
    if nothing ∈ par error("trying to derefine a root leaf, shame on you") end
    for c in par
        if issubset(Set(c.children), Set(cs))
            derefine!(c)
            @debug check(grid)
        else
            error("not all children of a cell marked for derefinement")
        end
    end
    grid = update_leafmesh!(grid)
    return grid
end


function getncinfo(grid::AdaptiveGrid{2})
    output = []
    for f in grid.faces
        ns, s = getinteriornodeinfo(f)
        for i in 1:length(ns)
            nodenr = findfirst(x -> x == ns[i], grid.nodes)
            masterfaceidx = getfaceidx(grid, f)
            slavefaceidxes = (x -> getfaceidx(grid, x)).(s[i])
            push!(output, [nodenr, masterfaceidx, slavefaceidxes])
        end
    end
    return output
end

function getedgeindex(grid::AdaptiveGrid, e::Element)
    i, j = getfaceidx(grid, e.boundaryof)
    k = findfirst(x -> x == e, e.boundaryof.boundaries)
    return (i, j, k)
end

getidx(grid::AdaptiveGrid{3}, e::Element{2,4,4}) = getfaceidx(grid, e)
getidx(grid::AdaptiveGrid{3}, e::Element{1,2,2}) = getedgeindex(grid, e)

function getncinfo(grid::AdaptiveGrid{3})
    output = []
    for f in grid.faces
        for n in getnodeneighborhood(f.boundaryof)
            if (node_in_el(f, n)) && !(n ∈ f.nodes)
                nodenr = findfirst(x -> x == n, grid.nodes)
                minfo = containmentinfo(f, n)
                relevantelements = (x -> containmentinfo(x, n)).(union(collect_children(f.attachedto), [f]))
                relevantelements = filter!(x -> length(x) != 0, relevantelements)
                relevantelements = collect(Iterators.flatten(relevantelements))
                slinfo = filter(x -> n ∈ vertices(x), relevantelements)
                minfo = filter(x -> !(n ∈ vertices(x)), relevantelements)
                push!(output, [nodenr, (x -> getidx(grid, x)).(minfo), (x -> getidx(grid, x)).(slinfo) ])
            end
        end
    end
    return output
end


# for testing purposes
@debug function check(grid::AdaptiveGrid)
    leafcells = collect(Iterators.flatten(collect_children.(grid.rootcells)))
    leaffacecandidates = collect(Iterators.flatten((x -> x.boundaries).(leafcells)))
    for l in leaffacecandidates
        get_master(l)
    end
end
