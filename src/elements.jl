mutable struct Element{dim,N,M} <: JuAFEM.AbstractCell{dim,N,M}
    nodes::Vector{Node}
    boundaries::Vector{Element}

    parent::Union{Nothing,Element}
    children::Vector{Element}

    boundaryof::Union{Nothing,Element}
    attachedto::Vector{Element}
    ismaster::Bool
    function Element(n, p)
        sideorders = [ [] ,[[1], [2]],  [[1,2],[2,3],[3,4],[4,1]],
            [[1,4,3,2], [1,2,6,5],[2,3,7,6], [3,4,8,7], [1,5,8,4], [5,6,7,8]]]
        N = length(n)
        dim = trunc(Int, log2(N))
        M = 2 * dim
        bds = (x -> Element(n[x], nothing)).(sideorders[dim + 1])
        v = Vector{Element}(bds)
        self = new{dim,N,M}(n, v, p, Vector{Element}(), nothing, Vector{Element}(), true)
        for x in self.boundaries
            x.boundaryof = self
        end
        return  self
    end
end

function ultimate_ancestor(e::Element)
    if (e.parent === nothing) return e end
    return ultimate_ancestor(e.parent)
end

function collect_children(e::Element)
    if length(e.children) == 0 return [e] end
    return collect(Iterators.flatten(collect_children.(e.children)))
end

function collect_children(v::Union{AbstractSet{C},AbstractVector{C}}) where C <: Element
    return collect(Iterators.flatten(collect_children.(v)))
end

# 0-dim things
function JuAFEM.vertices(e::Element)
    return tuple(e.nodes...)
end

# 1-codim things
function JuAFEM.faces(e::Element)
    return tuple((vertices.(e.boundaries))...)
end

# meant are the 2-codimensional ones, is there an order convention here?
function JuAFEM.edges(e::Element)
    edg = collect(Iterators.flatten((x -> x.boundaries).(e.boundaries)))
    edg = vertices.(edg)
    # all edges exist (in 3d) twice
    return unique!(x -> Set(x), edg) # exclamation mark or better without
end


JuAFEM.default_interpolation(::Type{Element{2,4,4}})  = Lagrange{2,RefCube,1}()
JuAFEM.default_interpolation(::Type{Element{3,8,6}})  = Lagrange{3,RefCube,1}()


function Base.show(io::IO, e::Element)
    print(io, (n -> n.x).(vertices(e)))
end

function eltocell(e::Element{dim,N,M}, n::Vector{Node{dim,T}}) where {dim,N,M,T} # formerly quadtocell
    nodenumbers = (x -> findfirst(y -> y == x, n)).(e.nodes)
    return JuAFEM.Cell{dim,N,M}(tuple(nodenumbers...))
end

# only works correctly on leaves
function get_master(e::Element) # slaves point to master, master points to parent[s!!!] of a list of slaves if several exists,
    if e.ismaster == true return e end
    i = findfirst(x -> x.ismaster == true, e.attachedto)
    return e.attachedto[i]
end

getdepth(n::Nothing) = -1
getdepth(e::Element) =  1 + getdepth(e.parent)



# makes the first one a master, does not consider children etc.
function attach_trivially!(e::Element, f::Element)
    masterface = nothing
    for i in e.boundaries
        for j in f.boundaries
            if issubset(vertices(i), vertices(j)) # better issetequal?
                i.attachedto = [j]
                j.attachedto = [i]
                i.ismaster = true
                j.ismaster = false
                masterface = i
            end
        end
    end
    return masterface
end


function getcontainednodes(f::Element{1,2,2})
    if !f.ismaster error("only applicable to master elements") end
    leaves = collect_children(union(f.attachedto, [f]))
    ns = unique!(collect(Iterators.flatten(vertices.(leaves))))
    return ns
end

function getinteriornodes(f::Element{1,2,2})
    ns = getcontainednodes(f)
    ns = deleteat!(ns, findall(x -> x ∈ vertices(f), ns))
    return ns
end

function getinteriornodeinfo(f::Element{1,2,2})
    ns = getinteriornodes(f)
    sfaces = collect_children(f.attachedto)
    sfaceinfo = (n -> findall(e -> n ∈ vertices(e), sfaces)).(ns)
    return ns, (x -> sfaces[x]).(sfaceinfo)
end


function getneighbourhood(v::Vector{E}) where E <: Element
    output = copy(v)
    for e in v
        for f in e.boundaries
            neigh = (x -> x.boundaryof).(collect_children(union(get_master(f).attachedto, [get_master(f)])))
            output = vcat(output, neigh)
        end
    end
    return convert(Vector{E}, unique!(output))
end

function getnodeneighborhood(e::Element)
    cells = getneighbourhood(getneighbourhood([e])) # sufficient for up to 3D
    return unique!(collect(Iterators.flatten(vertices.(cells))))
end


# only works in dim 3?
function node_in_el(e::Element{1,2,2}, n::Node)
    v1 = e.nodes[1].x - n.x
    v2 = e.nodes[2].x - n.x
    if length(v1) == 3 && !isapprox(norm(cross(v1, v2)), 0) return false end
    s = dot(v1, v2)
    if s > 0 return false end # isapprox version of that?
    return true
end

function getvec(e::Element{1,2,2})
    return e.nodes[2].x - e.nodes[1].x
end
# only works in dim 3?
function node_in_el(e::Element{2,4,4}, n::Node)
    s = getvec.(e.boundaries)
    v = (y -> y.x - n.x).(e.nodes)
    if !isapprox(dot(cross(s[1], s[2]), v[1]), 0) return false end
    rhs = (i -> dot(s[i], v[i]) / dot(s[i], s[i])).(collect(1:4))
    if (maximum(rhs) <= 1) && (minimum(rhs) >= 0) return true end
    if (maximum(rhs) <= 0) && (minimum(rhs) >= -1) return true end
    return false
end

function containmentinfo(e::Element{2,4,4}, n::Node)
    if !node_in_el(e, n) return [] end
    bds = findall(b -> node_in_el(b, n), e.boundaries)
    if length(bds) == 0 return [e] end
    return e.boundaries[bds]
end
