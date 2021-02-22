function getcenter(e::Element)
    dim = typeof(e.nodes[1]).parameters[1]
    T = typeof(e.nodes[1]).parameters[2]
    S = (n -> n.x).(e.nodes)
    v = sum(S) / length(S)
    return Node{dim,T}(v)
end

function substitute_nodes!(newnodes::Array{N}, nodes::Vector{N}) where N <: Node
    for i in 1:(length(newnodes))
        exnodes = findall(n -> isapprox(n.x, newnodes[i].x), nodes)
        if (length(exnodes)) == 1 newnodes[i] = nodes[exnodes[1]] end
        @debug if (length(exnodes)) > 1 error("node in nodelist several times") end
    end
    return newnodes
end

function findlowestancestorwithnodes(e::Element, n::Vector) # nodevector
    if issubset(n, e.nodes) return e end
    if e.parent === nothing return nothing end
    return findlowestancestorwithnodes(e.parent, n)
end

function have_identical_leaves(e::Element, f::Element)
    le = collect_children(e)
    fe = collect_children(f)
    return issetequal(le, fe)
end

function make_master!(e::Element)
    if length(e.attachedto) == 0 return e end
    @debug if (length(e.children) != 0) error("trying to make a non-leaf into a master") end
    oldmaster = get_master(e)
    if oldmaster == e return e end
    @debug if length(oldmaster.children) != 0 error("oldmaster has children") end
    @debug if !issubset(vertices(oldmaster), vertices(e)) error("nodes of old and new masters must coincide") end

    e.attachedto = copy(oldmaster.attachedto) # better with union and filter
    push!(e.attachedto, oldmaster)
    deleteat!(e.attachedto, findall(x -> have_identical_leaves(x, e), e.attachedto))
    e.attachedto = unique!(e.attachedto)
    e.ismaster = true
    slaveleaves = collect_children(e.attachedto)
    (x -> x.attachedto = [e]).(slaveleaves) # collate this?
    (x -> x.ismaster = false).(slaveleaves) # should suffice on oldmaster
    @debug if !e.ismaster error("somehow e seems to be one of its own slaveleaves") end
    return e
end

function descend_until_split(e::Element)
    if length(e.children) != 1 return e end
    return descend_until_split(e.children[1])
end

function attach_children!(newels::Vector{Element{dim,N,M}}, e::Element) where {dim,N,M} # used on face-level
    @debug if length(e.children) != 0 error("already has children") end
    @debug if length(newels) == 0 error("nothing to attach") end

    if (e.ismaster == false)
        (x -> x.attachedto = [get_master(e)]).(newels)
        (x -> x.ismaster = false).(newels)
        e.children = newels
        (x -> x.parent = e).(newels)
        return [get_master(newels[1])]
    end

    altmasters = findall(x -> length(x.children) == 0, descend_until_split.(e.attachedto))
    if length(altmasters) != 0
        newmaster = descend_until_split(e.attachedto[altmasters[1]])
        make_master!(newmaster)
        e.children = newels
        (x -> x.parent = e).(newels)
        (x -> x.attachedto = [newmaster]).(newels)
        (x -> x.ismaster = false).(newels)
        return [get_master(newels[1])]
    end

    e.children = newels
    (x -> x.parent = e).(newels)
    (x -> x.ismaster = true).(newels)

    for m in newels
        for oldslaves in e.attachedto
            for s in descend_until_split(oldslaves).children
                if vertices(s) ⊆ vertices(m)
                    (x -> x.attachedto = [m]).(collect_children(s))
                    push!(m.attachedto, s)
                end
            end
        end
    end

    newelattachments = collect(Iterators.flatten((x -> x.attachedto).(newels)))
    if !(collect_children(e.attachedto) ⊆ collect_children(newelattachments))
        error("not all of the attachments could be propagated the children of e")
    end

    return get_master.(newels)
end

function unattach_children!(e::Element)
    if length(e.children) == 0 return e end
    che = collect_children(e)
    toobigattachments = findall(x -> nothing == findlowestancestorwithnodes(x, e.nodes), get_master.(che))
    if length(toobigattachments) != 0
        e.children = []
        e.ismaster = false
        e.attachedto = [get_master(che[toobigattachments[1]])]
        return e
    end

    attach = Vector{Element}()
    for c in che
        if c.ismaster
            attach = vcat(attach, (c.attachedto))
        else
            attach = vcat(attach, collect_children(get_master(c).attachedto))
            attach = push!(attach, get_master(c))
        end
    end
    attach = (x -> findlowestancestorwithnodes(x, e.nodes)).(attach)
    attach = filter!(x -> length(intersect(collect_children(e), collect_children(x))) == 0, attach)
    e.children = []
    e.ismaster = true
    e.attachedto = unique!(convert(Vector{Element}, attach))
    (x -> x.ismaster = false).(collect_children(e.attachedto)) # elegantly collate this?
    (x -> x.attachedto = [e]).(collect_children(e.attachedto))
    get_master(e)
    return e
end

function derefine!(e::Element)
    @debug if length(e.children) == 0 error("can not derefine leaf") end
    for x in e.children
        if length(x.children) != 0 error("can not derefine cell whose children are not leaves") end
    end
    cs = e.children
    e.children = []
    for b in e.boundaries
        unattach_children!(b)
    end
    return e
end


function isotropic_split!(e::Element{2,4,4})
    return split!(e, [1,2])
end

function split!(e::Element{2,4,4}, axes::Vector{Int})
    if length(e.children) != 0 error("can not split non-leaf") end
    xstep = 1
    ystep = 1
    if 1 ∈ axes xstep = 2 end
    if 2 ∈ axes ystep = 2 end
    if xstep * ystep == 1 error("come on, don't split something into 1 piece, it's embarassing") end
    dim = typeof(e.nodes[1]).parameters[1]
    T = typeof(e.nodes[1]).parameters[2]
    nodemat = Matrix{Node{dim,T}}(undef, 3, 3)
    nodemat[1,1] = e.nodes[1]
    nodemat[3,1] = e.nodes[2]
    nodemat[3,3] = e.nodes[3]
    nodemat[1,3] = e.nodes[4]
    nodemat[2,1] = getcenter(e.boundaries[1])
    nodemat[3,2] = getcenter(e.boundaries[2])
    nodemat[2,3] = getcenter(e.boundaries[3])
    nodemat[1,2] = getcenter(e.boundaries[4])
    nodemat[2,2] = getcenter(e)
    nodesnearby = unique!(collect(Iterators.flatten(getcontainednodes.(get_master.(e.boundaries)))))
    nodemat = substitute_nodes!(nodemat, nodesnearby)
    cellmat = Matrix{Element{2,4,4}}(undef, (xstep, ystep))
    for i in 1:xstep
        for j in 1:ystep
            cellmat[i,j] = Element([nodemat[i,j],nodemat[i + (2 ÷ xstep),j],nodemat[i + (2 ÷ xstep),j + (2 ÷ ystep)],nodemat[i,j + (2 ÷ ystep)]], e)
        end
    end

    for i in 1:xstep
        for j in 1:ystep
            if i < xstep attach_trivially!(cellmat[i,j], cellmat[i + 1,j]) end
            if j < ystep attach_trivially!(cellmat[i,j], cellmat[i,j + 1]) end
        end
    end

    allfaces = collect(Iterators.flatten((x -> x.boundaries).(cellmat)))
    for i in 1:length(e.boundaries)
        side = filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[i] == x), allfaces)
        attach_children!(convert(Array{Element{1,2,2}}, side), e.boundaries[i])
    end

    e.children = reshape(cellmat, length(cellmat)) # [cellmat[1,1],cellmat[2,1],cellmat[2,2],cellmat[1,2]]
    return e
end


function isotropic_split!(e::Element{3,8,6})
    return split!(e, [1,2,3])
end

function split!(e::Element{3,8,6}, axes::Vector{Int})
    if length(e.children) != 0 error("can not split non-leaf") end
    xstep = 1
    ystep = 1
    zstep = 1
    if 1 ∈ axes xstep = 2 end
    if 2 ∈ axes ystep = 2 end
    if 3 ∈ axes zstep = 2 end
    if xstep * ystep * zstep == 1 error("come on, don't split something into 1 piece, it's embarassing") end
    dim = typeof(e.nodes[1]).parameters[1]
    T = typeof(e.nodes[1]).parameters[2]
    nodemat = Array{Node{dim,T}}(undef, (3, 3, 3))
    nodemat[1,1,1] = e.nodes[1]
    nodemat[3,1,1] = e.nodes[2]
    nodemat[3,3,1] = e.nodes[3]
    nodemat[1,3,1] = e.nodes[4]
    nodemat[1,1,3] = e.nodes[5]
    nodemat[3,1,3] = e.nodes[6]
    nodemat[3,3,3] = e.nodes[7]
    nodemat[1,3,3] = e.nodes[8]
    nodemat[2,1,1] = getcenter(e.boundaries[1].boundaries[4])
    nodemat[3,2,1] = getcenter(e.boundaries[1].boundaries[3])
    nodemat[2,3,1] = getcenter(e.boundaries[1].boundaries[2])
    nodemat[1,2,1] = getcenter(e.boundaries[1].boundaries[1])
    nodemat[2,2,1] = getcenter(e.boundaries[1])
    nodemat[2,1,3] = getcenter(e.boundaries[6].boundaries[1])
    nodemat[3,2,3] = getcenter(e.boundaries[6].boundaries[2])
    nodemat[2,3,3] = getcenter(e.boundaries[6].boundaries[3])
    nodemat[1,2,3] = getcenter(e.boundaries[6].boundaries[4])
    nodemat[2,2,3] = getcenter(e.boundaries[6])
    nodemat[3,1,2] = getcenter(e.boundaries[2].boundaries[2])
    nodemat[1,1,2] = getcenter(e.boundaries[2].boundaries[4])
    nodemat[2,1,2] = getcenter(e.boundaries[2])
    nodemat[1,3,2] = getcenter(e.boundaries[4].boundaries[2])
    nodemat[3,3,2] = getcenter(e.boundaries[4].boundaries[4])
    nodemat[2,3,2] = getcenter(e.boundaries[4])
    nodemat[3,2,2] = getcenter(e.boundaries[3])
    nodemat[1,2,2] = getcenter(e.boundaries[5])
    nodemat[2,2,2] = getcenter(e)
    nodesnearby = getnodeneighborhood(e)
    nodemat = substitute_nodes!(nodemat, nodesnearby)
    cellmat = Array{Element{3,8,6}}(undef, (xstep, ystep, zstep))
    for i in 1:xstep
        for j in 1:ystep
            for k in 1:zstep
                cellmat[i,j,k] = Element([nodemat[i,j,k],nodemat[i + (2 ÷ xstep),j,k],nodemat[i + (2 ÷ xstep),j + (2 ÷ ystep),k],nodemat[i,j + (2 ÷ ystep),k],
                nodemat[i,j,k + (2 ÷ zstep)],nodemat[i + (2 ÷ xstep),j,k + (2 ÷ zstep)],nodemat[i + (2 ÷ xstep),j + (2 ÷ ystep),k + (2 ÷ zstep)],nodemat[i,j + (2 ÷ ystep),k + (2 ÷ zstep)]], e)
            end
        end
    end
    for i in 1:xstep
        for j in 1:ystep
            for k in 1:zstep
                if i < xstep attach_trivially!(cellmat[i,j,k], cellmat[i + 1,j,k]) end
                if j < ystep attach_trivially!(cellmat[i,j,k], cellmat[i,j + 1,k]) end
                if k < zstep attach_trivially!(cellmat[i,j,k], cellmat[i,j,k + 1]) end
            end
        end
    end
    allfaces = collect(Iterators.flatten((x -> x.boundaries).(cellmat)))
    for i in 1:length(e.boundaries)
        side = filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[i] == x), allfaces)
        attach_children!(convert(Array{Element{2,4,4}}, side), e.boundaries[i])
    end

    e.children = reshape(cellmat, length(cellmat))
    # [cellmat[1,1,1],cellmat[2,1,1],cellmat[2,2,1],cellmat[1,2,1],
    # cellmat[1,1,2],cellmat[2,1,2],cellmat[2,2,2],cellmat[1,2,2]]

    return e
end
