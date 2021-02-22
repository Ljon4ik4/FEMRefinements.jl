function JuAFEM.MixedDofHandler(grid::AdaptiveGrid{dim,T}) where {dim,T}
    JuAFEM.MixedDofHandler{dim,T,typeof(grid)}(FieldHandler[], JuAFEM.CellVector(Int[],Int[],Int[]), JuAFEM.CellVector(Int[],Int[],Int[]), JuAFEM.CellVector(Vec{dim,T}[],Int[],Int[]), JuAFEM.ScalarWrapper(false), grid, JuAFEM.ScalarWrapper(-1))
end

function JuAFEM.__close!(dh::MixedDofHandler{dim,T,AdaptiveGrid{dim,T}}) where {dim, T}

    @assert !JuAFEM.isclosed(dh)
    field_names = JuAFEM.getfieldnames(dh)  # all the fields in the problem
    numfields =  length(field_names)

    # Create dicts that stores created dofs
    # Each key should uniquely identify the given type
    vertexdicts = [Dict{Int, Array{Int}}() for _ in 1:numfields]
    edgedicts = [Dict{Tuple{Int,Int}, Array{Int}}() for _ in 1:numfields]
    facedicts = [Dict{Tuple{Int,Int}, Array{Int}}() for _ in 1:numfields]
    celldicts = [Dict{Int, Array{Int}}() for _ in 1:numfields]

    # Set initial values
    nextdof = 1  # next free dof to distribute

    append!(dh.cell_dofs.offset, zeros(getncells(dh.grid)))
    append!(dh.cell_dofs.length, zeros(getncells(dh.grid)))

    @debug "\n\nCreating dofs\n"
    for fh in dh.fieldhandlers
        # sort the cellset since we want to loop through the cells in a fixed order
        cellnumbers = sort(collect(fh.cellset))
        nextdof = JuAFEM._close!(
            dh,
            cellnumbers,
            JuAFEM.getfieldnames(fh),
            JuAFEM.getfielddims(fh),
            JuAFEM.getfieldinterpolations(fh),
            nextdof,
            vertexdicts,
            edgedicts,
            facedicts,
            celldicts)
    end
    dh.ndofs[] = maximum(dh.cell_dofs.values)
    dh.closed[] = true

    #Create cell_nodes and cell_coords (similar to cell_dofs)
    push!(dh.cell_nodes.offset, 1)
    push!(dh.cell_coords.offset, 1)
    for cell in dh.grid.cells
        for nodeid in JuAFEM.vertices(cell, dh.grid)
            push!(dh.cell_nodes.values, nodeid)
            push!(dh.cell_coords.values, dh.grid.nodes[nodeid].x)
        end
        push!(dh.cell_nodes.offset, length(dh.cell_nodes.values)+1)
        push!(dh.cell_coords.offset, length(dh.cell_coords.values)+1)
        push!(dh.cell_nodes.length, length(cell.nodes))
        push!(dh.cell_coords.length, length(cell.nodes))
    end

    return dh, vertexdicts, edgedicts, facedicts

end

function JuAFEM._close!(dh::JuAFEM.MixedDofHandler{dim, T, AdaptiveGrid{dim, T}}, cellnumbers, field_names, field_dims, field_interpolations, nextdof, vertexdicts, edgedicts, facedicts, celldicts) where {dim, T}

    ip_infos = JuAFEM.InterpolationInfo[]
    for interpolation in field_interpolations
        ip_info = JuAFEM.InterpolationInfo(interpolation)
        # these are not implemented yet (or have not been tested)
        @assert(ip_info.nvertexdofs <= 1)
        @assert(ip_info.nedgedofs <= 1)
        @assert(ip_info.nfacedofs <= 1)
        @assert(ip_info.ncelldofs <= 1)  # not tested but probably works
        push!(ip_infos, ip_info)
    end

    # loop over all the cells, and distribute dofs for all the fields
    for ci in cellnumbers
        dh.cell_dofs.offset[ci] = length(dh.cell_dofs.values)+1

        cell = dh.grid.cells[ci]
        cell_dofs = Int[]  # list of global dofs for each cell
        @debug "Creating dofs for cell #$ci"

        for fi in 1:length(field_names)
            @debug "\tfield: $(field_names[fi])"
            ip_info = ip_infos[fi]

            if ip_info.nvertexdofs > 0
                vertices = JuAFEM.vertices(cell, dh.grid)
                nextdof = JuAFEM.add_vertex_dofs(cell_dofs, cell, vertexdicts[fi], field_dims[fi], ip_info.nvertexdofs, nextdof, vertices)
            end

            if ip_info.nedgedofs > 0 && dim == 3 #Edges only in 3d
                edges = JuAFEM.edges(cell, dh.grid)
                nextdof = JuAFEM.add_edge_dofs(cell_dofs, cell, edgedicts[fi], field_dims[fi], ip_info.nedgedofs, nextdof, edges)
            end

            if ip_info.nfacedofs > 0 && (ip_info.dim == dim)
                faces = JuAFEM.faces(cell, dh.grid)
                nextdof = JuAFEM.add_face_dofs(cell_dofs, cell, facedicts[fi], field_dims[fi], ip_info.nfacedofs, nextdof, faces)
            end

            if ip_info.ncelldofs > 0
                nextdof = JuAFEM.add_cell_dofs(cell_dofs, ci, celldicts[fi], field_dims[fi], ip_info.ncelldofs, nextdof)
            end

        end
        # after done creating dofs for the cell, push them to the global list
        push!(dh.cell_dofs.values, cell_dofs...)
        dh.cell_dofs.length[ci] = length(cell_dofs)

        @debug "Dofs for cell #$ci:\n\t$cell_dofs"
    end # cell loop
    return nextdof
end

function JuAFEM.add_vertex_dofs(cell_dofs, cell, vertexdict, field_dim, nvertexdofs, nextdof, vertices)
    for vertex in vertices
        @debug "\tvertex #$vertex"
        nextdof, dofs = JuAFEM.get_or_create_dofs!(nextdof, field_dim, dict=vertexdict, key=vertex)
        push!(cell_dofs, dofs...)
    end
    return nextdof
end

function JuAFEM.add_face_dofs(cell_dofs, cell, facedict, field_dim, nfacedofs, nextdof, faces)
    @assert nfacedofs == 1 "Currently only supports interpolations with nfacedofs = 1"

    for face in faces
        sface = JuAFEM.sortface(face)
        @debug "\tface #$sface"
        nextdof, dofs = JuAFEM.get_or_create_dofs!(nextdof, field_dim, dict=facedict, key=sface)
        push!(cell_dofs, dofs...)
    end
    return nextdof
end

function JuAFEM.add_edge_dofs(cell_dofs, cell, edgedict, field_dim, nedgedofs, nextdof, edges)
    @assert nedgedofs == 1 "Currently only supports interpolations with nedgedofs = 1"
    for edge in edges
        sedge, dir = JuAFEM.sortedge(edge)
        @debug "\tedge #$sedge"
        nextdof, dofs = JuAFEM.get_or_create_dofs!(nextdof, field_dim, dict=edgedict, key=sedge)
        push!(cell_dofs, dofs...)
    end
    return nextdof
end

##not sure if needed
#function add_cell_dofs(cell_dofs, cell, celldict, field_dim, ncelldofs, nextdof)
#    for celldof in 1:ncelldofs
#        @debug "\tcell #$cell"
#        nextdof, dofs = get_or_create_dofs!(nextdof, field_dim, dict=celldict, key=cell)
#        push!(cell_dofs, dofs...)
#    end
#    return nextdof
#end

function JuAFEM.add!(ch::ConstraintHandler, dbc::Dirichlet, ncdofs::Vector)
    JuAFEM.dbc_check(ch, dbc)
    field_idx = JuAFEM.find_field(ch.dh, dbc.field_name)
    # Extract stuff for the field
    interpolation = JuAFEM.getfieldinterpolation(ch.dh, field_idx)#ch.dh.field_interpolations[field_idx]
    field_dim = JuAFEM.getfielddim(ch.dh, field_idx)#ch.dh.field_dims[field_idx] # TODO: I think we don't need to extract these here ...
    bcvalue = JuAFEM.getbcvalue(ch.dh, field_idx)
    JuAFEM._add!(ch, dbc, dbc.faces, interpolation, field_dim, JuAFEM.field_offset(ch.dh, dbc.field_name), bcvalue, ncdofs)
    return ch
end

function JuAFEM._add!(ch::ConstraintHandler, dbc::Dirichlet, bcfaces::Set{Tuple{Int,Int}}, interpolation::Interpolation, field_dim::Int, offset::Int, bcvalue::JuAFEM.BCValues, ncdofs::Vector)
    # calculate which local dof index live on each face
    # face `i` have dofs `local_face_dofs[local_face_dofs_offset[i]:local_face_dofs_offset[i+1]-1]
    local_face_dofs = Int[]
    local_face_dofs_offset = Int[1]
    for (i, face) in enumerate(JuAFEM.faces(interpolation))
        for fdof in face, d in 1:field_dim
            if d ∈ dbc.components # skip unless this component should be constrained
                push!(local_face_dofs, (fdof-1)*field_dim + d + offset)
            end
        end
        push!(local_face_dofs_offset, length(local_face_dofs) + 1)
    end
    JuAFEM.copy!!(dbc.local_face_dofs, local_face_dofs)
    JuAFEM.copy!!(dbc.local_face_dofs_offset, local_face_dofs_offset)

    # loop over all the faces in the set and add the global dofs to `constrained_dofs`
    constrained_dofs = Int[]
    #_celldofs = fill(0, ndofs_per_cell(ch.dh))
    for (cellidx, faceidx) in bcfaces
        _celldofs = fill(0, ndofs_per_cell(ch.dh, cellidx))
        celldofs!(_celldofs, ch.dh, cellidx) # extract the dofs for this cell
        r = local_face_dofs_offset[faceidx]:(local_face_dofs_offset[faceidx+1]-1)
        dofs = _celldofs[local_face_dofs[r]] 
        for dof in dofs
            append!(constrained_dofs, dof - sum(ncdofs .<= dof )) # TODO: for-loop over r and simply push! to ch.prescribed_dofs
        end
        @debug println("adding dofs $(_celldofs[local_face_dofs[r]]) to dbc")
    end

    # save it to the ConstraintHandler
    push!(ch.dbcs, dbc)
    push!(ch.bcvalues, bcvalue)
    append!(ch.prescribed_dofs, constrained_dofs)
end

function JuAFEM.close!(ch::ConstraintHandler, ncdofs::Vector)
    fdofs = setdiff(1:ndofs(ch.dh), ch.prescribed_dofs)
    JuAFEM.copy!!(ch.free_dofs, fdofs)
    JuAFEM.copy!!(ch.prescribed_dofs, unique(ch.prescribed_dofs)) # for v0.7: unique!(ch.prescribed_dofs)
    sort!(ch.prescribed_dofs) # YOLO
    fill!(resize!(ch.values, length(ch.prescribed_dofs)), NaN)
    for i in 1:length(ch.prescribed_dofs)
        key = ch.prescribed_dofs[i]
        #ch.dofmapping[key-sum(ncdofs .<= key)] = i
        ch.dofmapping[key] = i
    end
    ch.closed[] = true
    return ch
end

# Updates the DBC's to the current time `time`
function JuAFEM.update!(ch::ConstraintHandler,ncdofs::Vector, time::Real=0.0)
    @assert ch.closed[]
    for (i,dbc) in enumerate(ch.dbcs)
        field_idx = JuAFEM.find_field(ch.dh, dbc.field_name)
        # Function barrier
        JuAFEM._update!(ch.values, dbc.f, dbc.faces, dbc.field_name, dbc.local_face_dofs, dbc.local_face_dofs_offset,
                 dbc.components, ch.dh, ch.bcvalues[i], ch.dofmapping, convert(Float64, time), ncdofs)
    end
end

# for faces
function JuAFEM._update!(values::Vector{Float64}, f::Function, faces::Set{Tuple{Int,Int}}, field::Symbol, local_face_dofs::Vector{Int}, local_face_dofs_offset::Vector{Int},
                  components::Vector{Int}, dh::JuAFEM.AbstractDofHandler, facevalues::JuAFEM.BCValues,
                  dofmapping::Dict{Int,Int}, time::T, ncdofs::Vector) where {T}

    dim = JuAFEM.getdim(dh.grid)
    _tmp_cellid = first(faces)[1]

    N = JuAFEM.nnodes_per_cell(dh.grid, _tmp_cellid)
    xh = zeros(Vec{dim, T}, N) # pre-allocate
    _celldofs = fill(0, ndofs_per_cell(dh, _tmp_cellid))

    for (cellidx, faceidx) in faces
        JuAFEM.cellcoords!(xh, dh, cellidx)
        JuAFEM.celldofs!(_celldofs, dh, cellidx) # update global dofs for this cell

        # no need to reinit!, enough to update current_face since we only need geometric shape functions M
        facevalues.current_face[] = faceidx

        # local dof-range for this face
        r = local_face_dofs_offset[faceidx]:(local_face_dofs_offset[faceidx+1]-1)
        counter = 1

        for location in 1:getnquadpoints(facevalues)
            x = spatial_coordinate(facevalues, location, xh)
            bc_value = f(x, time)
            @assert length(bc_value) == length(components)

            for i in 1:length(components)
                # find the global dof
                globaldof = _celldofs[local_face_dofs[r[counter]]]
                globaldof -= sum(ncdofs .<= globaldof)
                
                counter += 1

                dbc_index = dofmapping[globaldof]
                values[dbc_index] = bc_value[i]
                @debug println("prescribing value $(bc_value[i]) on global dof $(globaldof)")
            end
        end
    end
end

function facedofs(dh::JuAFEM.AbstractDofHandler, ele::Int, localface::Int)
    eledofs = celldofs(dh,ele)
    return eledofs[[JuAFEM.faces(Lagrange{2,RefCube,1}())[localface]...]], eledofs ## TODO decluster this 
end

function find_local_coordinate(interpolation, cell_coordinates, global_coordinate)
    dim = length(global_coordinate)
    local_guess = zero(Vec{dim})
    n_basefuncs = getnbasefunctions(interpolation)
    max_iters = 10
    tol_norm = 1e-10
    for iter in 1:10
        if iter == max_iters
            error("did not find a local coordinate")
        end
        N = JuAFEM.value(interpolation, local_guess)

        global_guess = zero(Vec{dim})
        for j in 1:n_basefuncs
            global_guess += N[j] * cell_coordinates[j]
        end
        residual = global_guess - global_coordinate
        if norm(residual) <= tol_norm
            break
        end
        dNdξ = JuAFEM.derivative(interpolation, local_guess)
        J = zero(Tensor{2, 2})
        for j in 1:n_basefuncs
            J += cell_coordinates[j] ⊗ dNdξ[j]
        end
        local_guess -= inv(J) ⋅ residual
    end
    return local_guess
end
