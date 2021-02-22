# maybe rewrite the 2d generator to the style of the 3d one
function generate_grid2D(cols::Int, rows::Int, v1::Vec{dim,T}, v2::Vec{dim,T}) where {dim,T <: Real}
    nodemat = Matrix{Node{dim,T}}(undef, (cols + 1, rows + 1))
    xstep = 1.0 / cols
    ystep = 1.0 / rows
    for i in 0:cols
        for j in 0:rows
            nodemat[i + 1,j + 1] = Node{dim,T}(i * xstep * v1 + j * ystep * v2)
        end
    end
    nodes = reshape(nodemat, length(nodemat))

    cellmat = Matrix{Element}(undef, (cols, rows))
    for i in 1:cols
        for j in 1:rows
            cellmat[i,j] = Element([nodemat[i, j] , nodemat[i + 1,j] ,nodemat[i + 1, j + 1], nodemat[i, j + 1]], nothing)
        end
    end
    cells = reshape(cellmat, length(cellmat))


    vertfacemat = Matrix{Element}(undef, (cols + 1, rows))
    for i in 2:cols
        for j in 1:rows
            vertfacemat[i,j] = attach_trivially!(cellmat[i - 1,j], cellmat[i,j])
        end
    end
    for j in 1:rows
        vertfacemat[1,j] = cellmat[1,j].boundaries[4]
        vertfacemat[cols + 1,j] = cellmat[cols,j].boundaries[2]
    end



    horfacemat = Matrix{Element}(undef, (cols, rows + 1))
    for i in 1:cols
        for j in 2:rows
            horfacemat[i,j] = attach_trivially!(cellmat[i,j - 1], cellmat[i,j])
        end
    end
    for i in 1:cols
        horfacemat[i,1] = cellmat[i,1].boundaries[1]
        horfacemat[i,rows + 1] = cellmat[i,rows].boundaries[3]
    end
    faces = vcat(reshape(vertfacemat, length(vertfacemat)), reshape(horfacemat, length(horfacemat)))

    cellsets = Dict{String,Set{Element}}()
    cellsets["all"] = Set(cells)
    nodesets = Dict{String,Set{Node}}()
    nodesets["corners"] = Set([nodemat[1,1], nodemat[cols + 1,1], nodemat[cols + 1,rows + 1],nodemat[1,rows + 1]])
    facesets = Dict{String,Set{Element}}()
    facesets["bottom"] = Set(horfacemat[:,1])
    facesets["left"] = Set(vertfacemat[1,:])
    facesets["top"] = Set(horfacemat[:,rows + 1])
    facesets["right"] = Set(vertfacemat[cols + 1,:])
    rootcells = cells
    grid = AdaptiveGrid{dim,T}(rootcells, cells, faces, nodes, cellsets, nodesets, facesets)
    return grid
end




function generate_grid3D(cols::Int, rows::Int, depth::Int, v1::Vec{dim,T}, v2::Vec{dim,T}, v3::Vec{dim,T}) where {dim,T <: Real}
    nodearr = Array{Node{dim,T}}(undef, (cols + 1, rows + 1, depth + 1))
    xstep = 1.0 / cols
    ystep = 1.0 / rows
    zstep = 1.0 / depth
    for i in 0:cols
        for j in 0:rows
            for k in 0:depth
                nodearr[i + 1,j + 1,k + 1] = Node{dim,T}(i * xstep * v1 + j * ystep * v2 + k * zstep * v3)
            end
        end
    end
    nodes = reshape(nodearr, length(nodearr))

    cellarr = Array{Element}(undef, (cols, rows, depth))
    cubenumberings = [[0, 0, 0], [1, 0, 0], [1,1,0], [0,1,0], [0, 0, 1], [1, 0, 1], [1,1,1], [0,1,1]]
    for i in 1:cols
        for j in 1:rows
            for k in 1:depth
                cellarr[i,j,k] = Element((x -> nodearr[CartesianIndex(Tuple([i,j,k] + x))]).(cubenumberings), nothing)
            end
        end
    end
    cells = reshape(cellarr, length(cellarr))
    for i in 1:cols
        for j in 1:rows
            for k in 1:depth
                if i < cols attach_trivially!(cellarr[i,j,k], cellarr[i + 1,j,k]) end
                if j < rows attach_trivially!(cellarr[i,j,k], cellarr[i,j + 1,k]) end
                if k < depth attach_trivially!(cellarr[i,j,k], cellarr[i,j,k + 1]) end
            end
        end
    end

    cellsets = Dict{String,Set{Element}}()
    cellsets["all"] = Set(cells)
    nodesets = Dict{String,Set{Node}}()
    nodesets["corners"] = Set((x -> nodearr[CartesianIndex(Tuple([cols + 1, rows + 1, depth + 1].^x))]).(cubenumberings))
    facesets = Dict{String,Set{Element}}()
    rootcells = cells
    grid = AdaptiveGrid{dim,T}(rootcells, cells, [], nodes, cellsets, nodesets, facesets)
    grid = update_leafmesh!(grid)
    grid.facesets["bottom"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[1] == x), grid.faces)) # labels possibly misleading
    grid.facesets["left"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[2] == x), grid.faces))
    grid.facesets["front"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[3] == x), grid.faces))
    grid.facesets["right"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[4] == x), grid.faces))
    grid.facesets["back"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[5] == x), grid.faces))
    grid.facesets["top"] = Set(filter(x -> (length(x.attachedto) == 0 && x.boundaryof.boundaries[6] == x), grid.faces))
    return grid
end
