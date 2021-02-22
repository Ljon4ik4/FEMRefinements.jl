function Makie.mesh(grid::AdaptiveGrid{2})
    nodes = getnodes(grid)
    coords = [node.x[i] for node in nodes, i in 1:JuAFEM.getdim(grid)]
    cells = getcells(grid)
    connectivity = getproperty.(cells, :nodes)
    elements = [element[i] for element in connectivity, i in 1:4]
    #scene = Makie.scatter(coords[:,1], coords[:,2], scale_plot=false)
    hangingnodes = getncinfo(grid) 
    hangingnodes = nodes[getindex.(hangingnodes,1)]
    hangingcoords = [node.x[i] for node in hangingnodes, i in 1:JuAFEM.getdim(grid)]
    scene = Makie.scatter(hangingcoords[:,1], hangingcoords[:,2], scale_plot=false,markersize=10)
    for (idx, element) in enumerate(eachrow(elements))
        Makie.lines!(scene, coords[element[[1,2,3,4,1]],1], coords[element[[1,2,3,4,1]],2], scale_plot=false)
        midpoint = ((coords[element[1],1] + coords[element[2],1]) / 2, (coords[element[2],2] + coords[element[3],2]) / 2)
        Makie.text!("cell $idx", position=midpoint, textsize=0.04)
    end
    return scene
end

function Makie.mesh(grid::AdaptiveGrid{3})
    nodes = getnodes(grid)
    coords = [node.x[i] for node in nodes, i in 1:JuAFEM.getdim(grid)]
    cells = getcells(grid)
    connectivity = getproperty.(cells, :nodes)
    elements = [element[i] for element in connectivity, i in 1:8]
    #scene = Makie.scatter(coords[:,1], coords[:,2], scale_plot=false)
    ncinfo = getncinfo(grid) 
    hangingnodes = nodes[getindex.(ncinfo,1)]
    hangingcoords = [node.x[i] for node in hangingnodes, i in 1:JuAFEM.getdim(grid)]
    scene = Makie.scatter(hangingcoords[:,1], hangingcoords[:,2],hangingcoords[:,3], scale_plot=false,markersize=10)
    for (idx,hangingcoord) in enumerate(eachrow(hangingcoords))
        number = getindex.(ncinfo,1)[idx]
        Makie.text!("node $number", position=(hangingcoord...,), textsize=0.04)
    end
    
    for (idx, element) in enumerate(eachrow(elements))
        Makie.lines!(scene, coords[element[[1,2,3,4,1]],1], coords[element[[1,2,3,4,1]],2], coords[element[[1,2,3,4,1]],3], scale_plot=false)
        Makie.lines!(scene, coords[element[[1,5]],1], coords[element[[1,5]],2], coords[element[[1,5]],3], scale_plot=false)
        Makie.lines!(scene, coords[element[[2,6]],1], coords[element[[2,6]],2], coords[element[[2,6]],3], scale_plot=false)
        Makie.lines!(scene, coords[element[[3,7]],1], coords[element[[3,7]],2], coords[element[[3,7]],3], scale_plot=false)        
        Makie.lines!(scene, coords[element[[4,8]],1], coords[element[[4,8]],2], coords[element[[4,8]],3], scale_plot=false)        
        Makie.lines!(scene, coords[element[[5,6,7,8,5]],1], coords[element[[5,6,7,8,5]],2], coords[element[[5,6,7,8,5]],3], scale_plot=false)        
        midpoint = ((coords[element[1],1] + coords[element[2],1]) / 2, (coords[element[2],2] + coords[element[3],2]) / 2, (coords[element[4],3] + coords[element[8],3]) / 2)
        Makie.text!("cell $idx", position=midpoint, textsize=0.04)
    end
    return scene
end

function Makie.mesh!(grid::AdaptiveGrid)
    nodes = getnodes(grid)
    coords = [node.x[i] for node in nodes, i in 1:JuAFEM.getdim(grid)]
    cells = getcells(grid)
    connectivity = getproperty.(cells, :nodes)
    elements = [element[i] for element in connectivity, i in 1:4]
    for (coordidx,coord) in enumerate(eachrow(coords))
        Makie.scatter!([coord[1]], [coord[2]], scale_plot=false)
        Makie.text!("$coordidx", position=(coord[1],coord[2], 0.3), textsize=0.04)
    end
    for element in eachrow(elements)
        Makie.lines!(coords[element[[1,2,3,4,1]],1], coords[element[[1,2,3,4,1]],2], scale_plot=false)
    end
end

function Makie.mesh(dh::JuAFEM.AbstractDofHandler, u::Array{T,1}, args...; field::Int=1, scale_plot=false, shading=false, kwargs...) where T
    # extract data
    grid = dh.grid
    nodes = getnodes(grid)
    coords = [node.x[i] for node in nodes, i in 1:JuAFEM.getdim(grid)]
    cells = getcells(grid)
    connectivity = getproperty.(cells, :nodes)
    elements = [element[i] for element in connectivity, i in 1:4]

    # nodes != dofs, hence we need to sort the solution vector node -> dof mapping
    solution = zeros(getnnodes(grid))
    for cell in CellIterator(dh)
        _celldofs = celldofs(cell)
        counter = 1
        for node in getnodes(cell)
            solution[node] = u[_celldofs[counter]]
            counter += 1
        end
    end

    # plot like Guttenberg in 15th century
    triangle_elements = vcat(elements[:,[1,2,3]], elements[:,[3,4,1]])
    return AbstractPlotting.mesh(coords, triangle_elements, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.mesh!(dh::JuAFEM.AbstractDofHandler, u::Array{T,1}, args...; field::Int=1, scale_plot=false, shading=false, kwargs...) where T
    # extract data
    grid = dh.grid
    nodes = getnodes(grid)
    coords = [node.x[i] for node in nodes, i in 1:JuAFEM.getdim(grid)]
    cells = getcells(grid)
    connectivity = getproperty.(cells, :nodes)
    elements = [element[i] for element in connectivity, i in 1:4]

    # nodes != dofs, hence we need to sort the solution vector node -> dof mapping
    solution = zeros(getnnodes(grid))
    for cell in CellIterator(dh)
        _celldofs = celldofs(cell)
        counter = 1
        for node in getnodes(cell)
            solution[node] = u[_celldofs[counter]]
            counter += 1
        end
    end

    # plot like Guttenberg in 15th century
    triangle_elements = vcat(elements[:,[1,2,3]], elements[:,[3,4,1]])
    return AbstractPlotting.mesh!(coords, triangle_elements, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface(dh::JuAFEM.AbstractDofHandler, u::Array{T,1}, args...; field::Int=1, scale_plot=false, shading=false, kwargs...) where T
    @assert JuAFEM.getdim(dh.grid) == 2 "Only 2D solutions supported!"
    C = getcelltype(dh.grid)
    nodes = getnodes(dh.grid)
    cells = getcells(dh.grid)
    coords = [node.x[i] for node in nodes, i in 1:2]
    connectivity = getproperty.(cells, :nodes)
    N = length(vertices(cells[1]))
    elements = [element[i] for element in connectivity, i in 1:N]
    solution = zeros(getnnodes(dh.grid))
    for cell in CellIterator(dh)
        _celldofs = celldofs(cell)
        counter = 1
        for node in getnodes(cell)
            solution[node] = u[_celldofs[counter]]
            counter += 1
        end
    end
    points = [AbstractPlotting.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(coords))]
    triangle_elements = vcat(elements[:,[1,2,3]], elements[:,[3,4,1]])
    return AbstractPlotting.mesh(points, triangle_elements, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface!(dh::JuAFEM.AbstractDofHandler, u::Array{T,1}, args...; field::Int=1, scale_plot=false, shading=false, kwargs...) where T
    @assert JuAFEM.getdim(dh.grid) == 2 "Only 2D solutions supported!"
    C = getcelltype(dh.grid)
    nodes = getnodes(dh.grid)
    cells = getcells(dh.grid)
    coords = [node.x[i] for node in nodes, i in 1:2]
    connectivity = getproperty.(cells, :nodes)
    N = length(vertices(cells[1]))
    elements = [element[i] for element in connectivity, i in 1:N]
    solution = zeros(getnnodes(dh.grid))
    for cell in CellIterator(dh)
        _celldofs = celldofs(cell)
        counter = 1
        for node in getnodes(cell)
            solution[node] = u[_celldofs[counter]]
            counter += 1
        end
    end
    points = [AbstractPlotting.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(coords))]
    triangle_elements = vcat(elements[:,[1,2,3]], elements[:,[3,4,1]])
    return AbstractPlotting.mesh!(points, triangle_elements, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function plot_solution2D(dh::JuAFEM.AbstractDofHandler, u; plotgrid=false)
    scene = Makie.mesh(dh, u)
    if plotgrid
        Makie.mesh!(dh.grid)
    end
    Makie.center!(scene)
    ls = Makie.colorlegend(:viridis, minimum(u):maximum(u) / 2:maximum(u), camera=Makie.campixel!, raw=true, width=(30, 300))
    scene_final = Makie.vbox(scene, ls)
end

function plot_solution3D(dh::JuAFEM.AbstractDofHandler, u; plotgrid=false)
    scene = Makie.surface(dh, u)
    if plotgrid
        Makie.mesh!(dh.grid)
    end
    Makie.center!(scene)
    ls = Makie.colorlegend(:viridis, minimum(u):maximum(u) / 2:maximum(u), camera=Makie.campixel!, raw=true, width=(30, 300))
    scene_final = Makie.vbox(scene, ls)
end
