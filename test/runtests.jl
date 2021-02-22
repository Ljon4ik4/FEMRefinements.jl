using Test
using JuAFEM
using JuAFEM:Vec

include("../src/FEMRefinements.jl")
using .FEMRefinements

@testset "grid2D_generator" begin
    v1 = Vec((1.0, 0.0))
    v2 = Vec((0.0, 1.0))
    grid = generate_grid2D(2, 2, v1, v2)
    @test getncells(grid) == 4
    @test JuAFEM.getcelltype(grid) == Element{2,4,4}
    @test getnnodes(grid) == 9
    @test JuAFEM.nnodes_per_cell(grid, 1) == 4
    @test getcells(grid)[1].nodes[1] == 1
    @test JuAFEM.n_faces_per_cell(grid) == 4
    @test getncells(grid) == length(grid.cells)
end

@testset "grid3D_generator" begin
    v1 = Vec((1.0, 0.0, 0.0))
    v2 = Vec((0.0, 1.0, 0.0))
    v3 = Vec((0.0, 0.0, 1.0))
    grid3D = generate_grid3D(2, 2, 2, v1, v2, v3)
    @test getncells(grid3D) == 8
    @test JuAFEM.getcelltype(grid3D) == Element{3,8,6}
    @test getnnodes(grid3D) == 27
    @test JuAFEM.nnodes_per_cell(grid3D, 1) == 8
    @test getcells(grid3D)[1].nodes[1] == 1
    @test JuAFEM.n_faces_per_cell(grid3D) == 6
    @test getncells(grid3D) == length(grid3D.cells)
end

@testset "JuAFEM Grid utils" begin
    v1 = Vec((1.0, 0.0))
    v2 = Vec((0.0, 1.0))
    grid = generate_grid2D(2, 2, v1, v2)
    addcellset!(grid, "cell_set", [1]);
    node_set = Set(1:getnnodes(grid))
    addnodeset!(grid, "node_set", node_set)
    @test getnodesets(grid) == Dict("corners" => Set([7, 9, 3, 1]), "node_set" => node_set)
    @test getnodes(grid, [1]) == [getnodes(grid, 1)]
    @test length(getnodes(grid, "node_set")) == 9
    @test collect(getcoordinates(getnodes(grid, 5)).data) ≈ [0.5, 0.5]
    @test getcells(grid, "cell_set") == [getcells(grid, 1)]
    f(x) = Tensor{1,1,Float64}((1 + x[1]^2 + 2x[2]^2,))
    values = compute_vertex_values(grid, f)
    @test f([0.0, 0.0]) == values[1]
    @test f([0.5, 0.5]) == values[5]
    @test f([1.0, 1.0]) == values[9]
    @test compute_vertex_values(grid, collect(1:9), f) == values
    @test length(compute_vertex_values(grid, "node_set", f)) == 9
    @test getcoordinates(grid, 1)[1] != getcoordinates(grid, 1)[2]
    @test getcoordinates(grid, 1)[1] != getcoordinates(grid, 1)[3]
    @test getcoordinates(grid, 1)[1] != getcoordinates(grid, 1)[4]
end

@testset "Refinements2D" begin
    v1 = Vec((1.0, 0.0))
    v2 = Vec((0.0, 1.0))
    grid = generate_grid2D(2, 2, v1, v2)
    grid_to_compare = generate_grid2D(2, 2, v1, v2)
    refine!(grid, [1])
    @test length(grid.cells) == (length(grid_to_compare.cells) + 3)
    @test length(grid.faces) == (length(grid_to_compare.faces) + 6)
    @test length(grid.nodes) == (length(grid_to_compare.nodes) + 5)
    @test JuAFEM.nnodes_per_cell(grid) == 4
    refine!(grid_to_compare, [1])
    @test getnodesets(grid) == getnodesets(grid_to_compare)
    @test getfacesets(grid) == getfacesets(grid_to_compare)
    @test getcellsets(grid) == getcellsets(grid_to_compare)
    i = findfirst(x -> ( (x.parent !== nothing) && (Vec((0.5, 0.5)) ∈ (n -> n.x).(vertices(x)))), grid.cells)
    refine!(grid, [i])
    @test length(grid.cells) == (length(grid_to_compare.cells) + 3)
    @test length(grid.faces) == (length(grid_to_compare.faces) + 4)
    @test length(grid.nodes) == (length(grid_to_compare.nodes) + 5)
    @test JuAFEM.nnodes_per_cell(grid) == 4
    for n in grid.nodes
        @test Base.unique(grid.nodes) == grid.nodes
    end

    v1 = Vec((2.0, 0.0))
    v2 = Vec((0.0, 2.0))
    gridd = generate_grid2D(1, 1, v1, v2);
    FEMRefinements.refine!(gridd, collect(1:getncells(gridd)))
    FEMRefinements.refine!(gridd, [1])
    a = [[1.0, 0.5],[0.5, 1.0]]
    z = 1
    ncinfo = getncinfo(gridd)
    for n in ncinfo
        @test gridd.nodes[n[1]].x == a[z]
        z = z + 1
    end
end

@testset "Derefinements2D" begin
    v1 = Vec((1.0, 0.0))
    v2 = Vec((0.0, 1.0))
    grid = generate_grid2D(2, 2, v1, v2)
    grid2compare = generate_grid2D(2, 2, v1, v2)
    refine!(grid, [1])
    a = findall((x -> x.parent !== nothing), grid.cells)
    derefine!(grid, a)
    @test length(grid.cells) == length(grid2compare.cells)
    @test length(grid.faces) == length(grid2compare.faces)
    @test length(grid.nodes) == length(grid2compare.nodes)
    @test JuAFEM.nnodes_per_cell(grid) == 4
    @test getcellsets(grid) == getcellsets(grid2compare)
end

@testset "Refinements3D" begin
    v1 = Vec((1.0, 0.0, 0.0))
    v2 = Vec((0.0, 1.0, 0.0))
    v3 = Vec((0.0, 0.0, 1.0))
    grid3D = generate_grid3D(2, 2, 2, v1, v2, v3)
    comparison = generate_grid3D(2, 2, 2, v1, v2, v3)
    refine!(grid3D, [1])
    @test length(grid3D.cells) == (length(comparison.cells) + 7)
    @test length(grid3D.faces) == (length(comparison.faces) + 21)
    @test length(grid3D.nodes) == (length(comparison.nodes) + 19)
    @test JuAFEM.nnodes_per_cell(grid3D) == 8
    for n in grid3D.nodes
        @test Base.unique(grid3D.nodes) == grid3D.nodes
    end
    refine!(comparison, [1])
    @test getnodesets(grid3D) == getnodesets(comparison)
    @test getfacesets(grid3D) == getfacesets(comparison)
    @test getcellsets(grid3D) == getcellsets(comparison)
    i = findfirst(x -> ( (x.parent !== nothing) && (Vec((0.5, 0.5, 0.5)) ∈ (n -> n.x).(vertices(x)))), grid3D.cells)
    refine!(grid3D, [i])
    @test length(grid3D.cells) == (length(comparison.cells) + 7)
    @test length(grid3D.faces) == (length(comparison.faces) + 12)
    @test length(grid3D.nodes) == (length(comparison.nodes) + 19)
    @test JuAFEM.nnodes_per_cell(grid3D) == 8

    v1 = Vec((2.0, 0.0, 0.0))
    v2 = Vec((0.0, 2.0, 0.0))
    v3 = Vec((0.0, 0.0, 2.0))
    grid = generate_grid3D(1, 1, 1, v1, v2, v3);
    FEMRefinements.refine!(grid, collect(1:getncells(grid)))
    FEMRefinements.refine!(grid, [1])
    a = [[1.0, 0.5, 0.0],[1.0, 0.0, 0.5],[1.0, 1.0, 0.5],[1.0, 0.5, 0.0],[1.0, 0.0, 0.5],[1.0, 0.5, 0.5],[1.0, 1.0, 0.5],[1.0, 0.5, 1.0],[1.0, 0.5, 1.0],[0.5, 1.0, 0.0],[0.5, 1.0, 0.0],[0.5, 1.0, 0.5],[0.0, 1.0, 0.5],[1.0, 1.0, 0.5],[0.5, 1.0, 1.0],[1.0, 1.0, 0.5],[0.0, 1.0, 0.5],[0.5, 1.0, 1.0],[0.5, 0.0, 1.0],[0.5, 0.5, 1.0],[0.0, 0.5, 1.0],[1.0, 0.5, 1.0],[0.5, 1.0, 1.0],[0.5, 0.0, 1.0],[1.0, 0.5, 1.0],[0.5, 1.0, 1.0],[0.0, 0.5, 1.0]]
    z = 1
    ncinfo = getncinfo(grid)
    for n in ncinfo
        @test grid.nodes[n[1]].x == a[z]
        z = z + 1
    end
end

@testset "Derefinements3D" begin
    v1 = Vec((1.0, 0.0, 0.0))
    v2 = Vec((0.0, 1.0, 0.0))
    v3 = Vec((0.0, 0.0, 1.0))
    grid3D = generate_grid3D(2, 2, 2, v1, v2, v3)
    compare_it = generate_grid3D(2, 2, 2, v1, v2, v3)
    refine!(grid3D, [1])
    a = findall((x -> x.parent !== nothing), grid3D.cells)
    derefine!(grid3D, a)
    @test length(grid3D.cells) == length(compare_it.cells)
    @test length(grid3D.faces) == length(compare_it.faces)
    @test length(grid3D.nodes) == length(compare_it.nodes)
    @test JuAFEM.nnodes_per_cell(grid3D) == 8
    @test getcellsets(grid3D) == getcellsets(compare_it)
end
