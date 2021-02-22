
using JuAFEM, SparseArrays
using JuAFEM:Vec
using StaticArrays


using FEMRefinements

v1 = Vec((2.0, 0.0))
v2 = Vec((0.0, 2.0))

grid = generate_grid2D(1, 2, v1, v2);
# FEMRefinements.refine!(grid, [2])
# FEMRefinements.refine!(grid, [1])
FEMRefinements.refine!(grid, collect(1:getncells(grid)))
FEMRefinements.refine!(grid, [1,4,8])
# FEMRefinements.refine!(grid, collect(1:getncells(grid)))
# FEMRefinements.refine!(grid, [1,4])
# FEMRefinements.refine!(grid, collect(1:getncells(grid)))
# FEMRefinements.refine!(grid, [getncells(grid) - 1,getncells(grid)])
# FEMRefinements.refine!(grid, [getncells(grid)])
# FEMRefinements.refine!(grid, [getncells(grid) - 1,getncells(grid)])
# FEMRefinements.refine!(grid, [15])


dim = 2
ip = Lagrange{dim,RefCube,1}()
qr = QuadratureRule{dim,RefCube}(2)
cellvalues = CellScalarValues(qr, ip);
field1 = Field(:u, ip, 1)
dh = MixedDofHandler(grid)
push!(dh, JuAFEM.FieldHandler([field1], Set((1:getncells(grid)))));
dh, vertexdict, _, _ = JuAFEM.__close!(dh);
K = create_sparsity_pattern(dh);

# using UnicodePlots
# fill!(K.nzval, 1.0)
# spy(K; height = 15)


ch = ConstraintHandler(dh);

∂Ω = union(getfaceset.((grid,), ["left"])...);

dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)

function doassemble(cellvalues::CellScalarValues{dim}, K::SparseMatrixCSC, dh::JuAFEM.MixedDofHandler) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)

    @inbounds for cell in CellIterator(dh)

        fill!(Ke, 0)
        fill!(fe, 0)


        reinit!(cellvalues, cell)


        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            for i in 1:n_basefuncs
                v  = shape_value(cellvalues, q_point, i)
                ∇v = shape_gradient(cellvalues, q_point, i)
                fe[i] += v * dΩ
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += (∇v ⋅ ∇u) * dΩ
                end
            end
        end

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
    return K, f
end

K, f = doassemble(cellvalues, K, dh);
ncinfo = getncinfo(grid)
hangingnodes = getindex.(ncinfo, 1)
mastercells = getindex.(ncinfo, 2)
slavecells = getindex.(ncinfo, 3)
ncdofs = get.((vertexdict[1],), hangingnodes, NaN) ## ist das mapping nicht andersrum?
if isa(dh, JuAFEM.MixedDofHandler)
    ncdofs = vcat(ncdofs...)
end
_ndofs = ndofs(dh)

I, J, V = Int64[], Int64[], Float64[]

ip = Lagrange{2,RefCube,1}()
tdof = 0 # TODO ++ if conforming dof
for dof in 1:_ndofs
    ncidx = findfirst(x -> x == dof, ncdofs)
    if ncidx !== nothing
        occurence = findall(x -> x == mastercells[ncidx], mastercells)
        # celldepth =
        if length(occurence) == 1
            mastercellidx = getindex.(mastercells, 1)[ncidx]
            mastercoords = getcoordinates(grid, mastercellidx)
            masterdofs = celldofs(dh, mastercellidx)
            hangingnode = grid.nodes[dof]
            facecoords = JuAFEM.faces(grid.cells[mastercellidx])
            masterfacedofs, masterdofs = facedofs(dh, mastercellidx, getindex.(mastercells, 2)[ncidx])
            local_coordinate = find_local_coordinate(Lagrange{2,RefCube,1}(), mastercoords, hangingnode.x)
            quad_rule = QuadratureRule{2,RefCube,Float64}([1], [local_coordinate])
            cv = CellScalarValues(quad_rule, ip)
            reinit!(cv, mastercoords)
            for i in 1:length(masterfacedofs)
                Pvalue = JuAFEM.shape_value(cv, 1, findfirst(x -> x == masterfacedofs[i], masterdofs))
                push!(I, dof)
                push!(J, masterfacedofs[i])
                push!(V, Pvalue)
            end
        end
    else
        push!(I, dof)
        push!(J, dof)
        push!(V, 1.0)
    end
end

P = sparse(I, J, V, _ndofs, _ndofs)

JuAFEM.add!(ch, dbc,ncdofs);
JuAFEM.close!(ch,ncdofs)
JuAFEM.update!(ch, ncdofs, 0.0);
P = P[:, 1:end .∉ [ncdofs]]
A = (P' * K) * P
f = P' * f
apply!(A, f, ch)
u = A \ f;

# vtk_grid("heat_equation", dh) do vtk
#    vtk_point_data(vtk, dh, P*u)
# end

# using Test                        #src
# @test norm(u) ≈ 3.307743912641305 #src

FEMRefinements.plot_solution3D(dh,P * u; plotgrid=true)
