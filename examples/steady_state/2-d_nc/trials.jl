using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")

fens1, fes1 = Q4block(0.5, 0.5, 2, 2)
fens2, fes2 = Q4block(0.5, 0.5, 2, 2)
fens2.xyz[:, 1] .+= 0.5
fens3, fes3 = Q4block(0.5, 0.5, 2, 2)
fens3.xyz[:, :] .+= 0.5
fens4, fes4 = Q4block(0.5, 0.5, 2, 2)
fens4.xyz[:, 2] .+= 0.5


boundary_boxes = [
    [0.0,0.0, 0.0,1.0], # left
    [1.0,1.0, 0.0,1.0], # right
    [0.0,1.0, 1.0,1.0], # top
    [0.0,1.0, 0.0,0.0], # bottom
]

edge_fes1 = meshboundary(fes1)
edge_fes2 = meshboundary(fes2)
edge_fes3 = meshboundary(fes3)
edge_fes4 = meshboundary(fes4)
fens_s = [fens1, fens2, fens3, fens4]
interface_fes = extract_interface_fes([edge_fes1, edge_fes2, edge_fes3, edge_fes4], fens_s,  boundary_boxes)

leftx = collect(range(0.5, stop=0.5, length=5))
lefty = collect(range(0.0, stop=1.0, length=5))
fens_left, fes_left = L2blockx2D(leftx, lefty)

bottomx = collect(range(0.0, stop=1.0, length=5))
bottomy = collect(range(0.5, stop=0.5, length=5))
fens_bottom, fes_bottom = L2blockx2D(bottomx, bottomy)

fens_i,fes_i,c = mergemeshes(fens_left, fes_left, fens_bottom, fes_bottom, 1e-8)
fes_i.conn = vcat(fes_i.conn, c.conn)
fes_i.label = vcat(fes_i.label, c.label)

make_union_mesh(interface_fes, fens_s,  fes_i, fens_i, 1; lam_order=0)