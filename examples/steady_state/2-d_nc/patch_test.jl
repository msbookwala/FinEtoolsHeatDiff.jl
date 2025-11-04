using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")

kappa = [1.0 0; 0 1.0] 
material = MatHeatDiff(kappa)

width1 = 1.0
height1 = 2.0
n1x = 1
n1y = 2
fens1, fes1 = T3block(width1, height1, n1x, n1y)

geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

box1 = [0.0,0.0,0.0,0.0]
dbc_nodes1 = selectnode(fens1; box=box1, inflate=1e-8)
for i in dbc_nodes1
    setebc!(T1, [i], 1, 0.0)
end

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, TriRule(1)), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]
F1 = zeros(size(K1, 1))
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]

edge_nodes1 = selectnode(fens1; box=[width1,width1, 0.0,height1], inflate=1e-8)
boundaryfes1 = meshboundary(fes1)
edge_fes1 = selectelem(fens1, boundaryfes1, withnodes=edge_nodes1)
#########################################################################################
width2 = 1.0
height2 = 1.0
n2x = 1
n2y = 3
fens2, fes2 = T3block(width2, height2, n2x, n2y)
# shift the second mesh to the right by 1.0
fens2.xyz[:, 1] .+= 1.0

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field
numberdofs!(T2)
femm2 = FEMMHeatDiff(IntegDomain(fes2, TriRule(1)), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]

edge_nodes2 = selectnode(fens2; box=[1.0,1.0, 0.0,height2], inflate=1e-8)
boundaryfes2 = meshboundary(fes2)
edge_fes2 = selectelem(fens2, boundaryfes2, withnodes=edge_nodes2)
##########################################################################################

xs_i = [1.0,1.0,1.0]
ys_i = [0.0,1.0, 2.0]
fens_i, fes_i = L2blockx2D(xs_i, ys_i)
geom_i = NodalField(fens_i.xyz)
u_i = NodalField(zeros(size(fens_i.xyz, 1), 2)) # displacement field
numberdofs!(u_i)
D1 = build_D_matrix(fens_i, fes_i, fens1, edge_fes1,
                     boundaryfes1; lam_order=1,tol=1e-8)
aaa = 0