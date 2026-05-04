using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("meshrefine.jl")


N_elem1 = 3
N_elem2 = 5
N_elem_i = min(N_elem1, N_elem2)
left_m = "h"
right_m = "h"
skew = 0.
lam_order = 0

kappa = [1.0 0.0 0.0; 0 1.0 0.0; 0.0 0.0 1.0] 
material = MatHeatDiff(kappa)


width2 = 1.0
height2 = 1.0
depth2 = 1.0
if right_m == "h"
    fens2, fes2 = H8block(width2, height2, depth2, 1, 1, 2)
    Rule2 = GaussRule(3,2)
else
    fens2, fes2 = T4block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = TetRule(4)
    
end 

fens2.xyz[5,3] -=0.25
fens2.xyz[7,3] +=0.25

boundaryfes2 = meshboundary(fes2)

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field

box2 = [1.0,1.0,0.0,0.0,0.0,0.0]
dbc_nodes2 = selectnode(fens2; box=box2, inflate=1e-8)
for i in dbc_nodes2
    setebc!(T2, [i], 1, 0.0)
end

applyebc!(T2)
numberdofs!(T2)

femm2 = FEMMHeatDiff(IntegDomain(fes2, Rule2), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]

l2 = selectelem(fens2, meshboundary(fes2), box = [0.0,1.0, 0.0,0.0, 0.0, 1.0], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), GaussRule(2,2)))


l1 = selectelem(fens2, meshboundary(fes2), box = [0.0,1.0, 1.0,1.0, 0.0, 1.0], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l1), GaussRule(2,2)))

fi2 = ForceIntensity(Float64[1.0])
fi1 = ForceIntensity(Float64[-1.0])

F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2 += distribloads(el1femm, geom2, T2, fi1, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]







X = K2_ff\ F2_ff

scattersysvec!(T2, X)

sol(x,y) = -y
err2 = L2_err(femm2, geom2, T2, sol)


File2 = "irq_test.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values), ("Err", err2.values)]
)

