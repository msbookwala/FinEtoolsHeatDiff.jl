using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")


N_elem1 = 2
N_elem2 = 3
N_elem_i = min(N_elem1, N_elem2)
left_m = "t"
right_m = "t"
skew = 0.
lam_order = 0

kappa = [1.0 0.0 0.0; 0 1.0 0.0; 0.0 0.0 1.0] 
material = MatHeatDiff(kappa)

#########################################################################################
width1 = 0.5
height1 = 1.0
depth1 = 1.0
if left_m == "h"
    fens1, fes1 = H8block(width1, height1, depth1, floor(Int, N_elem1), N_elem1, N_elem1)
    Rule1 = GaussRule(3,2)
else
    fens1, fes1 = T4block(width1, height1, depth1, floor(Int, N_elem1), N_elem1, N_elem1)
    Rule1 = TetRule(4)
end

boundaryfes1 = meshboundary(fes1)
edge_fes1 = subset(boundaryfes1, selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1, 0.0, depth1], inflate=1e-8))

# fens1.xyz[:, 1] .+= skew * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5)
geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, Rule1), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]

l1 = selectelem(fens1, meshboundary(fes1), box = [0.0,0.0, 0.0,height1, 0.0, depth1], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), TriRule(3)))
fi1 = ForceIntensity(Float64[-1.0])
F1 = distribloads(el1femm, geom1, T1, fi1, 2)
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]
##########################################################################################
width2 = 0.5
height2 = 1.0
depth2 = 1.0
if right_m == "h"
    fens2, fes2 = H8block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = GaussRule(3,2)
else
    fens2, fes2 = T4block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = TetRule(4)
    fens2.xyz[:,1] .+= 0.5
end 
boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2,  box=[0.5,0.5, 0.0,height2, 0.0, depth2], inflate=1e-8))

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
l2 = selectelem(fens2, meshboundary(fes2), box = [1.0,1.0, 0.0,height2, 0.0, depth2], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), TriRule(3)))
fi2 = ForceIntensity(Float64[1.0])
F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]
##########################################################################################

xs_i = 0.5
ys_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
zs_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
fens_i, fes_i = T3blockx(ys_i, zs_i, :b)
fens_i.xyz = hcat(xs_i*ones(size(fens_i.xyz, 1), 1), fens_i.xyz)

u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1))
femm_i = FEMMHeatDiff(IntegDomain(fes_i, TriRule(9)), material)
geom_i = NodalField(fens_i.xyz)
applyebc!(u_i)

numberdofs!(u_i)

M_i =mass(femm_i, geom_i, u_i)
D1, Pi_NC1, Pi_phi1, M_u1 = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
D2, Pi_NC2, Pi_phi2, M_u2 = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)

D2 = D2[:, setdiff(1:count(fens2), dbc_nodes2)]

A = [K1_ff          zeros(size(K1_ff,1), size(K2_ff,2))    D1';
     zeros(size(K2_ff,1), size(K1_ff,2))     K2_ff          -D2';
     D1               -D2               zeros(size(D1,1), size(D1,1))]
B = vcat(F1_ff, F2_ff, zeros(size(D1,1)))
X = A \ B

scattersysvec!(T1, X[1:size(K1_ff,1)])
scattersysvec!(T2, X[size(K1_ff,1)+1 : size(K1_ff,1)+size(K2_ff,1)])
scattersysvec!(u_i, X[size(K1_ff,1)+size(K2_ff,1)+1 : end])

sol(x,y) = x-1
err1 = L2_err(femm1, geom1, T1, sol)
err2 = L2_err(femm2, geom2, T2, sol)

File1 = "patch_test_left.vtk"
vtkexportmesh(
    File1,
    fens1, fes1,scalars = [("Temperature", T1.values), ("Err", err1.values)]
)
File2 = "patch_test_right.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values), ("Err", err2.values)]
)
# println(u_i.values)

