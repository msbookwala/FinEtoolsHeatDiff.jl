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

#########################################################################################
width1 = 0.5
height1 = 1.0
depth1 = 1.0
if left_m == "h"
    fens1, fes1 = H8block(width1, height1, depth1, 1, 1, 2)
    Rule1 = GaussRule(3,2)
else
    fens1, fes1 = T4block(width1, height1, depth1, floor(Int, N_elem1), N_elem1, N_elem1)
    Rule1 = TetRule(4)
end
fens1.xyz[6,3] +=0.25
fens1.xyz[8,3] -=0.25
# error("here")
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
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), GaussRule(2,2)))
fi1 = ForceIntensity(Float64[-1.0])
F1 = distribloads(el1femm, geom1, T1, fi1, 2)
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]

l1_ = selectelem(fens1, meshboundary(fes1), box = [0.5,0.5, 0.0,height1, 0.0, depth1], inflate=1e-8)
el1femm_ = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1_), GaussRule(2,2)))
fi1_ = ForceIntensity(Float64[1.0])
F1_ = distribloads(el1femm_, geom1, T1, fi1_, 2)
##########################################################################################
width2 = 0.5
height2 = 1.0
depth2 = 1.0
if right_m == "h"
    fens2, fes2 = H8block(width2, height2, depth2, 1, 1, 2)
    Rule2 = GaussRule(3,2)
else
    fens2, fes2 = T4block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = TetRule(4)
    
end 

fens2.xyz[:,1] .+= 0.5
fens2.xyz[5,3] -=0.25
fens2.xyz[7,3] +=0.25

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
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), GaussRule(2,2)))
fi2 = ForceIntensity(Float64[1.0])
F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]

l2_ = selectelem(fens2, meshboundary(fes2), box = [0.5,0.5, 0.0,height2, 0.0, depth2], inflate=1e-8)
el2femm_ = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2_), GaussRule(2,2)))
fi2_ = ForceIntensity(Float64[1.0])
F2_ = distribloads(el2femm_, geom2, T2, fi2_, 2)


# File1 = "left.vtk"
# vtkexportmesh(
#     File1,
#     fens1, fes1,scalars = []
# )
# File2 = "right.vtk"
# vtkexportmesh(
#     File2,
#     fens2, fes2,scalars = []
# )
# error("here")
##########################################################################################
fensixyz = [0.5 0.0 0.0;
             0.5 0.0 0.25;
             0.5 0.0 1.0 ;
             0.5 1.0 0.0;
             0.5 1.0 0.75;
             0.5 1.0 1.0]

fesiconn = [1 2 5 4;
                2 3 6 5]

fes_i =FESetQ4(fesiconn)
fens_i = FENodeSet(fensixyz)

if lam_order==1
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1))
else
    u_i = ElementalField(zeros(count(fes_i), 1))
end
numberdofs!(u_i)



D1, meta1 = common_refinement(fens1, edge_fes1, fens_i, fes_i; lam_order=lam_order, h=1.0, tri_order=2, triangulation_type = "naive" )
D2, meta2 = common_refinement(fens2, edge_fes2, fens_i, fes_i; lam_order=lam_order, h=1.0, tri_order=2, triangulation_type = "naive" )

# G1 = D1'*ones(size(D1,1))
# G2 = D2'*ones(size(D2,1))
# J1 = I(size(D1,2)) - G1*(F1_-G1)'/(G1'*G1)
# J2 = I(size(D2,2)) - G2*(F2_-G2)'/(G2'*G2)
# D1 = D1*J1
# D2 = D2*J2


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

File1 = "um_qqpatch_test_left.vtk"
vtkexportmesh(
    File1,
    fens1, fes1,scalars = [("Temperature", T1.values), ("Err", err1.values)]
)
File2 = "um_qqpatch_test_right.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values), ("Err", err2.values)]
)
file3 = "union.vtk"
vtkexportmesh(
    file3,
    meta1["fens_u"], meta1["fes_u"],scalars = []
)

println(u_i.values)

