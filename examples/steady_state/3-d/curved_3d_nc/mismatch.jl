using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
# include("utilities_old.jl")
include("meshrefine.jl")

m=1
N_elem1 = 7*m
N_elem2 = 6*m
N_elem_i = 6*m
left_m = "t"
right_m = "t"
skew = 0.
lam_order = 1
bend=0.0
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

fens1.xyz[:, 1] .+= bend * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5).^2#

geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, Rule1), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]

l1 = selectelem(fens1, meshboundary(fes1), box = [0.0,1.0, 0.0,0.0, 0.0, 1.0], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), TriRule(3)))
fi1 = ForceIntensity(Float64[1.0])
F1 = distribloads(el1femm, geom1, T1, fi1, 2)

l1 = selectelem(fens1, meshboundary(fes1), box = [0.0,1.0, 1.0,1.0, 0.0, 1.0], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), TriRule(3)))
fi1 = ForceIntensity(Float64[-1.0])
F1 += distribloads(el1femm, geom1, T1, fi1, 2)
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]
##########################################################################################
width2 = 0.5
height2 = 0.5
depth2 = 0.5
if right_m == "h"
    fens2, fes2 = H8block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = GaussRule(3,2)
else
    fens2, fes2 = T4block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = TetRule(4)
    
end 
fens2.xyz[:,1] .+= 0.5
fens2.xyz[:,2] .+=0.25
fens2.xyz[:,3] .+=0.25

boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2,  box=[0.5,0.5, 0.25,0.75, 0.25, 0.75], inflate=1e-8))

fens2.xyz[:, 1] .+= bend * (1.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 0.5).^2

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field

box2 = [1.0,1.0,0.25,0.25,0.25,0.25]
dbc_nodes2 = selectnode(fens2; box=box2, inflate=1e-8)
for i in dbc_nodes2
    setebc!(T2, [i], 1, 0.0)
end

applyebc!(T2)
numberdofs!(T2)

femm2 = FEMMHeatDiff(IntegDomain(fes2, Rule2), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]

l2 = selectelem(fens2, meshboundary(fes2), box = [0.0,1.0, 0.25,0.25, 0.0, 1], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), TriRule(3)))
fi2 = ForceIntensity(Float64[1.0])
F2 = distribloads(el2femm, geom2, T2, fi2, 2)

l2 = selectelem(fens2, meshboundary(fes2), box = [0.0,1.0, 0.75,0.75, 0.0, 1], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), TriRule(3)))
fi2 = ForceIntensity(Float64[-1.0])
F2+= distribloads(el2femm, geom2, T2, fi2, 2)

# F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]
##########################################################################################

xs_i = 0.5
ys_i = collect(linearspace(0.25, 0.75, N_elem_i+1))
zs_i = collect(linearspace(0.25, 0.75, N_elem_i+1))
fens_i, fes_i = T3blockx(ys_i, zs_i, :a)
fens_i.xyz = hcat(xs_i*ones(size(fens_i.xyz, 1), 1), fens_i.xyz)
fens_i.xyz[:, 1] .+= bend * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5).^2

if lam_order==1
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1))
else
    u_i  = ElementalField(zeros(size(fes_i.conn, 1), 1))
end

femm_i = FEMMHeatDiff(IntegDomain(fes_i, TriRule(9)), material)
geom_i = NodalField(fens_i.xyz)
applyebc!(u_i)

numberdofs!(u_i)

# D1, Pi_NC1, Pi_phi1, M_u1 = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
# D2, Pi_NC2, Pi_phi2, M_u2 = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)


# @time D1, meta1 = common_refinement(fens1, edge_fes1, fens_i, fes_i; lam_order=lam_order, h=0.03)
# @time D2, meta2 = common_refinement(fens2, edge_fes2, fens_i, fes_i; lam_order=lam_order, h=0.033)

@time D1, _ = common_refinement(fens1, edge_fes1, fens_i, fes_i; lam_order=lam_order, h=1.2/N_elem1)
@time D2, _ = common_refinement(fens2, edge_fes2, fens_i, fes_i; lam_order=lam_order, h=1.2/N_elem2)

# error()
# 
D2 = D2[:, setdiff(1:count(fens2), dbc_nodes2)]


A = [K1_ff          spzeros(size(K1_ff,1), size(K2_ff,2))    D1';
     spzeros(size(K2_ff,1), size(K1_ff,2))     K2_ff          -D2';
     D1               -D2               spzeros(size(D1,1), size(D1,1))]
B = vcat(F1_ff, F2_ff, zeros(size(D1,1)))
X = A \ B

scattersysvec!(T1, X[1:size(K1_ff,1)])
scattersysvec!(T2, X[size(K1_ff,1)+1 : size(K1_ff,1)+size(K2_ff,1)])
scattersysvec!(u_i, X[size(K1_ff,1)+size(K2_ff,1)+1 : end])

sol(x,y) = 0.25-y
err1 = L2_err(femm1, geom1, T1, sol)
err2 = L2_err(femm2, geom2, T2, sol)

File1 = "mismatch_left.vtk"
vtkexportmesh(
    File1,
    fens1, fes1,scalars = [("Temperature", T1.values), ("Err", err1.values)]
)
File2 = "mismatch_right.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values), ("Err", err2.values)]
)

# u1_fens = meta1["fens_u"]
# u2_fens = meta2["fens_u"]
# u1_fes = meta1["fes_u"]
# u2_fes = meta2["fes_u"]
# file3 = "mm_union_left.vtk"
# vtkexportmesh(
#     file3,
#     u1_fens, u1_fes,scalars = []
# )
# file4 = "mm_union_right.vtk"
# vtkexportmesh(
#     file4,
#     u2_fens, u2_fes,scalars = []
# )
println(u_i.values)

