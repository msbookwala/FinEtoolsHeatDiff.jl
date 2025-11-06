using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")


N_elem1 = 7
N_elem2 = 11
N_elem_i = min(N_elem1, N_elem2)
left_m = "t"
right_m = "q"
skew = 0.0
bend = 0.2
lam_order = 1

kappa = [1.0 0; 0 1.0] 
material = MatHeatDiff(kappa)

#########################################################################################
width1 = 1.0
height1 = 2.0
if left_m == "t"
    fens1, fes1 = T6block(width1, height1, floor(Int, N_elem1/2), N_elem1)
    Rule1 = TriRule(9)
else
    xs1 = collect(linearspace(0.0, width1, floor(Int, N_elem1/2)+1))
    ys1 = collect(linearspace(0.0, height1, N_elem1+1))
    fens1, fes1 = Q9blockx(xs1, ys1)
    Rule1 = GaussRule(2,4)
end

boundaryfes1 = meshboundary(fes1)
edge_fes1 = subset(boundaryfes1, selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1], inflate=1e-8))

p = maximum(length.(edge_fes1.conn)) - 1
#########################################################################################


width2 = 1.0
height2 = 2.0
if right_m == "t"
    fens2, fes2 = T6block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    Rule2 = TriRule(3)
else
    xs2 = collect(linearspace(1.0 ,1.0+ width2, floor(Int, N_elem2/2)+1))
    ys2 = collect(linearspace(0.0, height2, N_elem2+1))
    fens2, fes2 = Q9blockx(xs2, ys2)
    Rule2 = GaussRule(2,4)
end
# shift the second mesh to the right by 1.0
# fens2.xyz[:, 1] .+= 1.0

boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2, box=[1.0,1.0, 0.0,height2], inflate=1e-8))

##########################################################################################

xs_i = ones(N_elem_i+1)
ys_i = collect(linearspace(0.0, 2.0, N_elem_i+1))
fens_i, fes_i = L3blockx2D(xs_i, ys_i)

fens_u1, fes_u1, M_u1 = build_union_mesh(fens_i,fes_i, fens1, edge_fes1, p; lam_order=lam_order)
fens_u2, fes_u2, M_u2 = build_union_mesh(fens_i,fes_i, fens2, edge_fes2, p; lam_order=lam_order)
############################################################################################

# fens1.xyz[:, 1] .+= skew * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 1.0)
fens1.xyz[:, 1] .+= bend * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 1.0).^2
#

geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

# box1 = [0.0,0.0,0.0,0.0]
# dbc_nodes1 = selectnode(fens1; box=box1, inflate=1e-8)
# for i in dbc_nodes1
#     setebc!(T1, [i], 1, 0.0)
# end

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, Rule1), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]

l1 = selectelem(fens1, meshboundary(fes1), box = [0.0,0.0, 0.0,height1], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), GaussRule(1,2)))
fi1 = ForceIntensity(Float64[-1.0])
F1 = distribloads(el1femm, geom1, T1, fi1, 2)
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]


################################################################################

# fens2.xyz[:, 1] .+= skew * (2.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 1.0)
fens2.xyz[:, 1] .+= bend * (2.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 1.0).^2

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field

box2 = [2.0,2.0,0.0,0.0]
dbc_nodes2 = selectnode(fens2; box=box2, inflate=1e-8)
for i in dbc_nodes2
    setebc!(T2, [i], 1, 0.0)
end

applyebc!(T2)

numberdofs!(T2)
femm2 = FEMMHeatDiff(IntegDomain(fes2, Rule2), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]

l2 = selectelem(fens2, meshboundary(fes2), box = [2.0,2.0, 0.0,height2], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), GaussRule(1,2)))
fi2 = ForceIntensity(Float64[1.0])
F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]


#################################################################################
# fens_i.xyz[:, 1] .+= skew * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 1.0)
fens_i.xyz[:, 1] .+= bend * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 1.0).^2
fens_u1.xyz[:, 1] .+= bend * fens_u1.xyz[:, 1].*(fens_u1.xyz[:, 2] .- 1.0).^2
fens_u2.xyz[:, 1] .+= bend * fens_u2.xyz[:, 1].*(fens_u2.xyz[:, 2] .- 1.0).^2

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 1)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1)) # Lagrange multipliers field
end
numberdofs!(u_i)
femm_i = FEMMHeatDiff(IntegDomain(fes_i, GaussRule(1,2)), material)


D1,Pi_NC1,Pi_phi1 = build_D_matrix(fens_u1, fes_u1, M_u1, fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
D2,Pi_NC2,Pi_phi2 = build_D_matrix(fens_u2, fes_u2, M_u2, fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)
# (error("oh no!"))


# D1 = D1[:, setdiff(1:count(fens1), dbc_nodes1)]
D2 = D2[:, setdiff(1:count(fens2), dbc_nodes2)]

A = [K1_ff          zeros(size(K1_ff,1), size(K2_ff,2))    D1';
     zeros(size(K2_ff,1), size(K1_ff,2))     K2_ff          -D2';
     D1               -D2               zeros(size(D1,1), size(D1,1))]
# A = cholesky(A)
# show(Matrix(A))
B = vcat(F1_ff, F2_ff, zeros(size(D1,1)))
X = A \ B

scattersysvec!(T1, X[1:size(K1_ff,1)])
scattersysvec!(T2, X[size(K1_ff,1)+1 : size(K1_ff,1)+size(K2_ff,1)])
scattersysvec!(u_i, X[size(K1_ff,1)+size(K2_ff,1)+1 : end])

z1 = zeros(size(fens1.xyz, 1))
z2 = zeros(size(fens2.xyz, 1))
# geom1vals = hcat(geom1.values, z1)
# geom2vals = hcat(geom2.values, z2)
geom1vals = geom1.values
geom2vals = geom2.values

File1 = "patch_test_left.vtk"
vtkexportmesh(
    File1,
    fens1, fes1,scalars = [("Temperature", T1.values)]
)
File2 = "patch_test_right.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values)]
)
println(u_i.values)
