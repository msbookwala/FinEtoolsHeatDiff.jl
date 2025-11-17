using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
using Plots
include("utilities.jl")

# r = 0
N_elem1 = 2*(2^r)
N_elem2 = 3*(2^r)
N_elem_i = min(N_elem1, N_elem2)
left_m = "q"
right_m = "q"
skew = 0.
bend = 0.5
# lam_order = 2
kappa = [1.0 0; 0 1.0] 
material = MatHeatDiff(kappa)
Q = -6.0
sol(x,y) = 1 .+ x.^2 + 2*y.^2
q(x,y) = 4*y

function topflux!(forceout, XYZ, tangents, feid, qpid)
    forceout[1] = q(XYZ[1], XYZ[2]) #heat source
    return forceout
end
function botflux!(forceout, XYZ, tangents, feid, qpid)
    forceout[1] = -q(XYZ[1], XYZ[2]) #heat source
    return forceout
end

#########################################################################################
width1 = 0.5
height1 = 1.0
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



#########################################################################################


width2 = 0.5
height2 = 1.0
if right_m == "t"
    fens2, fes2 = T6block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    Rule2 = TriRule(9)
    fens2.xyz[:, 1] .+= 0.5
else
    xs2 = collect(linearspace(0.5 ,0.5+ width2, floor(Int, N_elem2/2)+1))
    ys2 = collect(linearspace(0.0, height2, N_elem2+1))
    fens2, fes2 = Q9blockx(xs2, ys2)
    Rule2 = GaussRule(2,4)
end


boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2, box=[0.5,0.5, 0.0,height2], inflate=1e-8))

##########################################################################################

p = maximum(length.(edge_fes1.conn)) - 1
p_sd = maximum(length.(edge_fes1.conn)) - 1
p_i = min(p_sd, lam_order)

xs_i = 0.5*ones(N_elem_i+1)
ys_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
if p_i < 2 && bend == 0
    fens_i, fes_i = L2blockx2D(xs_i, ys_i)
else
    fens_i, fes_i = L3blockx2D(xs_i, ys_i)
end

fens_u1, fes_u1, _ = build_union_mesh(fens_i,fes_i, fens1, edge_fes1, p; lam_order=lam_order)
fens_u2, fes_u2, _ = build_union_mesh(fens_i,fes_i, fens2, edge_fes2, p; lam_order=lam_order)
############################################################################################



fens1.xyz[:, 1] .+= skew * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5)
fens1.xyz[:, 1] .+= bend * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5).^2
#

geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

boxleft = [0.0,0.0,0.0,1.0]
dbc_nodes1 = selectnode(fens1; box=boxleft, inflate=1e-8)
for i in dbc_nodes1
    setebc!(T1, [i], 1, sol(fens1.xyz[i,1], fens1.xyz[i,2]))
end

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, Rule1), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]
K1_fd = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:fd]

l1top = selectelem(fens1, meshboundary(fes1), box = [0.0,1., height1,height1], inflate=1e-8)
l1bot = selectelem(fens1, meshboundary(fes1), box = [0.0,1., 0.0,0.0], inflate=1e-8)
el1femmbot = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1bot), GaussRule(1,2)))
el1femmtop = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1top), GaussRule(1,2)))

fi1 = ForceIntensity(Float64, 1, topflux!)
fi2 = ForceIntensity(Float64, 1, botflux!)
fis = ForceIntensity(Float64[Q])
F1 = distribloads(el1femmtop, geom1, T1, fi1, 2)
F1 += distribloads(el1femmbot, geom1, T1, fi2, 2)
F1 += distribloads(femm1, geom1, T1, fis, 3) 
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f] - K1_fd * gathersysvec(T1, :d)



################################################################################

fens2.xyz[:, 1] .+= skew * (1.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 0.5)
fens2.xyz[:, 1] .+= bend * (1.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 0.5).^2

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field

box2 = [1.0,1.0,0.0,1.0]
dbc_nodes2 = selectnode(fens2; box=box2, inflate=1e-8)
for i in dbc_nodes2
    setebc!(T2, [i], 1, sol(fens2.xyz[i,1], fens2.xyz[i,2]))
end

applyebc!(T2)

numberdofs!(T2)
femm2 = FEMMHeatDiff(IntegDomain(fes2, Rule2), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]
K2_fd = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:fd]


l2top = selectelem(fens2, meshboundary(fes2), box = [0.,1.0, height2,height2], inflate=1e-8)
l2bot = selectelem(fens2, meshboundary(fes2), box = [0.,1.0, 0.0,0.0], inflate=1e-8)
el2femmtop = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2top), GaussRule(1,2)))
el2femmbot = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2bot), GaussRule(1,2)))


F2 = distribloads(el2femmtop, geom2, T2, fi1, 2)
F2 += distribloads(el2femmbot, geom2, T2, fi2, 2)
F2 += distribloads(femm2, geom2, T2, fis, 3)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f] - K2_fd * gathersysvec(T2, :d)


#################################################################################
fens_i.xyz[:, 1] .+= skew * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5)
fens_u1.xyz[:, 1] .+= skew * fens_u1.xyz[:, 1].*(fens_u1.xyz[:, 2] .- 0.5)
fens_u2.xyz[:, 1] .+= skew * fens_u2.xyz[:, 1].*(fens_u2.xyz[:, 2] .- 0.5)

fens_i.xyz[:, 1] .+= bend * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5).^2
fens_u1.xyz[:, 1] .+= bend * fens_u1.xyz[:, 1].*(fens_u1.xyz[:, 2] .- 0.5).^2
fens_u2.xyz[:, 1] .+= bend * fens_u2.xyz[:, 1].*(fens_u2.xyz[:, 2] .- 0.5).^2

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 1)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1)) # Lagrange multipliers field
end
numberdofs!(u_i)
femm_i = FEMMHeatDiff(IntegDomain(fes_i, GaussRule(1,2)), MatHeatDiff(reshape([1.0], 1, 1)))


D1,Pi_NC1,Pi_phi1 = build_D_matrix(fens_u1, fes_u1, fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
D2,Pi_NC2,Pi_phi2 = build_D_matrix(fens_u2, fes_u2, fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)
# (error("oh no!"))

f_lam = -D1[:, dbc_nodes1] * gathersysvec(T1, :d) + D2[:, dbc_nodes2] * gathersysvec(T2, :d)
D1 = D1[:, setdiff(1:count(fens1), dbc_nodes1)]
D2 = D2[:, setdiff(1:count(fens2), dbc_nodes2)]

A = [K1_ff          zeros(size(K1_ff,1), size(K2_ff,2))    D1';
     zeros(size(K2_ff,1), size(K1_ff,2))     K2_ff          -D2';
     D1               -D2               zeros(size(D1,1), size(D1,1))]
# A = cholesky(A)
# show(Matrix(A))
B = vcat(F1_ff, F2_ff, f_lam)
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

l2err1 = L2_err(femm1, geom1, T1, sol)
l2err2 = L2_err(femm2, geom2, T2, sol)


File1 = "quadratic_test_left.vtk"
vtkexportmesh(
    File1,
    fens1, fes1,scalars = [("Temperature", T1.values), ("Err", l2err1.values)]
)
File2 = "quadratic_test_right.vtk"
vtkexportmesh(
    File2,
    fens2, fes2,scalars = [("Temperature", T2.values), ("Err", l2err2.values)]
)
println(u_i.values)
# # plot(geom_i.values[:,1], u_i.values, seriestype=:scatter, title="Lagrange Multipliers", xlabel="Node Number", ylabel="Multiplier Value")
# plot( u_i.values, seriestype=:scatter, title="Lagrange Multipliers", xlabel="Node Number", ylabel="Multiplier Value")

tot_l2 = sqrt(sum(l2err1.values.^2) + sum(l2err2.values.^2))

# exact_lagrange(x,y) = -2*x*cos(atan(0.5*skew)) + 4*y*sin(atan(0.5*skew))
exact_lagrange(x,y) = -2*x/sqrt(1+(bend*(y-0.5))^2) + 4*y*(bend*(y-0.5))/sqrt(1+(bend*(y-0.5))^2)


lag_err = L2_err(femm_i, geom_i, u_i, exact_lagrange)
tot_lag_err = sqrt(sum(lag_err.values.^2))