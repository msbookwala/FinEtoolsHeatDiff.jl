using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")

lam_order = 0

function f(x,y)
    return 1+x.^2 + 2*y.^2
end
function lam1(x,y)
    return 1
end
function lam2(x,y)
    return 2
end
exact_lam= [lam1, lam2]
function a(x,y)
    if x<0.5 && y<0.5
        return 0.000025
    elseif x>=0.5 && y<0.5
        return 1.0
    elseif x<0.5 && y>=0.5
        return 1.0
    elseif x>=0.5 && y>=0.5
        return 0.00025
    end
end
function Q!(forceout, XYZ, tangents, feid, qpid)
    x = XYZ[1]
    y = XYZ[2]
    forceout[1] = -6.0
    return forceout
end
boundary_boxes = [
    [0.0,0.0, 0.0,1.0], # left
    [1.0,1.0, 0.0,1.0], # right
    [0.0,1.0, 1.0,1.0], # top
    [0.0,1.0, 0.0,0.0], # bottom
] 

ctr = 1
# r = 5
nelems = [4,2,2,1]*(2^r)
fes_all = []
fens_all = []
K_all =[]
F_all =[]
femm_all = []
T_all = []
geom_all = []
dbc_nodes_all = []
material_all = []
bfes_all = []
for i in 1:2
    for j in 1:2
        n = nelems[2*(i-1)+j]
        fens_local,fes_local = T3block(0.5,0.5,n,n)
        fens_local.xyz[:,1] .+= (i-1)*0.5
        fens_local.xyz[:,2] .+= (j-1)*0.5
        bfes_local = meshboundary(fes_local)

        geom_local = NodalField(fens_local.xyz)
        T_local = NodalField(zeros(size(fens_local.xyz, 1), 1)) # displacement field

        dbc_nodes_local = []
        for box in boundary_boxes
            nodes_in_box = selectnode(fens_local; box=box, inflate=1e-8)
            append!(dbc_nodes_local, nodes_in_box)
        end

        sort!(unique!(dbc_nodes_local))
        for i_ in dbc_nodes_local
            setebc!(T_local, [i_], 1, f(fens_local.xyz[i_,1], fens_local.xyz[i_,2]))
        end
        applyebc!(T_local)
        numberdofs!(T_local)

        kappa = [1.0 0; 0 1.0] 
        material_local = MatHeatDiff(kappa)

        femm_local = FEMMHeatDiff(IntegDomain(fes_local, TriRule(3)), material_local)
        K_local = conductivity(femm_local, geom_local, T_local)
        K_ff_local = matrix_blocked(K_local, nfreedofs(T_local), nfreedofs(T_local))[:ff]
        K_fd_local = matrix_blocked(K_local, nfreedofs(T_local), nfreedofs(T_local))[:fd]

        fis_local = ForceIntensity(Float64,1, Q!)
        F_local = distribloads(femm_local, geom_local, T_local, fis_local, 2)
        F_ff_local = vector_blocked(F_local, nfreedofs(T_local))[:f] - K_fd_local * gathersysvec(T_local, :d)

        push!(fes_all, fes_local)
        push!(fens_all, fens_local)
        push!(K_all, K_ff_local)
        push!(F_all, F_ff_local)
        push!(femm_all, femm_local)
        push!(T_all, T_local)
        push!(geom_all, geom_local)
        push!(dbc_nodes_all, dbc_nodes_local)
        push!(material_all, material_local)
        push!(bfes_all, bfes_local)
    end
end
#assemble Ks, push K_all[i] at diagonal
K = []
for i in 1:4
    row = []
    for j in 1:4
        if i == j
            push!(row, K_all[i])
        else
            push!(row, zeros(Float64, size(K_all[i],1), size(K_all[j],2)))

        end
    end
    row = hcat(row...)
    push!(K, row)
end
K = vcat(K...)
F = vcat(F_all...)

edge_fes_all = extract_interface_fes(bfes_all, fens_all, boundary_boxes)

N_elem_i = 2*nelems[4]

xs_i1 = 0.5*ones(N_elem_i+1)
ys_i1 = collect(linearspace(0.0, 1.0, N_elem_i+1))
fensi_1, fesi_1 = L2blockx2D(xs_i1, ys_i1)


xs_i2 = collect(linearspace(0.0, 1.0, N_elem_i+1))
ys_i2 = 0.5*ones(N_elem_i+1)
fensi_2, fesi_2 = L2blockx2D(xs_i2, ys_i2)
fens_i = [fensi_1; fensi_2]
fes_i = [fesi_1; fesi_2]
u_i = []
F_lam = []
geom_i = []
femm_i = []
for j in 1:2
    if lam_order == 0
        push!(u_i, ElementalField(zeros(count(fes_i[j]), 1))) # Lagrange multipliers field
    else
        push!(u_i, NodalField(zeros(size(fens_i[j].xyz, 1), 1))) # Lagrange multipliers field
    end
    numberdofs!(u_i[j])
    femm = FEMMHeatDiff(IntegDomain(fes_i[j], GaussRule(1,4)), MatHeatDiff(reshape([1.0], 1, 1)))
    geom = NodalField(fens_i[j].xyz)
    push!(femm_i, femm)
    push!(geom_i, geom)
    push!(F_lam, zeros(nfreedofs(u_i[j])))
end

fens_us = []

D = []
multipliers = [-1 -1; -1 1; 1 -1; 1 1]
for i in 1:4
    D_ = []
    for j in 1:2
        fens_u, fes_u, _ = build_union_mesh(fens_i[j],fes_i[j], fens_all[i], edge_fes_all[i], 1; lam_order=lam_order, to_trim =true, dir = 3-j)
        Dij, _,_ = build_D_matrix(fens_u, fes_u, fens_i[j], fes_i[j], fens_all[i], edge_fes_all[i]; lam_order=lam_order)
        
        Dij = multipliers[i,j] * Dij
        F_lam[j] += - Dij[:, dbc_nodes_all[i]] * gathersysvec(T_all[i], :d)

        Dij = Dij[:, setdiff(1:count(fens_all[i]), dbc_nodes_all[i])]
        push!(D_, Dij)
        push!(fens_us, fens_u)
    end
    D_ = vcat(D_...)
    push!(D, D_)
end
D_mat = hcat(D...)

A = [K D_mat'; D_mat zeros(Float64, size(D_mat,1), size(D_mat,1))]
b = vcat(F, vcat(F_lam...))
x = A\b

# distibute solution back to fields
offset = 0
for i in 1:4
    ndofs_j = size(K_all[i],1)
    scattersysvec!(T_all[i], x[offset+1:offset+ndofs_j])
    global offset += ndofs_j
end
for j in 1:2
    ndofs_j = nents(u_i[j])
    scattersysvec!(u_i[j], x[offset+1:offset+ndofs_j])
    # println("Lagrange multipliers field ", j, ":", u_i[j].values)
    global offset += ndofs_j
end
errs = []
T_l2sq = 0.0
for i in 1:4
    global T_l2sq
    err = L2_err(femm_all[i], geom_all[i], T_all[i], f)
    push!(errs, err)
    T_l2sq += sum(err.values.^2)
end
T_l2 = sqrt(T_l2sq)

l_l2sq = 0.0
for j in 1:2
    global l_l2sq
    err =L2_err(femm_i[j], geom_i[j], u_i[j], exact_lam[j])
    l_l2sq += sum(err.values.^2)
end
l_l2 = sqrt(l_l2sq)
println("T_l2 = ", T_l2)
println("l_l2 = ", l_l2)

for i in 1:4
    filename = "quad_4_$i.vtk"
    vtkexportmesh(
        filename,
        fens_all[i],
        fes_all[i];
        scalars = [
            ("Temperature", T_all[i].values),
            ("Err", errs[i].values)
        ],
    )
end
