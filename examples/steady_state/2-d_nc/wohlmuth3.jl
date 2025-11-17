using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
include("utilities.jl")

lam_order = 0

function f(x,y)
    return (x-0.5)*(y-0.5)*exp(-10*((x-0.5)^2 + (y-0.5)^2))/a(x,y)
end
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
    forceout[1] = (20*exp(- 10*(x - 1/2)^2 - 10*(y - 1/2)^2)*(2*x - 1)*(2*y - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 1))/a(x,y)
    return forceout
end
boundary_boxes = [
    [0.0,0.0, 0.0,1.0], # left
    [1.0,1.0, 0.0,1.0], # right
    [0.0,1.0, 1.0,1.0], # top
    [0.0,1.0, 0.0,0.0], # bottom
] 

ctr = 1
nelems = [4,1,1,4]
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

        kappa = [a(0.25+0.5*(i-1),0.25+0.5*(j-1)) 0; 0 a(0.25+0.5*(i-1),0.25+0.5*(j-1))] 
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
edge_fes_all = extract_interface_fes(bfes_all, fens_all, boundary_boxes)

N_elem_i = min(nelems[1], nelems[4])

xs_i1 = 0.5*ones(N_elem_i+1)
ys_i1 = collect(linearspace(0.0, 1.0, N_elem_i+1))
fensi_1, fesi_1 = L2blockx2D(xs_i1, ys_i1)


xs_i2 = collect(linearspace(0.0, 1.0, N_elem_i+1))
ys_i2 = 0.5*ones(N_elem_i+1)
fensi_2, fesi_2 = L2blockx2D(xs_i2, ys_i2)
fens_i = [fensi_1; fensi_2]
fes_i = [fesi_1; fesi_2]

fens_us = []
F_lam =[]
D = []

for i in 1:4
    D_ = []
    F_lam_ = []
    for j in 1:2

        fens_u, fes_u, _ = build_union_mesh(fens_i[j],fes_i[j], fens_all[i], edge_fes_all[i], 1; lam_order=lam_order, to_trim =true, dir = 3-j)
        Dij, _,_ = build_D_matrix(fens_u, fes_u, fens_i[j], fes_i[j], fens_all[i], edge_fes_all[i]; lam_order=lam_order)
        F_lam_ij = Dij

        push!(D_, Dij)
        push!(fens_us, fens_u)
    end
    push!(D, D_)
end

