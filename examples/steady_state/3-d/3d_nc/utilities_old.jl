using SparseArrays
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule
using Triangulate



function extract_interface_fes(edge_fe_s, fen_s, boxes)
    # returns indices of edge fes that are on the interface 
    # Extract boundary elements within given boxes and take the complement
    interface_fes = []
    for i in 1:length(edge_fe_s)
        boundary_fes = []
        for box in boxes
            boundary_fes = union(boundary_fes, selectelem(fen_s[i], edge_fe_s[i], box=box, inflate=1e-8))
        end
        interface_fes_idxi =  setdiff(range(1, count(edge_fe_s[i])), boundary_fes)
        push!(interface_fes, subset(edge_fe_s[i], interface_fes_idxi))
    end
    return interface_fes
end

# function make_union_mesh(sd_interface_fe_s, sd_fens_s,  fes_i, fens_i, p; lam_order=0)
#     # make a common list of points within the interface elements and then do unique
#     all_points = []
#     for (fe_s, fens_s) in zip(sd_interface_fe_s, sd_fens_s)
#         for fe in fe_s.conn
#             for node in fe
#                 push!(all_points, fens_s.xyz[node,:])
#             end
#         end
#     end
#     points = unique(all_points)
#     refine_1d_mesh(fens_i, fes_i, points; tol=1e-12)
#     a=1

# end

function build_D_matrix(fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8)

    nodes_i = unique(collect(Iterators.flatten(fes_i.conn[:])))
    nodes_sd = unique(collect(Iterators.flatten(edge_fes.conn[:])))

    pts_i  = fens_i.xyz[nodes_i, :]
    pts_sd = fens_sd.xyz[nodes_sd, :]
    all_points = vcat(pts_i, pts_sd)

    unique_points = unique(all_points, dims=1)
    unique_points_yz = unique_points[:, 2:3]
    tri_in = Triangulate.TriangulateIO()
    tri_in.pointlist = unique_points_yz' 
    tri_out, vor_out = Triangulate.triangulate("Q", tri_in)  



    # fens_u = FENodeSet(hcat(unique_points[:,1],tri_out.pointlist'))
    # fes_u = FESetT3(tri_out.trianglelist')

    fens_u = fens_i
    fes_u = fes_i

    File2 = "union_mesh.vtk"
    vtkexportmesh(
    File2,
    fens_u, fes_u,scalars = []
        )

    kappa = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] 
    material = MatHeatDiff(kappa)
    geom_u = NodalField(fens_u.xyz)
    u_u = NodalField(zeros(size(fens_u.xyz, 1), 1))
    numberdofs!(u_u)
    femm_u = FEMMHeatDiff(IntegDomain(fes_u, TriRule(9)), material)
    if lam_order == 0
        M_u = mass_like(femm_u, geom_u, u_u)
    else
        M_u = mass(femm_u, geom_u, u_u)
    end
    # print(length(fes_u.conn), " elements in the union mesh.\n")
    # print(size(M_u, 1), " DOFs in the union mesh.\n")

    p=1
    Pi_NC = Lagrange_interpolation_matrix(fens_u.xyz, fens_sd.xyz, edge_fes.conn, p; dim_u=1, tol=tol)
    Pi_phi = Lagrange_interpolation_matrix(fens_u.xyz, fens_i.xyz, fes_i.conn, p; dim_u=1, tol=tol)
    D = Pi_phi' * M_u * Pi_NC
    return D, Pi_NC, Pi_phi, M_u
end


function barycentric_coords_3D(P, A, B, C)
    v0 = B .- A
    v1 = C .- A
    nvec = cross(v0, v1)
    n_norm = norm(nvec)
    n_norm == 0.0 && error("barycentric_coords_3D: degenerate triangle.")
    n_unit = nvec ./ n_norm
    AP = P .- A
    dist = dot(AP, n_unit)
    Pproj = P .- dist .* n_unit
    APp = Pproj .- A

    d00 = dot(v0, v0)
    d01 = dot(v0, v1)
    d11 = dot(v1, v1)
    d20 = dot(APp, v0)
    d21 = dot(APp, v1)
    denom = d00*d11 - d01*d01
    abs(denom) < eps(Float64) && error("barycentric_coords_3D: triangle basis nearly singular.")

    v = (d11*d20 - d01*d21) / denom
    w = (d00*d21 - d01*d20) / denom
    u = 1.0 - v - w
    lambdas = [u, v, w]
    return lambdas, abs(dist)
end

function point_in_element(P, Y,
                          nodes, p::Int; tol)
    # 3D triangle; p is ignored for now (p=1 only)
    Pvec = Float64[P[1], P[2], P[3]]
    A = @view Y[nodes[1], :]
    B = @view Y[nodes[2], :]
    C = @view Y[nodes[3], :]
    lambdas, dist = barycentric_coords_3D(Pvec, A, B, C)
    inside = all(lambdas .>= -tol) && all(lambdas .<= 1.0 + tol)
    return inside, lambdas, dist
end



function Lagrange_interpolation_matrix(X, Y, conn, p::Int;
                                       dim_u::Int=1, tol::Float64=1e-8)
    p != 1 && error("Lagrange_interpolation_matrix (3D): only p=1 (linear triangles) implemented.")
    npts = size(X, 1)
    nnds = size(Y, 1)
    nels = size(conn, 1)

    I = Int[]
    J = Int[]
    V = Float64[]

    for r in 1:npts
        Pr = @view X[r, :]
        found = false
        for elem_idx in 1:nels
            nodes = conn[elem_idx][:]
            length(nodes) == 3 || error("Lagrange_interpolation_matrix: base elements must be 3-node triangles.")
            inside, lambdas, dist = point_in_element((Pr[1], Pr[2], Pr[3]),
                                                     Y, nodes, p; tol=tol)
            inside || continue
            for a in 1:3
                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[a]-1)*dim_u + k + 1)
                    push!(V, lambdas[a])
                end
            end
            found = true
            break
        end
        if !found
            error("Lagrange_interpolation_matrix: union node $r not found in base interface mesh.")
        end
    end

    return sparse(I, J, V, dim_u*npts, dim_u*nnds)
end