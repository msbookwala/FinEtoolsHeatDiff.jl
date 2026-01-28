using SparseArrays
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule
using Triangulate
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule

# ---------------------------
# Connectivity helpers
# ---------------------------
@inline _nnodes_conn(conn) = conn isa AbstractMatrix ? size(conn, 2) : length(conn[1])

@inline function _elem_nodes(conn, e::Int)
    if conn isa AbstractMatrix
        return collect(@view conn[e, :])
    else
        return collect(conn[e])  # tuple -> Vector
    end
end

@inline function _all_nodes(conn)
    if conn isa AbstractMatrix
        return unique(vec(conn))
    else
        return unique(collect(Iterators.flatten(conn)))
    end
end


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
# =========================================================================================
# 3D surface utilities: Q4 (quad) faces + rectangular fast path (no Newton)
# =========================================================================================

function _unique_sorted_tol(vals::AbstractVector{<:Real}; tol::Float64=1e-12)
    v = sort(Float64.(vals))
    out = Float64[]
    for x in v
        if isempty(out) || abs(x - out[end]) > tol
            push!(out, x)
        end
    end
    return out
end

@inline function q4_shapes(xi::Float64, eta::Float64)
    return 0.25 .* [
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ]
end

@inline function q4_shapes_deriv(xi::Float64, eta::Float64)
    # returns 4x2 matrix: [dN/dxi  dN/deta]
    dNdxi  = 0.25 .* [-(1-eta), +(1-eta), +(1+eta), -(1+eta)]
    dNdeta = 0.25 .* [-(1-xi),  -(1+xi), +(1+xi), +(1-xi)]
    return hcat(dNdxi, dNdeta)
end

function _is_parallelogram_quad(A::AbstractVector, B::AbstractVector, C::AbstractVector, D::AbstractVector; tol::Float64=1e-12)
    # Parallelogram condition: A + C == B + D
    L = max(norm(B-A), norm(C-B), norm(D-C), norm(A-D), 1.0)
    return norm(A .+ C .- B .- D) <= tol * L
end

function _plane_dist(P::AbstractVector, A::AbstractVector, B::AbstractVector, D::AbstractVector)
    n = cross(B .- A, D .- A)
    nn = norm(n)
    nn == 0.0 && return Inf
    nhat = n ./ nn
    return abs(dot(P .- A, nhat))
end

function invmap_q4_affine(P::AbstractVector, A::AbstractVector, B::AbstractVector, C::AbstractVector, D::AbstractVector)
    # For parallelogram quads, the standard bilinear map reduces to affine:
    x0 = (A .+ B .+ C .+ D) ./ 4.0
    a  = (-A .+ B .+ C .- D) ./ 4.0
    b  = (-A .- B .+ C .+ D) ./ 4.0
    M = hcat(a, b)  # 3x2
    uv = M \ (P .- x0)  # least squares (exact if in plane)
    return Float64(uv[1]), Float64(uv[2])
end

function invmap_q4_newton(P::AbstractVector, A::AbstractVector, B::AbstractVector, C::AbstractVector, D::AbstractVector;
                          xi0::Float64=0.0, eta0::Float64=0.0, tol::Float64=1e-13, maxiter::Int=50, alpha::Float64=1e-12)
    xi, eta = xi0, eta0
    X = (A, B, C, D)
    for _ in 1:maxiter
        N  = q4_shapes(xi, eta)
        dN = q4_shapes_deriv(xi, eta)

        x      = zeros(3)
        dx_dxi = zeros(3)
        dx_deta= zeros(3)
        @inbounds for a in 1:4
            x       .+= N[a]      .* X[a]
            dx_dxi  .+= dN[a, 1]  .* X[a]
            dx_deta .+= dN[a, 2]  .* X[a]
        end

        r = x .- P
        if norm(r) <= tol
            return xi, eta, true
        end

        J = hcat(dx_dxi, dx_deta)            # 3x2
        d = (J' * J .+ alpha * I(2)) \ (J' * (-r))  # 2x2 solve
        xi  += d[1]
        eta += d[2]
    end
    return xi, eta, false
end

function point_in_quad_Q4(P::AbstractVector, Y::AbstractMatrix, nodes::AbstractVector{<:Integer};
                          tol::Float64=1e-10, rectangles::Bool=true)
    length(nodes) == 4 || error("point_in_quad_Q4: need 4 nodes for Q4.")
    A = vec(@view Y[nodes[1], :])
    B = vec(@view Y[nodes[2], :])
    C = vec(@view Y[nodes[3], :])
    D = vec(@view Y[nodes[4], :])

    dist = _plane_dist(P, A, B, D)
    dist <= tol || return (false, zeros(4), dist)

    xi = 0.0; eta = 0.0
    if rectangles && _is_parallelogram_quad(A,B,C,D; tol=tol)
        xi, eta = invmap_q4_affine(P, A, B, C, D)
        ok = true
    else
        xi, eta, ok = invmap_q4_newton(P, A, B, C, D; tol=tol)
    end

    inside = ok && (abs(xi) <= 1.0 + 10*tol) && (abs(eta) <= 1.0 + 10*tol)
    N = q4_shapes(xi, eta)
    return inside, N, dist
end

# --- triangle helper (kept separate to avoid name collisions with your 1D "point_in_element") ---
function point_in_triangle_3D(P::AbstractVector, Y::AbstractMatrix, nodes::AbstractVector{<:Integer};
                              tol::Float64=1e-10)
    length(nodes) == 3 || error("point_in_triangle_3D: need 3 nodes.")
    A = @view Y[nodes[1], :]
    B = @view Y[nodes[2], :]
    C = @view Y[nodes[3], :]
    lambdas, dist = barycentric_coords_3D(P, A, B, C)
    inside = (dist <= tol) && all(lambdas .>= -tol) && all(lambdas .<= 1.0 + tol)
    return inside, lambdas, dist
end

function Lagrange_interpolation_matrix_surface(X::AbstractMatrix, Y::AbstractMatrix, conn;
                                               dim_u::Int=1, tol::Float64=1e-8, rectangles::Bool=true)
    npts = size(X, 1)
    nnds = size(Y, 1)
    nels = (conn isa AbstractMatrix) ? size(conn, 1) : length(conn)

    I = Int[]; J = Int[]; V = Float64[]

    for r in 1:npts
        Pr = vec(@view X[r, :])
        found = false
        for e in 1:nels
            nodes = _elem_nodes(conn, e)
            nn = length(nodes)

            if nn == 3
                inside, w, _ = point_in_triangle_3D(Pr, Y, nodes; tol=tol)
                inside || continue
                for a in 1:3, k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[a]-1)*dim_u + k + 1)
                    push!(V, w[a])
                end
                found = true
                break

            elseif nn == 4
                inside, N, _ = point_in_quad_Q4(Pr, Y, nodes; tol=tol, rectangles=rectangles)
                inside || continue
                for a in 1:4, k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[a]-1)*dim_u + k + 1)
                    push!(V, N[a])
                end
                found = true
                break
            else
                error("Lagrange_interpolation_matrix_surface: unsupported surface element with $nn nodes.")
            end
        end
        found || error("Lagrange_interpolation_matrix_surface: union node $r not found in base surface mesh.")
    end

    return sparse(I, J, V, dim_u*npts, dim_u*nnds)
end


function build_S_from_elements_surface(fens_u_xyz::AbstractMatrix, conn_u,
                                      fens_frame_xyz::AbstractMatrix, conn_frame;
                                      tol::Float64=1e-8, rectangles::Bool=true, dim_u::Int=1)
    n_u_e = (conn_u isa AbstractMatrix) ? size(conn_u, 1) : length(conn_u)
    n_f_e = (conn_frame isa AbstractMatrix) ? size(conn_frame, 1) : length(conn_frame)

    I = Int[]; J = Int[]; V = Float64[]

    for q in 1:n_u_e
        nodes_q = _elem_nodes(conn_u, q)
        P = vec(sum(fens_u_xyz[nodes_q, :], dims=1) ./ length(nodes_q))

        found = false
        for e in 1:n_f_e
            nodes_e = _elem_nodes(conn_frame, e)
            nn = length(nodes_e)

            inside = false
            if nn == 3
                inside, _, _ = point_in_triangle_3D(P, fens_frame_xyz, nodes_e; tol=tol)
            elseif nn == 4
                inside, _, _ = point_in_quad_Q4(P, fens_frame_xyz, nodes_e; tol=tol, rectangles=rectangles)
            else
                error("build_S_from_elements_surface: unsupported frame element with $nn nodes.")
            end

            if inside
                for k in 0:(dim_u-1)
                    push!(I, (q-1)*dim_u + k + 1)
                    push!(J, (e-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end
                found = true
                break
            end
        end

        found || @warn "build_S_from_elements_surface: union element $q not matched to any frame element (tol=$tol)"
    end

    return sparse(I, J, V, n_u_e*dim_u, n_f_e*dim_u)
end

function build_union_surface_Q4rect(fens_i, fes_i, fens_sd, face_fes; tol::Float64=1e-12)
    nodes0 = _elem_nodes(fes_i.conn, 1)
    length(nodes0) == 4 || error("build_union_surface_Q4rect: frame elements must be Q4.")

    A = vec(@view fens_i.xyz[nodes0[1], :])
    B = vec(@view fens_i.xyz[nodes0[2], :])
    D = vec(@view fens_i.xyz[nodes0[4], :])

    v1 = B .- A
    v2 = D .- A
    M  = hcat(v1, v2)  # 3x2

    nodes_frame = _all_nodes(fes_i.conn)
    nodes_sd    = _all_nodes(face_fes.conn)

    pts = vcat(fens_i.xyz[nodes_frame, :], fens_sd.xyz[nodes_sd, :])

    uv = zeros(size(pts, 1), 2)
    for i in 1:size(pts, 1)
        uv[i, :] .= (M \ (vec(@view pts[i, :]) .- A))  # <-- no transpose
    end

    us = _unique_sorted_tol(uv[:, 1]; tol=tol)
    vs = _unique_sorted_tol(uv[:, 2]; tol=tol)

    length(us) >= 2 || error("build_union_surface_Q4rect: not enough unique u breakpoints.")
    length(vs) >= 2 || error("build_union_surface_Q4rect: not enough unique v breakpoints.")

    fens_uv, fes_u = Q4blockx(us, vs)

    xyz = zeros(size(fens_uv.xyz, 1), 3)
    for i in 1:size(fens_uv.xyz, 1)
        u = fens_uv.xyz[i, 1]
        v = fens_uv.xyz[i, 2]
        xyz[i, :] .= (A .+ u .* v1 .+ v .* v2)          # <-- no transpose
    end

    fens_u3 = FENodeSet(xyz)
    return fens_u3, fes_u
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

# function build_D_matrix(fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8)

#     nodes_i = unique(collect(Iterators.flatten(fes_i.conn[:])))
#     nodes_sd = unique(collect(Iterators.flatten(edge_fes.conn[:])))

#     pts_i  = fens_i.xyz[nodes_i, :]
#     pts_sd = fens_sd.xyz[nodes_sd, :]
#     all_points = vcat(pts_i, pts_sd)

#     unique_points = unique(all_points, dims=1)
#     unique_points_yz = unique_points[:, 2:3]
#     tri_in = Triangulate.TriangulateIO()
#     tri_in.pointlist = unique_points_yz' 
#     tri_out, vor_out = Triangulate.triangulate("Q", tri_in)  



#     # fens_u = FENodeSet(hcat(unique_points[:,1],tri_out.pointlist'))
#     # fes_u = FESetT3(tri_out.trianglelist')

#     fens_u = fens_i
#     fes_u = fes_i

#     File2 = "union_mesh.vtk"
#     vtkexportmesh(
#     File2,
#     fens_u, fes_u,scalars = []
#         )

#     kappa = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] 
#     material = MatHeatDiff(kappa)
#     geom_u = NodalField(fens_u.xyz)
#     u_u = NodalField(zeros(size(fens_u.xyz, 1), 1))
#     numberdofs!(u_u)
#     femm_u = FEMMHeatDiff(IntegDomain(fes_u, TriRule(9)), material)
#     if lam_order == 0
#         M_u = mass_like(femm_u, geom_u, u_u)
#     else
#         M_u = mass(femm_u, geom_u, u_u)
#     end
#     # print(length(fes_u.conn), " elements in the union mesh.\n")
#     # print(size(M_u, 1), " DOFs in the union mesh.\n")

#     p=1
#     Pi_NC = Lagrange_interpolation_matrix(fens_u.xyz, fens_sd.xyz, edge_fes.conn, p; dim_u=1, tol=tol)
#     Pi_phi = Lagrange_interpolation_matrix(fens_u.xyz, fens_i.xyz, fes_i.conn, p; dim_u=1, tol=tol)
#     D = Pi_phi' * M_u * Pi_NC
#     return D, Pi_NC, Pi_phi, M_u
# end

function build_D_matrix(fens_i, fes_i, fens_sd, face_fes;
                        lam_order::Int=0, tol::Float64=1e-8,
                        rectangles::Bool=true,
                        export_union_vtk::Bool=false)

   nn_frame = _nnodes_conn(fes_i.conn)
    nn_sd    = _nnodes_conn(face_fes.conn)


    # ---------------------------------------------------------------------
    # Q4 (quad) skeleton + Q4 (quad) subdomain faces  ---> structured Q4 union
    # ---------------------------------------------------------------------
    if nn_frame == 4 && nn_sd == 4
        fens_u, fes_u = build_union_surface_Q4rect(fens_i, fes_i, fens_sd, face_fes; tol=tol)

        kappa = [1.0 0.0 0.0;
                 0.0 1.0 0.0;
                 0.0 0.0 1.0]
        material = MatHeatDiff(kappa)

        geom_u = NodalField(fens_u.xyz)
        u_u    = NodalField(zeros(size(fens_u.xyz, 1), 1))
        numberdofs!(u_u)

        # Use sufficiently accurate surface quadrature
        femm_u = FEMMHeatDiff(IntegDomain(fes_u, GaussRule(2, 3)), material)

        M_u = (lam_order == 0) ? mass_like(femm_u, geom_u, u_u) : mass(femm_u, geom_u, u_u)

        Pi_NC = Lagrange_interpolation_matrix_surface(fens_u.xyz, fens_sd.xyz, face_fes.conn;
                                                      dim_u=1, tol=tol, rectangles=rectangles)

        if lam_order == 0
            S = build_S_from_elements_surface(fens_u.xyz, fes_u.conn, fens_i.xyz, fes_i.conn;
                                              tol=tol, rectangles=rectangles, dim_u=1)
            D = S' * M_u * Pi_NC
            return D, Pi_NC, S, M_u, fens_u, fes_u
        else
            Pi_phi = Lagrange_interpolation_matrix_surface(fens_u.xyz, fens_i.xyz, fes_i.conn;
                                                           dim_u=1, tol=tol, rectangles=rectangles)
            D = Pi_phi' * M_u * Pi_NC
            return D, Pi_NC, Pi_phi, M_u, fens_u, fes_u
        end
    end

    # ---------------------------------------------------------------------
    # Fallback: your old triangle-based union mesh path (kept)
    # Works when the base surface meshes are triangles (T3).
    # ---------------------------------------------------------------------
    nodes_i  = unique(vec(fes_i.conn))
    nodes_sd = unique(vec(face_fes.conn))

    pts_i  = fens_i.xyz[nodes_i, :]
    pts_sd = fens_sd.xyz[nodes_sd, :]
    all_points = vcat(pts_i, pts_sd)

    unique_points = unique(all_points, dims=1)

    # Existing assumption: interface is parameterized by columns 2:3
    unique_points_yz = unique_points[:, 2:3]
    tri_in = Triangulate.TriangulateIO()
    tri_in.pointlist = unique_points_yz'
    tri_out, _ = Triangulate.triangulate("Q", tri_in)

    fens_u = FENodeSet(hcat(unique_points[:, 1], tri_out.pointlist'))
    fes_u  = FESetT3(tri_out.trianglelist')

    if export_union_vtk
        vtkexportmesh("union_mesh.vtk", fens_u, fes_u, scalars=[])
    end

    kappa = [1.0 0.0 0.0;
             0.0 1.0 0.0;
             0.0 0.0 1.0]
    material = MatHeatDiff(kappa)

    geom_u = NodalField(fens_u.xyz)
    u_u    = NodalField(zeros(size(fens_u.xyz, 1), 1))
    numberdofs!(u_u)

    femm_u = FEMMHeatDiff(IntegDomain(fes_u, TriRule(9)), material)
    M_u = (lam_order == 0) ? mass_like(femm_u, geom_u, u_u) : mass(femm_u, geom_u, u_u)

    # Use surface interpolation (supports triangles and quads, but this branch is mainly for triangles)
    Pi_NC  = Lagrange_interpolation_matrix_surface(fens_u.xyz, fens_sd.xyz, face_fes.conn;
                                                   dim_u=1, tol=tol, rectangles=rectangles)
    Pi_phi = Lagrange_interpolation_matrix_surface(fens_u.xyz, fens_i.xyz,  fes_i.conn;
                                                   dim_u=1, tol=tol, rectangles=rectangles)

    D = Pi_phi' * M_u * Pi_NC
    return D, Pi_NC, Pi_phi, M_u, fens_u, fes_u
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