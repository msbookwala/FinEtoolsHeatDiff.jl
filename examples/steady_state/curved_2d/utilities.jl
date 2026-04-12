using SparseArrays
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule



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

function build_D_matrix(fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8, give_m = false)
    # edge_nodes_sd = unique(collect(Iterators.flatten(edge_fes.conn[:])))
    p = maximum(length.(edge_fes.conn)) - 1
    fens_u, fes_u, M_u = build_union_mesh(fens_i,fes_i, fens_sd, edge_fes, p; lam_order=lam_order)
    X = fens_u.xyz[ :, 1:2]
    
    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[:, 1:2], edge_fes.conn, p)
    if lam_order != 0
        Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, p)
        D = Pi_phi' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, Pi_phi, M_u
        else
            return D, Pi_NC, Pi_phi
        end
    else 
    # R = build_R_from_node_ids(edge_nodes_sd, count(fens_sd); dim_u=1)
        S = build_S_from_elements(fens_u.xyz[:, 1:2], fes_u.conn,
                              fens_i.xyz[:, 1:2], fes_i.conn, 1; tol=tol, dim_u=1)
        D = S' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, S, M_u
        else
            return D, Pi_NC, S
        end
    end
end

function build_D_matrix(fens_u, fes_u, fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8, give_m = false)
    p_sd = maximum(length.(edge_fes.conn)) - 1
    p_i = maximum(length.(fes_i.conn)) - 1
    X = fens_u.xyz[ :, 1:2]

    kappa = [1.0 0; 0 1.0] 
    material = MatHeatDiff(kappa)
    geom_u = NodalField(fens_u.xyz)
    u_u = NodalField(zeros(size(fens_u.xyz, 1), 1))
    numberdofs!(u_u)
    femm_u = FEMMHeatDiff(IntegDomain(fes_u, GaussRule(1, 4)), material)
    if lam_order == 0
        M_u = mass_like(femm_u, geom_u, u_u)
    else
        M_u = mass(femm_u, geom_u, u_u)
    end

    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[:, 1:2], edge_fes.conn, p_sd;partition_of_unity=:inverse)
    if lam_order != 0
        Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, p_i;partition_of_unity=true)
        D = Pi_phi' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, S, M_u
        else
            return D, Pi_NC, Pi_phi
        end
    else 
    # R = build_R_from_node_ids(edge_nodes_sd, count(fens_sd); dim_u=1)
        S = build_S_from_elements(fens_u.xyz[:, 1:2], fes_u.conn,
                              fens_i.xyz[:, 1:2], fes_i.conn, p_i; tol=tol, dim_u=1)
        D = S' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, S, M_u
        else
            return D, Pi_NC, S
        end
    end
end

function trim(points, mask_points;tol = 1e-12, dir = 2)
    q = sortperm(mask_points[:,dir])
    mask_points = mask_points[q, :]
    endpts = [mask_points[1,:], mask_points[end,:]]
    trimmed_points = []
    for i in 1:size(points, 1)
        p = points[i, :]
        if (p[1] >= endpts[1][1] - tol) && (p[1] <= endpts[2][1] + tol) &&
           (p[2] >= endpts[1][2] - tol) && (p[2] <= endpts[2][2] + tol)
            push!(trimmed_points, p)
        end
    end
    return trimmed_points
end

# TODO: parameterise sort and make it agnostic to the direction. here it is in y direction only
function build_union_mesh(fens_i,fes_i, fens_sd, edge_fes, p; lam_order = 0, curved=true, to_trim = false, dir = 2)
        corner_nodes_sd = unique(stack(edge_fes.conn, dims=1)[:,1:2])
        corner_nodes_i = unique(stack(fes_i.conn, dims=1)[:,1:2])
        if to_trim
            trimmed_sd = trim(fens_sd.xyz[corner_nodes_sd, :], fens_i.xyz[corner_nodes_i, :]; tol=1e-10, dir=dir)
            trimmed_i = trim(fens_i.xyz[corner_nodes_i, :], fens_sd.xyz[corner_nodes_sd, :]; tol=1e-10, dir=dir)
            endpoints = unique(vcat(trimmed_i, trimmed_sd), dims=1)
            endpoints = hcat(endpoints...)'
        else
            endpoints = unique(vcat(fens_i.xyz[corner_nodes_i, :], fens_sd.xyz[corner_nodes_sd, :]), dims=1)
        end
        q = sortperm(endpoints[:, dir])
        endpoints = endpoints[q, :]

    if p==1
        fens_u, fes_u = L2blockx2D(endpoints[:, 1], endpoints[:, 2])
    elseif p==2
        fens_u, fes_u = L3blockx2D(endpoints[:, 1], endpoints[:, 2])
    else
        error("build_union_mesh: p=$p not implemented")
    end
    kappa = [1.0 0; 0 1.0] 
    material = MatHeatDiff(kappa)
    geom_u = NodalField(fens_u.xyz)
    u_u = NodalField(zeros(size(fens_u.xyz, 1), 1))
    numberdofs!(u_u)
    femm_u = FEMMHeatDiff(IntegDomain(fes_u, GaussRule(1, 2)), material)
    if lam_order == 0
        M_u = mass_like(femm_u, geom_u, u_u)
    else
        M_u = mass(femm_u, geom_u, u_u)
    end
    return fens_u, fes_u, M_u
end

# function Lagrange_interpolation_matrix(X, Y, conn, p; dim_u=1, tol = 1e-8)
#     eps_end = tol
#     npts = size(X, 1)
#     nnds = size(Y, 1)
#     nels = size(conn, 1)
#     I = Int[]
#     J = Int[]
#     V = Float64[]

#     if p == 1
#         perm = [1, 2]
#     elseif p == 2
#         perm = [1, 3, 2]
#     else
#         perm = collect(1:(p+1))
#     end

#     for r in 1:npts
#         best_elem = 0
#         best_xi = 0.0
#         best_dist = Inf

#         # Pick the element with the minimum projection distance
#         for elem_idx in 1:nels
#             nodes = conn[elem_idx][:]
#             xi, dist, ok = get_xi((Float64(X[r,1]), Float64(X[r,2])), Y, nodes, p; tol=tol)

#             if dist < best_dist
#                 best_dist = dist
#                 best_xi = xi
#                 best_elem = elem_idx
#             end
#         end

#         best_elem == 0 && continue

#         nodes = conn[best_elem][:]
#         xi = best_xi

#         if abs(xi + 1.0) <= eps_end
#             for k in 0:(dim_u-1)
#                 push!(I, (r-1)*dim_u + k + 1)
#                 push!(J, (nodes[1]-1)*dim_u + k + 1)
#                 push!(V, 1.0)
#             end

#         elseif abs(xi - 1.0) <= eps_end
#             # For your connectivity ordering, last endpoint is nodes[2] for p=1,
#             # and nodes[2] also corresponds to xi=+1 after applying perm below.
#             for k in 0:(dim_u-1)
#                 push!(I, (r-1)*dim_u + k + 1)
#                 push!(J, (nodes[2]-1)*dim_u + k + 1)
#                 push!(V, 1.0)
#             end

#         else
#             N = lagrange_1d(xi, p)
#             N = N[perm]

#             for a in 1:(p+1)
#                 for k in 0:(dim_u-1)
#                     push!(I, (r-1)*dim_u + k + 1)
#                     push!(J, (nodes[a]-1)*dim_u + k + 1)
#                     push!(V, N[a])
#                 end
#             end
#         end
#     end

#     return sparse(I, J, V, dim_u*npts, dim_u*nnds)
# end

function Lagrange_interpolation_matrix(X, Y, conn, p;
                                       dim_u=1,
                                       tol=1e-8,
                                       partition_of_unity=true)

    eps_end = tol
    npts = size(X, 1)
    nnds = size(Y, 1)
    nels = size(conn, 1)

    I = Int[]
    J = Int[]
    V = Float64[]

    if p == 1
        perm = [1, 2]
    elseif p == 2
        perm = [1, 3, 2]
    else
        perm = collect(1:(p+1))
    end

    for r in 1:npts
        best_elem = 0
        best_xi = 0.0
        best_dist = Inf

        for elem_idx in 1:nels
            nodes = conn[elem_idx][:]
            xi, dist, ok = get_xi((Float64(X[r,1]), Float64(X[r,2])), Y, nodes, p; tol=tol)

            if dist < best_dist
                best_dist = dist
                best_xi = xi
                best_elem = elem_idx
            end
        end

        best_elem == 0 && continue

        nodes = conn[best_elem][:]
        xi = best_xi

        if best_dist <= tol
            if abs(xi + 1.0) <= eps_end
                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[1]-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end

            elseif abs(xi - 1.0) <= eps_end
                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[2]-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end

            else
                N = lagrange_1d(xi, p)
                N = N[perm]

                for a in 1:(p+1)
                    for k in 0:(dim_u-1)
                        push!(I, (r-1)*dim_u + k + 1)
                        push!(J, (nodes[a]-1)*dim_u + k + 1)
                        push!(V, N[a])
                    end
                end
            end

        else
            if p == 1
                x1, y1 = Y[nodes[1],1], Y[nodes[1],2]
                x2, y2 = Y[nodes[2],1], Y[nodes[2],2]
                px, py = X[r,1], X[r,2]

                d1 = hypot(px - x1, py - y1)
                d2 = hypot(px - x2, py - y2)

                tproj = 0.5 * (xi + 1.0)
                xproj = x1 + tproj * (x2 - x1)
                yproj = y1 + tproj * (y2 - y1)

                d1p = hypot(xproj - x1, yproj - y1)
                d2p = hypot(xproj - x2, yproj - y2)

                if partition_of_unity == true
                    denom = d1p + d2p
                    if denom <= eps(Float64)
                        N1 = 0.5
                        N2 = 0.5
                    else
                        N1 = d2p / denom
                        N2 = d1p / denom
                    end

                elseif partition_of_unity == false
                    denom = d1 + d2
                    if denom <= eps(Float64)
                        N1 = 0.0
                        N2 = 0.0
                    else
                        N1 = d2p / denom
                        N2 = d1p / denom
                    end

                elseif partition_of_unity == :inverse
                    denom = d1p + d2p
                    if denom <= eps(Float64)
                        N1 = 0.0
                        N2 = 0.0
                    else
                        N1 = d2 / denom
                        N2 = d1 / denom
                    end

                else
                    error("partition_of_unity must be true, false, or :inverse")
                end

                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[1]-1)*dim_u + k + 1)
                    push!(V, N1)

                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[2]-1)*dim_u + k + 1)
                    push!(V, N2)
                end

            else
                N = lagrange_1d(xi, p)
                N = N[perm]

                if partition_of_unity == true
                    s = sum(N)
                    if abs(s) > eps(Float64)
                        N ./= s
                    end
                end

                for a in 1:(p+1)
                    for k in 0:(dim_u-1)
                        push!(I, (r-1)*dim_u + k + 1)
                        push!(J, (nodes[a]-1)*dim_u + k + 1)
                        push!(V, N[a])
                    end
                end
            end
        end
    end

    return sparse(I, J, V, dim_u*npts, dim_u*nnds)
end

function lagrange_1d(xi, p)
    xi_n = range(-1.0, 1.0; length = p+1)
    N = ones(p+1)
    for a in 1:(p+1)
        for b in 1:(p+1)
            if a != b
                N[a] *= (xi - xi_n[b]) / (xi_n[a] - xi_n[b])
            end
        end
    end
    return N
end

function lagrange_shapes_deriv_1d(p::Int, xi::Float64)
    if p == 1
        # nodes at [-1, 1]: dN1/dxi = -1/2, dN2/dxi = 1/2
        return [-0.5, 0.5]
    elseif p == 2
        # nodes at [-1, 0, 1]:
        # N0 = xi*(xi-1)/2, N1 = 1 - xi^2, N2 = xi*(xi+1)/2
        # dN = [xi - 1/2, -2*xi, xi + 1/2]
        return [xi - 0.5, -2.0*xi, xi + 0.5]
    else
        xin = range(-1.0, 1.0, length=p+1)
        dN = zeros(p+1)
        @inbounds for a in 1:p+1
            xia = xin[a]
            s = 0.0
            for k in 1:p+1
                k == a && continue
                num = 1.0; den = 1.0
                for b in 1:p+1
                    b == a && continue
                    if b == k
                        den *= (xia - xin[b])
                    else
                        num *= (xi - xin[b])
                        den *= (xia - xin[b])
                    end
                end
                s += num/den
            end
            dN[a] = s
        end
        return dN
    end
end

function curve_map(xi, Y, nodes, p)
    N   = lagrange_1d(xi, p)
    dN  = lagrange_shapes_deriv_1d(p, xi)
    x = 0.0; y = 0.0
    dx = 0.0; dy = 0.0
    # temporary change in order of nodes
    if p ==1
        perm = [1,2]
    elseif p ==2
        perm = [1,3,2]
    end
    N = N[perm]
    dN = dN[perm]
    for a in 1:p+1
        j = nodes[a]
        Xja, Yja = Y[j,1], Y[j,2]
        Na, dNa = N[a], dN[a]
        x   += Na  * Xja;    y   += Na  * Yja
        dx  += dNa * Xja;    dy  += dNa * Yja
    end
    return x, y, dx, dy
end

function get_xi(P, Y, nodes, p;
                xi0=0.0, tol=1e-15, maxiter=200,
                alpha=1e-12, maxstep=0.5)

    # Exact orthogonal projection for linear elements
    if p == 1
        x1, y1 = Y[nodes[1], 1], Y[nodes[1], 2]
        x2, y2 = Y[nodes[2], 1], Y[nodes[2], 2]

        vx = x2 - x1
        vy = y2 - y1
        wx = P[1] - x1
        wy = P[2] - y1

        L2 = vx*vx + vy*vy
        if L2 <= eps(Float64)
            dist = hypot(P[1] - x1, P[2] - y1)
            return 0.0, dist, false
        end

        # Orthogonal projection onto the infinite line
        t = (wx*vx + wy*vy) / L2

        # Clamp to the segment
        tproj = clamp(t, 0.0, 1.0)

        xproj = x1 + tproj * vx
        yproj = y1 + tproj * vy
        dist = hypot(P[1] - xproj, P[2] - yproj)

        # Map [0,1] -> [-1,1]
        xi = 2.0 * tproj - 1.0

        return xi, dist, true
    end

    # Initial guess from chord projection for curved elements
    x1, y1 = Y[nodes[1], 1], Y[nodes[1], 2]
    x2, y2 = Y[nodes[end], 1], Y[nodes[end], 2]
    vx = x2 - x1
    vy = y2 - y1
    wx = P[1] - x1
    wy = P[2] - y1
    L2 = vx*vx + vy*vy
    xi = L2 > eps(Float64) ? (2.0 * clamp((wx*vx + wy*vy) / L2, 0.0, 1.0) - 1.0) : xi0

    for _ in 1:maxiter
        x, y, dx, dy = curve_map(xi, Y, nodes, p)
        rx, ry = x - P[1], y - P[2]
        dist = hypot(rx, ry)

        # Closest-point condition
        b = dx*rx + dy*ry
        if abs(b) <= tol
            return xi, dist, true
        end

        A = dx*dx + dy*dy
        dxi = -b / (A + alpha)

        if abs(dxi) > maxstep
            dxi = sign(dxi) * maxstep
        end

        xi = clamp(xi + dxi, -1.0, 1.0)
    end

    x, y, _, _ = curve_map(xi, Y, nodes, p)
    dist = hypot(x - P[1], y - P[2])
    return xi, dist, false
end

function point_in_element(P, Y, nodes, p; tol)
    xi, d, ok = get_xi(P, Y, nodes, p)
    inside = (ok && d <= tol && xi >= -1.0 - 1e-12 && xi <= 1.0 + 1e-12)
    return inside, xi, d
end
function build_R_from_node_ids(interface_nodes, n_total_nodes; dim_u)
    n_iface = length(interface_nodes)
    nrows = n_iface * dim_u
    ncols = n_total_nodes * dim_u
    I = Int[]; J = Int[]; V = Float64[]
    for (i, node) in enumerate(interface_nodes)
        for k in 0:(dim_u-1)
            push!(I, (i-1)*dim_u + k + 1)
            push!(J, (node-1)*dim_u + k + 1)
            push!(V, 1.0)
        end
    end
    return sparse(I, J, V, nrows, ncols)
end

function build_S_from_elements(fens_u_xyz, conn_u, fens_frame_xyz, conn_frame, p_frame; tol=1e-4, dim_u=1)

    n_u_e = size(conn_u, 1)
    n_f_e = size(conn_frame, 1)
    I = Int[]; J = Int[]; V = Float64[]

    for q in 1:n_u_e
        nodes_q = conn_u[q][:]

        # midpoint of the union element in physical space
        if length(nodes_q) == 2
            A = (Float64(fens_u_xyz[nodes_q[1],1]), Float64(fens_u_xyz[nodes_q[1],2]))
            B = (Float64(fens_u_xyz[nodes_q[2],1]), Float64(fens_u_xyz[nodes_q[2],2]))
            mid = ((A[1] + B[1]) / 2.0, (A[2] + B[2]) / 2.0)
        elseif length(nodes_q) == 3
            # for quadratic union element, use the geometric midpoint on the chord
            A = (Float64(fens_u_xyz[nodes_q[1],1]), Float64(fens_u_xyz[nodes_q[1],2]))
            B = (Float64(fens_u_xyz[nodes_q[2],1]), Float64(fens_u_xyz[nodes_q[2],2]))
            mid = ((A[1] + B[1]) / 2.0, (A[2] + B[2]) / 2.0)
        else
            error("Unsupported union element with $(length(nodes_q)) nodes")
        end

        found = false
        best_e = 0
        best_dist = Inf

        for e in 1:n_f_e
            nodes_e = conn_frame[e][:]

            # project midpoint onto skeleton element and get xi, distance
            xi, d, ok = get_xi(mid, fens_frame_xyz, nodes_e, p_frame; tol=tol)

            # midpoint belongs to this skeleton element if its projection lies in [-1,1]
            # and projected distance is small
            inside = (ok && d <= tol && xi >= -1.0 - 1e-12 && xi <= 1.0 + 1e-12)

            if inside
                for k in 0:(dim_u-1)
                    push!(I, (q-1)*dim_u + k + 1)
                    push!(J, (e-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end
                found = true
                break
            end

            if d < best_dist
                best_dist = d
                best_e = e
            end
        end

        if !found
            # optional fallback: assign to closest skeleton element
            if best_e != 0
                for k in 0:(dim_u-1)
                    push!(I, (q-1)*dim_u + k + 1)
                    push!(J, (best_e-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end
            else
                @warn "build_S_from_elements: union element $q not matched to any frame element (tol=$tol)"
            end
        end
    end

    return sparse(I, J, V, n_u_e*dim_u, n_f_e*dim_u)
end

function build_dual_basis_matrix(interface_nodes, n_whole; dim_u=1)
    n_int = length(interface_nodes)
    n_int >= 2 || error("Need at least 2 interface nodes.")

    I = Int[]
    J = Int[]
    V = Float64[]

    touched = falses(n_whole)

    # Fill interface rows with dual-basis coefficients
    for a in 1:n_int
        gi = interface_nodes[a]   # global row node id
        touched[gi] = true

        if a == 1
            # left end: ψ1 = 2N1 - N2
            nbrs = ((gi, 2.0), (interface_nodes[2], -1.0))
        elseif a == n_int
            # right end: ψn = -N_{n-1} + 2N_n
            nbrs = ((interface_nodes[n_int-1], -1.0), (gi, 2.0))
        else
            # interior: ψi = -N_{i-1} + 2N_i - N_{i+1}
            nbrs = (
                (interface_nodes[a-1], -1.0),
                (gi, 2.0),
                (interface_nodes[a+1], -1.0),
            )
        end

        for k in 0:(dim_u-1)
            row = (gi-1)*dim_u + k + 1
            for (gj, val) in nbrs
                col = (gj-1)*dim_u + k + 1
                push!(I, row)
                push!(J, col)
                push!(V, val)
            end
        end
    end

    # Identity on all non-interface rows
    for node in 1:n_whole
        if !touched[node]
            for k in 0:(dim_u-1)
                idx = (node-1)*dim_u + k + 1
                push!(I, idx)
                push!(J, idx)
                push!(V, 1.0)
            end
        end
    end

    return sparse(I, J, V, dim_u*n_whole, dim_u*n_whole)
end

function dual_basis_matrix_1d(n::Int)
    n >= 2 || error("Need at least 2 interface nodes.")

    I = Int[]
    J = Int[]
    V = Float64[]

    # for i in 2:n-1
    for i in 1:n
        # diagonal
        push!(I, i); push!(J, i); push!(V, 2.0)

        # left neighbor
        if i > 1
            push!(I, i); push!(J, i-1); push!(V, -1.0)
        end

        # right neighbor
        if i < n
            push!(I, i); push!(J, i+1); push!(V, -1.0)
        end
    end

    return sparse(I, J, V, n, n)
end