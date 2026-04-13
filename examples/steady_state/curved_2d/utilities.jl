using SparseArrays
using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule

function extract_interface_fes(edge_fe_s, fen_s, boxes)
    interface_fes = []
    for i in 1:length(edge_fe_s)
        boundary_fes = []
        for box in boxes
            boundary_fes = union(boundary_fes, selectelem(fen_s[i], edge_fe_s[i], box=box, inflate=1e-8))
        end
        interface_fes_idxi = setdiff(range(1, count(edge_fe_s[i])), boundary_fes)
        push!(interface_fes, subset(edge_fe_s[i], interface_fes_idxi))
    end
    return interface_fes
end

function build_D_matrix(fens_i, fes_i, fens_sd, edge_fes; lam_order = 0, tol = 1e-8, give_m = false)
    p = maximum(length.(edge_fes.conn)) - 1
    fens_u, fes_u, M_u = build_union_mesh(fens_i, fes_i, fens_sd, edge_fes, p; lam_order=lam_order)
    X = fens_u.xyz[:, 1:2]

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
        S = build_S_from_elements(
            fens_u.xyz[:, 1:2], fes_u.conn,
            fens_i.xyz[:, 1:2], fes_i.conn, 1;
            tol=tol, dim_u=1
        )
        D = S' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, S, M_u
        else
            return D, Pi_NC, S
        end
    end
end

function build_D_matrix(fens_u, fes_u, fens_i, fes_i, fens_sd, edge_fes; lam_order = 0, tol = 1e-8, give_m = false)
    p_sd = maximum(length.(edge_fes.conn)) - 1
    p_i = maximum(length.(fes_i.conn)) - 1
    X = fens_u.xyz[:, 1:2]

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

    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[:, 1:2], edge_fes.conn, p_sd)

    if lam_order != 0
        Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, p_i)
        D = Pi_phi' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, Pi_phi, M_u
        else
            return D, Pi_NC, Pi_phi
        end
    else
        S = build_S_from_elements(
            fens_u.xyz[:, 1:2], fes_u.conn,
            fens_i.xyz[:, 1:2], fes_i.conn, p_i;
            tol=tol, dim_u=1
        )
        D = S' * M_u * Pi_NC
        if give_m
            return D, Pi_NC, S, M_u
        else
            return D, Pi_NC, S
        end
    end
end

function trim(points, mask_points; tol = 1e-12, dir = 2)
    q = sortperm(mask_points[:, dir])
    mask_points = mask_points[q, :]
    endpts = [mask_points[1, :], mask_points[end, :]]
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

function build_union_mesh(fens_i, fes_i, fens_sd, edge_fes, p; lam_order = 0, curved = true, to_trim = false, dir = 2)
    p == 1 || error("This build_union_mesh version is only for p=1")

    tol_merge = 1e-10

    function unique_sorted(vals; tol=1e-10)
        v = sort(collect(vals))
        out = Float64[]
        for x in v
            if isempty(out) || abs(x - out[end]) > tol
                push!(out, x)
            end
        end
        return out
    end

    function ordered_points_from_conn(xyz, conn; dir=2)
        isempty(conn) && error("build_union_mesh: empty connectivity")

        ids = Int[]
        for c in conn
            append!(ids, Int.(collect(c)))
        end
        ids = unique(ids)

        pts = Array{Float64}(undef, length(ids), 2)
        for (k, id) in enumerate(ids)
            pts[k, 1] = xyz[id, 1]
            pts[k, 2] = xyz[id, 2]
        end

        q = sortperm(pts[:, dir])
        return pts[q, :]
    end

    # skeleton points in physical space, ordered
    skel_pts = ordered_points_from_conn(fens_i.xyz, fes_i.conn; dir=dir)
    size(skel_pts, 1) >= 2 || error("build_union_mesh: skeleton needs at least 2 points")

    # cumulative chord parameter on skeleton
    s_skel = zeros(size(skel_pts, 1))
    for a in 2:length(s_skel)
        s_skel[a] = s_skel[a-1] + norm(skel_pts[a, :] .- skel_pts[a-1, :])
    end

    # sub-edge points to project to skeleton
    sd_pts = ordered_points_from_conn(fens_sd.xyz, edge_fes.conn; dir=dir)

    if to_trim
        trimmed_sd = trim(sd_pts, skel_pts; tol=1e-10, dir=dir)
        if !isempty(trimmed_sd)
            pts = zeros(length(trimmed_sd), 2)
            for i in 1:length(trimmed_sd)
                pts[i, 1] = trimmed_sd[i][1]
                pts[i, 2] = trimmed_sd[i][2]
            end
            sd_pts = pts
        end
    end

    # project each sub-edge node onto the skeleton polyline
    s_proj = Float64[]
    for a in 1:size(sd_pts, 1)
        P = (Float64(sd_pts[a,1]), Float64(sd_pts[a,2]))

        best_dist = Inf
        best_s = 0.0

        for e in 1:length(fes_i.conn)
            nodes_e = Int.(collect(fes_i.conn[e]))
            xi, d, ok = get_xi(P, fens_i.xyz[:, 1:2], nodes_e, 1; tol=1e-12)

            if d < best_dist
                best_dist = d

                n1 = nodes_e[1]
                n2 = nodes_e[2]

                A1 = fens_i.xyz[n1, 1]
                A2 = fens_i.xyz[n1, 2]
                B1 = fens_i.xyz[n2, 1]
                B2 = fens_i.xyz[n2, 2]

                i1 = findfirst(k -> abs(skel_pts[k,1] - A1) <= 1e-8 && abs(skel_pts[k,2] - A2) <= 1e-8,
                               1:size(skel_pts,1))
                i2 = findfirst(k -> abs(skel_pts[k,1] - B1) <= 1e-8 && abs(skel_pts[k,2] - B2) <= 1e-8,
                               1:size(skel_pts,1))

                if i1 === nothing || i2 === nothing
                    error("build_union_mesh: could not match skeleton endpoints")
                end

                s1 = s_skel[i1]
                s2 = s_skel[i2]

                t = 0.5 * (xi + 1.0)
                best_s = s1 + t * (s2 - s1)
            end
        end

        push!(s_proj, best_s)
    end

    # merged breakpoints in skeleton parameter
    s_union = unique_sorted(vcat(s_skel, s_proj); tol=tol_merge)

    # reconstruct union mesh points ON the skeleton polyline
    xy_union = zeros(length(s_union), 2)
    for a in 1:length(s_union)
        s = s_union[a]

        if s <= s_skel[1] + tol_merge
            xy_union[a, :] .= skel_pts[1, :]
            continue
        elseif s >= s_skel[end] - tol_merge
            xy_union[a, :] .= skel_pts[end, :]
            continue
        end

        placed = false
        for j in 1:length(s_skel)-1
            if s_skel[j] - tol_merge <= s <= s_skel[j+1] + tol_merge
                ds = s_skel[j+1] - s_skel[j]
                t = ds <= eps(Float64) ? 0.0 : (s - s_skel[j]) / ds
                xy_union[a, 1] = (1.0 - t) * skel_pts[j, 1] + t * skel_pts[j+1, 1]
                xy_union[a, 2] = (1.0 - t) * skel_pts[j, 2] + t * skel_pts[j+1, 2]
                placed = true
                break
            end
        end

        placed || error("build_union_mesh: failed to place union point on skeleton")
    end

    fens_u, fes_u = L2blockx2D(xy_union[:, 1], xy_union[:, 2])

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

function lagrange_1d(xi, p)
    xi_n = range(-1.0, 1.0; length = p + 1)
    N = ones(p + 1)
    for a in 1:(p + 1)
        for b in 1:(p + 1)
            if a != b
                N[a] *= (xi - xi_n[b]) / (xi_n[a] - xi_n[b])
            end
        end
    end
    return N
end

function lagrange_shapes_deriv_1d(p::Int, xi::Float64)
    if p == 1
        return [-0.5, 0.5]
    elseif p == 2
        return [xi - 0.5, -2.0 * xi, xi + 0.5]
    else
        xin = range(-1.0, 1.0, length = p + 1)
        dN = zeros(p + 1)
        @inbounds for a in 1:(p + 1)
            xia = xin[a]
            s = 0.0
            for k in 1:(p + 1)
                k == a && continue
                num = 1.0
                den = 1.0
                for b in 1:(p + 1)
                    b == a && continue
                    if b == k
                        den *= (xia - xin[b])
                    else
                        num *= (xi - xin[b])
                        den *= (xia - xin[b])
                    end
                end
                s += num / den
            end
            dN[a] = s
        end
        return dN
    end
end

function curve_map(xi, Y, nodes, p)
    N = lagrange_1d(xi, p)
    dN = lagrange_shapes_deriv_1d(p, xi)

    if p == 1
        perm = [1, 2]
    elseif p == 2
        perm = [1, 3, 2]
    else
        perm = collect(1:(p + 1))
    end

    N = N[perm]
    dN = dN[perm]

    x = 0.0
    y = 0.0
    dx = 0.0
    dy = 0.0

    for a in 1:(p + 1)
        j = nodes[a]
        Xja, Yja = Y[j, 1], Y[j, 2]
        Na, dNa = N[a], dN[a]
        x += Na * Xja
        y += Na * Yja
        dx += dNa * Xja
        dy += dNa * Yja
    end

    return x, y, dx, dy
end

function get_xi(P, Y, nodes, p; xi0 = 0.0, tol = 1e-15, maxiter = 200, alpha = 1e-12, maxstep = 0.5)
    if p == 1
        x1, y1 = Y[nodes[1], 1], Y[nodes[1], 2]
        x2, y2 = Y[nodes[2], 1], Y[nodes[2], 2]

        vx = x2 - x1
        vy = y2 - y1
        wx = P[1] - x1
        wy = P[2] - y1

        L2 = vx * vx + vy * vy
        if L2 <= eps(Float64)
            dist = hypot(P[1] - x1, P[2] - y1)
            return 0.0, dist, false
        end

        t = (wx * vx + wy * vy) / L2
        tproj = clamp(t, 0.0, 1.0)

        xproj = x1 + tproj * vx
        yproj = y1 + tproj * vy
        dist = hypot(P[1] - xproj, P[2] - yproj)

        xi = 2.0 * tproj - 1.0
        return xi, dist, true
    end

    x1, y1 = Y[nodes[1], 1], Y[nodes[1], 2]
    x2, y2 = Y[nodes[end], 1], Y[nodes[end], 2]
    vx = x2 - x1
    vy = y2 - y1
    wx = P[1] - x1
    wy = P[2] - y1
    L2 = vx * vx + vy * vy
    xi = L2 > eps(Float64) ? (2.0 * clamp((wx * vx + wy * vy) / L2, 0.0, 1.0) - 1.0) : xi0

    for _ in 1:maxiter
        x, y, dx, dy = curve_map(xi, Y, nodes, p)
        rx, ry = x - P[1], y - P[2]
        dist = hypot(rx, ry)

        b = dx * rx + dy * ry
        if abs(b) <= tol
            return xi, dist, true
        end

        A = dx * dx + dy * dy
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

function Lagrange_interpolation_matrix(X, Y, conn, p; dim_u = 1, tol = 1e-8, partition_of_unity = true)
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
        perm = collect(1:(p + 1))
    end

    for r in 1:npts
        best_elem = 0
        best_xi = 0.0
        best_dist = Inf

        for elem_idx in 1:nels
            nodes = conn[elem_idx][:]
            xi, dist, ok = get_xi((Float64(X[r, 1]), Float64(X[r, 2])), Y, nodes, p; tol=tol)
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
                for k in 0:(dim_u - 1)
                    push!(I, (r - 1) * dim_u + k + 1)
                    push!(J, (nodes[1] - 1) * dim_u + k + 1)
                    push!(V, 1.0)
                end
            elseif abs(xi - 1.0) <= eps_end
                for k in 0:(dim_u - 1)
                    push!(I, (r - 1) * dim_u + k + 1)
                    push!(J, (nodes[2] - 1) * dim_u + k + 1)
                    push!(V, 1.0)
                end
            else
                N = lagrange_1d(xi, p)
                N = N[perm]
                for a in 1:(p + 1)
                    for k in 0:(dim_u - 1)
                        push!(I, (r - 1) * dim_u + k + 1)
                        push!(J, (nodes[a] - 1) * dim_u + k + 1)
                        push!(V, N[a])
                    end
                end
            end
        else
            if p == 1
                x1, y1 = Y[nodes[1], 1], Y[nodes[1], 2]
                x2, y2 = Y[nodes[2], 1], Y[nodes[2], 2]
                px, py = X[r, 1], X[r, 2]

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

                for k in 0:(dim_u - 1)
                    push!(I, (r - 1) * dim_u + k + 1)
                    push!(J, (nodes[1] - 1) * dim_u + k + 1)
                    push!(V, N1)

                    push!(I, (r - 1) * dim_u + k + 1)
                    push!(J, (nodes[2] - 1) * dim_u + k + 1)
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

                for a in 1:(p + 1)
                    for k in 0:(dim_u - 1)
                        push!(I, (r - 1) * dim_u + k + 1)
                        push!(J, (nodes[a] - 1) * dim_u + k + 1)
                        push!(V, N[a])
                    end
                end
            end
        end
    end

    return sparse(I, J, V, dim_u * npts, dim_u * nnds)
end

function build_R_from_node_ids(interface_nodes, n_total_nodes; dim_u)
    n_iface = length(interface_nodes)
    nrows = n_iface * dim_u
    ncols = n_total_nodes * dim_u
    I = Int[]
    J = Int[]
    V = Float64[]
    for (i, node) in enumerate(interface_nodes)
        for k in 0:(dim_u - 1)
            push!(I, (i - 1) * dim_u + k + 1)
            push!(J, (node - 1) * dim_u + k + 1)
            push!(V, 1.0)
        end
    end
    return sparse(I, J, V, nrows, ncols)
end

function build_S_from_elements(fens_u_xyz, conn_u, fens_frame_xyz, conn_frame, p_frame; tol = 1e-4, dim_u = 1)
    n_u_e = size(conn_u, 1)
    n_f_e = size(conn_frame, 1)

    I = Int[]
    J = Int[]
    V = Float64[]

    for q in 1:n_u_e
        nodes_q = conn_u[q][:]

        if length(nodes_q) == 2
            A = (Float64(fens_u_xyz[nodes_q[1], 1]), Float64(fens_u_xyz[nodes_q[1], 2]))
            B = (Float64(fens_u_xyz[nodes_q[2], 1]), Float64(fens_u_xyz[nodes_q[2], 2]))
            mid = ((A[1] + B[1]) / 2.0, (A[2] + B[2]) / 2.0)
        elseif length(nodes_q) == 3
            A = (Float64(fens_u_xyz[nodes_q[1], 1]), Float64(fens_u_xyz[nodes_q[1], 2]))
            B = (Float64(fens_u_xyz[nodes_q[2], 1]), Float64(fens_u_xyz[nodes_q[2], 2]))
            mid = ((A[1] + B[1]) / 2.0, (A[2] + B[2]) / 2.0)
        else
            error("Unsupported union element with $(length(nodes_q)) nodes")
        end

        found = false
        best_e = 0
        best_dist = Inf

        for e in 1:n_f_e
            nodes_e = conn_frame[e][:]
            xi, d, ok = get_xi(mid, fens_frame_xyz, nodes_e, p_frame; tol=tol)
            inside = (ok && d <= tol && xi >= -1.0 - 1e-12 && xi <= 1.0 + 1e-12)

            if inside
                for k in 0:(dim_u - 1)
                    push!(I, (q - 1) * dim_u + k + 1)
                    push!(J, (e - 1) * dim_u + k + 1)
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
            if best_e != 0
                for k in 0:(dim_u - 1)
                    push!(I, (q - 1) * dim_u + k + 1)
                    push!(J, (best_e - 1) * dim_u + k + 1)
                    push!(V, 1.0)
                end
            else
                @warn "build_S_from_elements: union element $q not matched to any frame element (tol=$tol)"
            end
        end
    end

    return sparse(I, J, V, n_u_e * dim_u, n_f_e * dim_u)
end

function build_dual_basis_matrix(interface_nodes, n_whole; dim_u = 1)
    n_int = length(interface_nodes)
    n_int >= 2 || error("Need at least 2 interface nodes.")

    I = Int[]
    J = Int[]
    V = Float64[]

    touched = falses(n_whole)

    for a in 1:n_int
        gi = interface_nodes[a]
        touched[gi] = true

        if a == 1
            nbrs = ((gi, 2.0), (interface_nodes[2], -1.0))
        elseif a == n_int
            nbrs = ((interface_nodes[n_int - 1], -1.0), (gi, 2.0))
        else
            nbrs = (
                (interface_nodes[a - 1], -1.0),
                (gi, 2.0),
                (interface_nodes[a + 1], -1.0),
            )
        end

        for k in 0:(dim_u - 1)
            row = (gi - 1) * dim_u + k + 1
            for (gj, val) in nbrs
                col = (gj - 1) * dim_u + k + 1
                push!(I, row)
                push!(J, col)
                push!(V, val)
            end
        end
    end

    for node in 1:n_whole
        if !touched[node]
            for k in 0:(dim_u - 1)
                idx = (node - 1) * dim_u + k + 1
                push!(I, idx)
                push!(J, idx)
                push!(V, 1.0)
            end
        end
    end

    return sparse(I, J, V, dim_u * n_whole, dim_u * n_whole)
end

function dual_basis_matrix_1d(n::Int)
    n >= 2 || error("Need at least 2 interface nodes.")

    I = Int[]
    J = Int[]
    V = Float64[]

    for i in 1:n
        push!(I, i)
        push!(J, i)
        push!(V, 2.0)

        if i > 1
            push!(I, i)
            push!(J, i - 1)
            push!(V, -1.0)
        end

        if i < n
            push!(I, i)
            push!(J, i + 1)
            push!(V, -1.0)
        end
    end

    return sparse(I, J, V, n, n)
end