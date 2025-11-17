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

function build_D_matrix(fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8)
    # edge_nodes_sd = unique(collect(Iterators.flatten(edge_fes.conn[:])))
    p = maximum(length.(edge_fes.conn)) - 1
    fens_u, fes_u, M_u = build_union_mesh(fens_i,fes_i, fens_sd, edge_fes, p; lam_order=lam_order)
    X = fens_u.xyz[ :, 1:2]
    
    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[:, 1:2], edge_fes.conn, p)
    if lam_order != 0
        Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, p)
        D = Pi_phi' * M_u * Pi_NC
        return D, Pi_NC, Pi_phi
    else 
    # R = build_R_from_node_ids(edge_nodes_sd, count(fens_sd); dim_u=1)
        S = build_S_from_elements(fens_u.xyz[:, 1:2], fes_u.conn,
                              fens_i.xyz[:, 1:2], fes_i.conn, 1; tol=tol, dim_u=1)
        D = S' * M_u * Pi_NC
        return D, Pi_NC, S
    end
end

function build_D_matrix(fens_u, fes_u, fens_i, fes_i, fens_sd, edge_fes; lam_order = 0,tol=1e-8)
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

    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[:, 1:2], edge_fes.conn, p_sd)
    if lam_order != 0
        Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, p_i)
        D = Pi_phi' * M_u * Pi_NC
        return D, Pi_NC, Pi_phi
    else 
    # R = build_R_from_node_ids(edge_nodes_sd, count(fens_sd); dim_u=1)
        S = build_S_from_elements(fens_u.xyz[:, 1:2], fes_u.conn,
                              fens_i.xyz[:, 1:2], fes_i.conn, p_i; tol=tol, dim_u=1)
        D = S' * M_u * Pi_NC
        return D, Pi_NC, S
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

function Lagrange_interpolation_matrix(X, Y, conn, p; dim_u=1, tol = 1e-8)
    eps_end = tol
    npts = size(X, 1)
    nnds = size(Y, 1)
    nels = size(conn, 1)
    I = Int[]
    J = Int[]
    V = Float64[]
    if p ==1
        perm = [1,2]
    elseif p ==2
        perm = [1,3,2]
    end
    for r in 1:npts
        for elem_idx in 1:nels
            nodes = conn[elem_idx][:]
            in_elem, xi, dist = point_in_element((Float64(X[r,1]), Float64(X[r,2])), Y, nodes, p; tol=tol)
            in_elem||continue
            if abs(xi + 1.0) <= eps_end
                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[1]-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end
                break
            elseif abs(xi - 1.0) <= eps_end
                for k in 0:(dim_u-1)
                    push!(I, (r-1)*dim_u + k + 1)
                    push!(J, (nodes[2]-1)*dim_u + k + 1)
                    push!(V, 1.0)
                end
                break
            else
                N = lagrange_1d(xi, p)
                N = N[perm]
                for a in 1:(p+1)
                    for k in 0:(dim_u-1)
                        push!(I, (r -1)*dim_u + k + 1)
                        push!(J, (nodes[a]-1)*dim_u + k + 1)
                        push!(V, N[a])
                        # if a ==1
                        #     push!(V, N[a])
                        # elseif a ==2
                        #     push!(V, N[end])
                        # else
                        #     push!(V, N[a-1])
                        # end
                    end
                end
                break
            end
        end
    end
    # if max value is >1, give error
    # if maximum(V) > 1.0 + tol
    #     error("Lagrange interpolation matrix has values greater than 1.")
    # end

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
    xi = xi0
    for _ in 1:maxiter
        x, y, dx, dy = curve_map(xi, Y, nodes, p)
        rx, ry = x - P[1], y - P[2]
        dist = hypot(rx, ry)
        b = dx*rx + dy*ry
        if abs(b) <= tol
            return xi, dist, true
        end
        A = dx*dx + dy*dy
        dxi = - b / (A + alpha)
        if abs(dxi) > maxstep
            dxi = sign(dxi) * maxstep
        end
        xi = clamp(xi + dxi, -1.0, 1.0)
    end
    x, y, _, _= curve_map(xi, Y, nodes, p)
    dist = hypot(x - P[1], y - P[2])
    return xi, dist, false
end

function point_in_element(P, Y,
                          nodes, p; tol)
    xi, d, ok = get_xi(P, Y, nodes, p)
    return (ok && d <= tol && xi >= -1.0 - 1e-12 && xi <= 1.0 + 1e-12), xi, d
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
        nA = nodes_q[1]; nB = nodes_q[2]
        A = (Float64(fens_u_xyz[nA,1]), Float64(fens_u_xyz[nA,2]))
        B = (Float64(fens_u_xyz[nB,1]), Float64(fens_u_xyz[nB,2]))
        mid = ((A[1]+B[1])/2.0, (A[2]+B[2])/2.0)
        if p_frame == 2
            mid = (Float64(fens_u_xyz[nodes_q[3],1]), Float64(fens_u_xyz[nodes_q[3],2]))
        end
        found = false

        for e in 1:n_f_e
            nodes_e = conn_frame[e][:]
            inside, xi, d = point_in_element(mid, fens_frame_xyz, nodes_e, p_frame; tol=tol)
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
        if !found
            @warn "build_S_from_elements: union element $q not matched to any frame element (tol=$tol)"
        end
    end

    return sparse(I, J, V, n_u_e*dim_u, n_f_e*dim_u)
end
