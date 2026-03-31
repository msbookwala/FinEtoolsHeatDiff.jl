using LinearAlgebra
using SparseArrays
using WriteVTK
using Meshes
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
# ###############################################################################

function clip_polygon(subject::Vector{Vector{Float64}},
                      clip::Vector{Vector{Float64}}; eps=1e-12)

    isempty(subject) && return Vector{Vector{Float64}}[]
    isempty(clip)    && return Vector{Vector{Float64}}[]

    same(p, q) = norm(p - q) <= eps

    function poly_normal(poly)
        n = zeros(3)
        m = length(poly)
        for i in 1:m
            p = poly[i]
            q = poly[mod1(i + 1, m)]
            n += cross(p, q)
        end
        nn = norm(n)
        nn <= eps && error("Degenerate clip polygon")
        return n / nn
    end

    function side(P, A, B, n)
        return dot(cross(B - A, P - A), n)
    end

    function intersect_seg_clipedge(S, E, A, B, n)
        sS = side(S, A, B, n)
        sE = side(E, A, B, n)
        den = sS - sE
        abs(den) <= eps && return nothing
        t = sS / den
        return S + t * (E - S)
    end

    function push_unique!(out, p)
        for q in out
            same(p, q) && return
        end
        push!(out, p)
    end

    function clean(poly)
        out = Vector{Vector{Float64}}()
        for p in poly
            if isempty(out) || !same(out[end], p)
                push!(out, p)
            end
        end
        if length(out) > 1 && same(out[1], out[end])
            pop!(out)
        end

        # remove any remaining repeated vertices
        out2 = Vector{Vector{Float64}}()
        for p in out
            push_unique!(out2, p)
        end
        return out2
    end

    n = poly_normal(clip)
    output = copy(subject)

    for i in 1:length(clip)
        A = clip[i]
        B = clip[mod1(i + 1, length(clip))]
        input = output
        output = Vector{Vector{Float64}}()
        isempty(input) && break

        S = input[end]
        sS = side(S, A, B, n)
        Sin = sS >= -eps

        for E in input
            sE = side(E, A, B, n)
            Ein = sE >= -eps

            if Ein
                if !Sin
                    X = intersect_seg_clipedge(S, E, A, B, n)
                    X !== nothing && push!(output, X)
                end
                push!(output, E)
            elseif Sin
                X = intersect_seg_clipedge(S, E, A, B, n)
                X !== nothing && push!(output, X)
            end

            S = E
            sS = sE
            Sin = Ein
        end

        output = clean(output)
    end

    return clean(output)
end

# ##############################################################################
function get_node_id(x::Vector{Float64}, node_map, XU,
                     ai, ax, IA, JA, VA, 
                     bi, bx, IB, JB, VB; order=1)
    # TODO: cannot use anything in the negative
    return get!(node_map, abs.(round.(x; digits=5))) do
        # create new node
        push!(XU, x)
        
        baryA = barycentre(x, ax)
        nids_a = [size(XU,1) for k in 1:length(ai)]
        nids_b = [size(XU,1) for k in 1:length(bi)]
        push!(IA, nids_a...)
        push!(JA, ai...)
        push!(VA, baryA...)
        if order == 1
            baryB = barycentre(x, bx)
            push!(IB, nids_b...)
            push!(JB, bi...)
            push!(VB, baryB...)
        end
        size(XU,1)

    end
end
# ###############################################################################
struct Grid
    origin::Vector{Float64}
    h::Float64
    cells::Dict{NTuple{3,Int}, Vector{Int}}
end

cell_index(g::Grid, x::Vector{Float64}) = (
    floor(Int, (x[1]-g.origin[1])/g.h),
    floor(Int, (x[2]-g.origin[2])/g.h),
    floor(Int, (x[3]-g.origin[3])/g.h)
)

function aabb(X::Matrix{Float64}, poly_conn::Vector{Int})
    pts = X[poly_conn, :]
    mn = vec(minimum(pts, dims=1))
    mx = vec(maximum(pts, dims=1))
    return mn, mx
end

function build_grid(X::Matrix{Float64}, Conn::Matrix{Int}; h::Float64, pad::Float64=0.0)
    mn = vec(minimum(X, dims=1)) .- pad
    g = Grid(mn, h, Dict{NTuple{3,Int}, Vector{Int}}())
    for k in 1:size(Conn,1)
        tri = Conn[k, :]
        aabb_mn, aabb_mx = aabb(X, tri)
        aabb_mn .-= pad
        aabb_mx .+= pad
        i0 = cell_index(g, aabb_mn)
        i1 = cell_index(g, aabb_mx)
        for i in i0[1]:i1[1], j in i0[2]:i1[2], l in i0[3]:i1[3]
            key = (i,j,l)
            push!(get!(g.cells, key, Int[]), k)
        end
    end
    return g
end

function query_grid(g::Grid, aabb_mn::Vector{Float64}, aabb_mx::Vector{Float64})
    i0 = cell_index(g, aabb_mn)
    i1 = cell_index(g, aabb_mx)

    cand = Int[]
    for i in i0[1]:i1[1], j in i0[2]:i1[2], l in i0[3]:i1[3]
        key = (i,j,l)
        if haskey(g.cells, key)
            append!(cand, g.cells[key])
        end
    end
    sort!(cand)
    unique!(cand)
    return cand
end
# ###############################################################################

function barycentre(x, xas::Matrix)
    if size(xas, 1) == 3
        a = xas[1, :]
        b = xas[2, :]
        c = xas[3, :]

        v0 = b - a
        v1 = c - a

        d00 = dot(v0, v0)
        d01 = dot(v0, v1)
        d11 = dot(v1, v1)

        denom = d00 * d11 - d01 * d01
        @assert denom != 0 "Triangle is degenerate"
        v2 = x - a

        d20 = dot(v2, v0)
        d21 = dot(v2, v1)

        λ2 = (d11 * d20 - d01 * d21) / denom
        λ3 = (d00 * d21 - d01 * d20) / denom
        λ1 = 1.0 - λ2 - λ3

        bary = [λ1, λ2, λ3]

        return bary
    elseif size(xas, 1) == 4
        a = xas[1, :]
        b = xas[2, :]
        c = xas[3, :]
        d = xas[4, :]

        # Bilinear map:
        # X(s,t) = (1-s)(1-t)a + s(1-t)b + st c + (1-s)t d
        #
        # Equivalently:
        # X(s,t) = a + s(b-a) + t(d-a) + s*t(c-b-d+a)

        e1 = b - a
        e2 = d - a
        e3 = c - b - d + a

        # good default for clipped points inside the element
        s = 0.5
        t = 0.5
        for iter in 1:20
            # residual in R^3
            r = a + s*e1 + t*e2 + (s*t)*e3 - x

            # tangent vectors dX/ds and dX/dt in R^3
            Js = e1 + t*e3
            Jt = e2 + s*e3

            # Solve least-squares Newton step:
            # [Js Jt] * δ ≈ r
            J = hcat(Js, Jt)   # 3 x 2
            δ = J \ r          # least-squares, no singular planar 3x3 solve

            s -= δ[1]
            t -= δ[2]

            if norm(δ) < 1e-12 || norm(r) < 1e-12
                break
            end
        end

        λ1 = (1 - s) * (1 - t)
        λ2 = s * (1 - t)
        λ3 = s * t
        λ4 = (1 - s) * t

        return [λ1, λ2, λ3, λ4]
    else
        error("Unsupported element with $(size(xas,1)) nodes")
    end
end 


# ################################################################################
function common_refinement(fensA, fesA, fensB, fesB; h = 0.1, lam_order = 1, tri_order = 1, triangulation_type = "naive" )
                           # connA and connB must be for quads or tri completely for now. hence matrix.
                           # TODO: for mixed quad+tri, add vector of vectors

    XA = fensA.xyz
    connA = stack(fesA.conn, dims=1)
    XB = fensB.xyz
    connB = stack(fesB.conn, dims=1)
    nA = size(connA,1)
    nB = size(connB,1)
    gridB = build_grid(XB, connB; h=h, pad=0.0)

    parentA = Int[]
    parentB = Int[]

    # containers for barycentric points
    IA = Int[]; JA = Int[]; VA = Float64[]
    IB = Int[]; JB = Int[]; VB = Float64[]

    SIA = Int[]; SJA = Int[]; SVA = Float64[]
    SIB = Int[]; SJB = Int[]; SVB = Float64[]


    XU = Vector{Vector{Float64}}()
    connU = Array{Int}[]
    node_map = Dict{Vector{Float64},Int}()
    for i in 1:nA
        count = 0
        ai = connA[i,:]
        aabb_mn, aabb_mx = aabb(XA, ai)
        pad = 0
        aabb_mn .-= pad
        aabb_mx .+= pad
        cands = query_grid(gridB, aabb_mn, aabb_mx)
        # println("aabb_min: $aabb_mn, aabb_mx: $aabb_mx")
        # @info "Cands for element $i in A: $(length(cands))"
        for j in cands
            count += 1
           
            
            bi = connB[j,:]
            ai = connA[i,:]
            ax = XA[connA[i,:], :]
            bx = XB[connB[j,:], :]



            # Clipping: change to custom implementation is needed
            # clipped = clip(PolyArea([(XA[connA[i, k], 1], XA[connA[i, k], 2]) for k in 1:size(connA,2)]), 
            #     PolyArea([(XB[connB[j, k], 1], XB[connB[j, k], 2]) for k in 1:size(connB,2)]), 
            #     SutherlandHodgmanClipping())
            
            clipped = clip_polygon([vec(XA[connA[i, k], :]) for k in 1:size(connA,2)],
                                  [vec(XB[connB[j, k], :]) for k in 1:size(connB,2)])
                             
            if isnothing(clipped)
                continue
            end
            

            # nv = length(clipped.rings[1].vertices)
            nv = length(clipped)
            conn = Vector{Int}()
            for k in 1:nv
                # vs = [clipped.rings[1].vertices[k].coords.x.val, clipped.rings[1].vertices[k].coords.y.val, 0.0]
                
                vs = [clipped[k][1], clipped[k][2], clipped[k][3]]
                push!(conn, get_node_id(vs, node_map, XU,
                                        ai, ax, IA, JA, VA, 
                                        bi, bx, IB, JB, VB; order=lam_order))
            end
            conn = unique(conn)
            nv = length(conn)
            # # Triangulation
            if triangulation_type=="naive"
                if length(conn)>=3
                    for k in 3:nv
                        conn_cuurent = [conn[1], conn[k-1], conn[k]]

                        if tri_order == 2
                            # add midpoints of edges
                            mid12 = (XU[conn[1]] + XU[conn[k-1]]) / 2
                            mid23 = (XU[conn[k-1]] + XU[conn[k]]) / 2
                            mid31 = (XU[conn[k]] + XU[conn[1]]) / 2

                            mid12_id = get_node_id(mid12, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order)
                            mid23_id = get_node_id(mid23, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order)
                            mid31_id = get_node_id(mid31, node_map, XU,
                                                    ai, ax, IA, JA, VA, 
                                                    bi, bx, IB, JB, VB; order=lam_order)

                            conn_cuurent = [conn[1], conn[k-1], conn[k], mid12_id, mid23_id, mid31_id]
                        end
                        push!(connU, conn_cuurent )
                        push!(parentA, i)
                        push!(parentB, j)
                        if lam_order==0
                            push!(IB, size(connU,1))
                            push!(JB, j)
                            push!(VB, 1.0)
                        end
                    end
                end
            elseif triangulation_type=="cp"
                # get midpoint of clipped polygon for use as Steiner point in triangulation
                centroid = zeros(3)
                for k in 1:nv
                    centroid .+= XU[conn[k]]
                end
                centroid ./= nv
                cnid = get_node_id(centroid, node_map, XU,
                            ai, ax, IA, JA, VA, 
                            bi, bx, IB, JB, VB; order=lam_order)
                # triangulate using centroid as Steiner 
                for k in 1:nv
                    conn_cuurent = [conn[k], conn[mod1(k+1, nv)], cnid]

                    if tri_order == 2
                        # add midpoints of edges
                        mid12 = (XU[conn[k]] + XU[conn[mod1(k+1, nv)]]) / 2
                        mid23 = (XU[conn[mod1(k+1, nv)]] + XU[cnid]) / 2
                        mid31 = (XU[cnid] + XU[conn[k]]) / 2

                        mid12_id = get_node_id(mid12, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order)
                        mid23_id = get_node_id(mid23, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order)
                        mid31_id = get_node_id(mid31, node_map, XU,
                                                ai, ax, IA, JA, VA, 
                                                bi, bx, IB, JB, VB; order=lam_order)

                        conn_cuurent = [conn[k], conn[mod1(k+1,nv)], cnid , mid12_id, mid23_id, mid31_id]
                    end
                    push!(connU, conn_cuurent )
                    push!(parentA, i)
                    push!(parentB, j)
                    if lam_order==0
                        push!(IB, size(connU,1))
                        push!(JB, j)
                        push!(VB, 1.0)
                    end
                end
            end
        end
    end
    # making dimensions consistent for C/D
    if tri_order ==1
        fesu = FESetT3(stack(connU, dims=1))
    elseif tri_order == 2
        fesu = FESetT6(stack(connU, dims=1))
    end
    fensu = FENodeSet(stack(XU, dims=1))
    # if tri_order==2
    #     fensu, fesu = T3toT6(fensu, fesu)
    # end

    push!(IA, 1)
    push!(JA, size(XA, 1))
    push!(VA, 0.0)


    PiA = sparse(IA, JA, VA)
    PiB = sparse(IB, JB, VB)

    # PiA = 0
    # PiB = 0

    
    kappa = [1.0 0; 0 1.0] 
    material = MatHeatDiff(kappa)
    geomu = NodalField(fensu.xyz)
    uu = NodalField(zeros(size(fensu.xyz, 1), 1))
    numberdofs!(uu)
    femmu = FEMMHeatDiff(IntegDomain(fesu, TriRule(6)), material)

    if lam_order ==1
        M = mass(femmu, geomu, uu)
    elseif lam_order == 0
        M = mass_like(femmu, geomu, uu)
    end
    C = PiB' * M * PiA

    

    meta = Dict(
        "XA" => XA,
        "connA" => connA,
        "XB" => XB,
        "connB" => connB,
        "fes_u" => fesu,
        "fens_u" => fensu,
        "parentA" => parentA,
        "parentB" => parentB,
        "PiA" => PiA,
        "PiB" => PiB,
        "M" => M,
    )

    return C, meta
end

# ###############################################################################
function write_vtk_tri_mesh(filename::String,
                            X::Matrix{Float64},
                            T::Matrix{Int};
                            cellfields::AbstractDict{<:AbstractString,<:AbstractVector}=Dict{String,AbstractVector}())
    points = permutedims(X)  # 3 x n
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, vec(T[i,:])) for i in 1:size(T,1)]

    vtk_grid(filename, points, cells) do vtk
        for (name, data) in cellfields
            vtk[String(name), VTKCellData()] = data
        end
    end
    return nothing
end
# ###############################################################################


# ------------------------------------------------------------
# Trial
# -------------------------------------------------------------

# XA = [
#     0.0 0.0 0.0;
#     1.0 0.0 0.0;
#     1.0 1.0 0.0;
#     0.0 1.0 0.0;
#     0.25 0.0 0.0;
#     0.75 1.0 0.0;
# ]
# connA = [
#     1 5 6 4;
#     5 2 3 6
# ]
# fes_1 = FESetQ4(connA)
# fens_1 = FENodeSet(XA)

# # # Mesh B: 4 triangles (fan from center)
# XB = [
#     0.0 0.0 0.0;
#     1.0 0.0 0.0;
#     1.0 1.0 0.0;
#     0.0 1.0 0.0;
#     0.75 0.0 0.0;
#     0.25 1.0 0.0;
    
# ]
# connB = [
#     1 5 6 4;
#     5 2 3 6
# ]

# fes_2 = FESetQ4(connB)
# fens_2 = FENodeSet(XB)

# N_elem_1 = 1
# N_elem_2 = 3
# zs_1 = 0.
# ys_1 = collect(linearspace(0.0, 1.0, N_elem_1+1))
# xs_1 = collect(linearspace(0.0, 1.0, N_elem_1+1))
# fens_1, fes_1 = T3blockx(xs_1, ys_1, :a)
# fens_1.xyz = hcat( fens_1.xyz,zs_1*ones(size(fens_1.xyz, 1), 1))

# ys_2 = collect(linearspace(0.0, 1.0, N_elem_2+1))
# xs_2 = collect(linearspace(0.0, 1.0, N_elem_2+1))
# zs_2 = 0.
# fens_2, fes_2 = T3blockx(xs_2, ys_2, :b)
# fens_2.xyz = hcat( fens_2.xyz,zs_2*ones(size(fens_2.xyz, 1), 1))

# XA = fens_1.xyz
# connA = stack(fes_1.conn, dims=1)
# XB = fens_2.xyz
# connB = stack(fes_2.conn, dims=1)

# # error("stop")

# C, meta = common_refinement(fens_1, fes_1, fens_2, fes_2; order=1, h=0.13)

# XU = stack(meta["XU"], dims=1)
# connU = stack(meta["connU"], dims=1)
# write_vtk_tri_mesh("meshU", XU, connU, cellfields=Dict("parentA" => meta["parentA"], "parentB" => meta["parentB"]))

# fA = "meshA"
# write_vtk_tri_mesh(fA, XA, connA, cellfields=Dict("id" => 1:size(connA,1)))
# fB = "meshB"
# write_vtk_tri_mesh(fB, XB, connB, cellfields=Dict("id" => 1:size(connB,1)))