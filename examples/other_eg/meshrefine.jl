using LinearAlgebra
using SparseArrays
using WriteVTK
using Meshes

function get_node_id(x::Vector{Float64}, node_map, XU)
    return get!(node_map, round.(x; digits=5)) do
        # create new node
        push!(XU, x)
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


# ################################################################################
function common_refinement(XA, connA, XB, connB)
    nA = size(connA,1)
    nB = size(connB,1)
    gridB = build_grid(XB, connB; h=0.5, pad=0.1)

    XU = Vector{Vector{Float64}}()
    connU = Array{Int}[]
    node_map = Dict{Vector{Float64},Int}()
    clips = []
    for i in 1:nA
        for j in 1:nB
            ai = connA[i,:]
            bi = connB[j,:]
            ax = XA[connA[i,:], :]
            bx = XB[connB[j,:], :]
            # Clipping: change to custom implementation is needed
            clipped = clip(PolyArea([(XA[connA[i, k], 1], XA[connA[i, k], 2]) for k in 1:size(connA,2)]), 
                PolyArea([(XB[connB[j, k], 1], XB[connB[j, k], 2]) for k in 1:size(connB,2)]), 
                SutherlandHodgmanClipping())
            nv = length(clipped.rings[1].vertices)
            conn = Vector{Int}()
            for k in 1:nv
                vs = [clipped.rings[1].vertices[k].coords.x.val, clipped.rings[1].vertices[k].coords.y.val, 0.0]
                push!(conn, get_node_id(vs, node_map, XU))
            end
            conn = unique(conn)
            nv = length(conn)
            if length(conn) < 3
                continue

            elseif length(conn) == 3
                push!(connU, conn)

            # Triangulation
            elseif  length(conn)>3
                for k in 3:nv
                    push!(connU, [conn[1], conn[k-1], conn[k]])
                end
            end
            
            push!(clips, clipped)
        end
    end
    return XU, connU, gridB
end

# ------------------------------------------------------------
# Trial
# -------------------------------------------------------------

XA = [
    0.0 0.0 0.0;
    1.0 0.0 0.0;
    1.0 1.0 0.0;
    0.0 1.0 0.0;
    0.25 0.0 0.0;
    0.75 1.0 0.0;
]
connA = [
    1 5 6 4;
    5 2 3 6
]

# Mesh B: 4 triangles (fan from center)
XB = [
    0.0 0.0 0.0;
    1.0 0.0 0.0;
    1.0 1.0 0.0;
    0.0 1.0 0.0;
    0.75 0.0 0.0;
    0.25 1.0 0.0;
    
]
connB = [
    1 5 6 4;
    5 2 3 6
]
XU, connU, grid = common_refinement(XA, connA, XB, connB)