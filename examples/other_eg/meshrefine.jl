using LinearAlgebra
using SparseArrays
using WriteVTK
using Meshes

function get_node_id(x::Vector{Float64}, node_map, XU)
    # key = qkey(x)
    return get!(node_map, round.(x; digits=5)) do
        # create new node
        # XU = vcat(XU, reshape(x, 1, 3))
        push!(XU, x)
        # Averts = vcat(Averts, zeros(Int, 1, 3))
        # Bverts = vcat(Bverts, zeros(Int, 1, 3))
        # Abary  = vcat(Abary,  zeros(Float64, 1, 3))
        # Bbary  = vcat(Bbary,  zeros(Float64, 1, 3))
        # push!(hasA, false)
        # push!(hasB, false)
        # print(XU)
        size(XU,1) 
    end
end


function common_refinement(XA, connA, XB, connB)
    nA = size(connA,1)
    nB = size(connB,1)
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
    return XU, connU
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
XU, connU = common_refinement(XA, connA, XB, connB)