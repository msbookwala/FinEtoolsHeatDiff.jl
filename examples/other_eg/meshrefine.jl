using LinearAlgebra
using SparseArrays
using WriteVTK
using Meshes
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule

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



# ################################################################################
function common_refinement(XA::Matrix{Float64}, connA::Matrix{Int},
                           XB::Matrix{Float64}, connB::Matrix{Int})
                           # connA and connB must be for quads or tri completely for now. hence matrix.
                           # TODO: for mixed quad+tri, add vector of vectors
    nA = size(connA,1)
    nB = size(connB,1)
    gridB = build_grid(XB, connB; h=0.25, pad=0.1)

    parentA = Int[]
    parentB = Int[]

    XU = Vector{Vector{Float64}}()
    connU = Array{Int}[]
    node_map = Dict{Vector{Float64},Int}()
    for i in 1:nA
        ai = connA[i,:]
        aabb_mn, aabb_mx = aabb(XA, ai)
        pad = 0.1
        aabb_mn .-= pad
        aabb_mx .+= pad
        cands = query_grid(gridB, aabb_mn, aabb_mx)
        for j in cands
            print("Processing A:$i, B:$j\n")
            
            bi = connB[j,:]
            ax = XA[connA[i,:], :]
            bx = XB[connB[j,:], :]



            # Clipping: change to custom implementation is needed
            clipped = clip(PolyArea([(XA[connA[i, k], 1], XA[connA[i, k], 2]) for k in 1:size(connA,2)]), 
                PolyArea([(XB[connB[j, k], 1], XB[connB[j, k], 2]) for k in 1:size(connB,2)]), 
                SutherlandHodgmanClipping())
            
            if isnothing(clipped)
                continue
            end
            

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
                push!(parentA, i)
                push!(parentB, j)

            # Triangulation
            elseif  length(conn)>3
                for k in 3:nv
                    push!(connU, [conn[1], conn[k-1], conn[k]])
                    push!(parentA, i)
                    push!(parentB, j)
                end
            end
        end
    end
    return XU, connU, parentA, parentB
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

# # Mesh B: 4 triangles (fan from center)
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

N_elem_1 = 4
N_elem_2 = 5
zs_1 = 0.
ys_1 = collect(linearspace(0.0, 1.0, N_elem_1+1))
xs_1 = collect(linearspace(0.0, 1.0, N_elem_1+1))
fens_1, fes_1 = T3blockx(xs_1, ys_1, :a)
fens_1.xyz = hcat( fens_1.xyz,zs_1*ones(size(fens_1.xyz, 1), 1))

ys_2 = collect(linearspace(0.0, 1.0, N_elem_2+1))
xs_2 = collect(linearspace(0.0, 1.0, N_elem_2+1))
zs_2 = 0.
fens_2, fes_2 = T3blockx(xs_2, ys_2, :a)
fens_2.xyz = hcat( fens_2.xyz,zs_2*ones(size(fens_2.xyz, 1), 1))

XA = fens_1.xyz
connA = stack(fes_1.conn, dims=1)
XB = fens_2.xyz
connB = stack(fes_2.conn, dims=1)

# error("stop")

XU, connU, parentA, parentB = common_refinement(XA, connA, XB, connB)
XU = stack(XU, dims=1)
connU = stack(connU, dims=1)
write_vtk_tri_mesh("meshU", XU, connU, cellfields=Dict("parentA" => parentA, "parentB" => parentB))