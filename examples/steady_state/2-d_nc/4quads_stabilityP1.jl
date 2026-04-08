using SparseArrays
using LinearAlgebra
using Plots

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors

include("utilities.jl")

function f(x,y)
    return 1 + x^2 + 2*y^2
end

function lam1(x,y)
    return 1
end

function lam2(x,y)
    return 2
end

function Q!(forceout, XYZ, tangents, feid, qpid)
    forceout[1] = -6.0
    return forceout
end

function p0_mass_matrix_1d(fens, fes; ind=2)
    ne = count(fes)
    m = zeros(ne)

    for (e, conn) in enumerate(fes.conn)
        n1, n2 = conn
        x1 = fens.xyz[n1, ind]
        x2 = fens.xyz[n2, ind]
        m[e] = abs(x2 - x1)
    end

    return spdiagm(0 => m)
end

function blockdiag_sparse(blocks::Vector)
    out = spzeros(0, 0)
    for A in blocks
        out = blockdiag(out, sparse(A))
    end
    return out
end

function compute_beta_for_refinement(mult; lam_order=0)
    println("\n===== mult = $mult =====")

    boundary_boxes = [
        [0.0,0.0, 0.0,1.0], # left
        [1.0,1.0, 0.0,1.0], # right
        [0.0,1.0, 1.0,1.0], # top
        [0.0,1.0, 0.0,0.0], # bottom
    ]

    r = mult
    nelems = [4,2,2,1] .* (2^r)

    fes_all = Any[]
    fens_all = Any[]
    K_all = SparseMatrixCSC{Float64,Int}[]
    F_all = Vector{Float64}[]
    femm_all = Any[]
    T_all = Any[]
    geom_all = Any[]
    dbc_nodes_all = Vector{Int}[]
    material_all = Any[]
    bfes_all = Any[]

    for i in 1:2
        for j in 1:2
            n = ceil(Int, nelems[2*(i-1)+j])
            fens_local, fes_local = T3block(0.5, 0.5, n, n)
            fens_local.xyz[:,1] .+= (i-1)*0.5
            fens_local.xyz[:,2] .+= (j-1)*0.5
            bfes_local = meshboundary(fes_local)

            geom_local = NodalField(fens_local.xyz)
            T_local = NodalField(zeros(size(fens_local.xyz, 1), 1))

            dbc_nodes_local = Int[]
            for box in boundary_boxes
                nodes_in_box = selectnode(fens_local; box=box, inflate=1e-8)
                append!(dbc_nodes_local, nodes_in_box)
            end
            sort!(unique!(dbc_nodes_local))

            for inode in dbc_nodes_local
                setebc!(T_local, [inode], 1, f(fens_local.xyz[inode,1], fens_local.xyz[inode,2]))
            end
            applyebc!(T_local)
            numberdofs!(T_local)

            kappa = [1.0 0; 0 1.0]
            material_local = MatHeatDiff(kappa)

            femm_local = FEMMHeatDiff(IntegDomain(fes_local, TriRule(3)), material_local)
            K_local = conductivity(femm_local, geom_local, T_local)
            K_ff_local = sparse(matrix_blocked(K_local, nfreedofs(T_local), nfreedofs(T_local))[:ff])
            K_fd_local = matrix_blocked(K_local, nfreedofs(T_local), nfreedofs(T_local))[:fd]

            fis_local = ForceIntensity(Float64, 1, Q!)
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

    K = blockdiag_sparse(K_all)
    F = vcat(F_all...)

    edge_fes_all = extract_interface_fes(bfes_all, fens_all, boundary_boxes)

    N_elem_i = 2 * nelems[4]-1

    xs_i1 = 0.5 .* ones(N_elem_i + 1)
    ys_i1 = collect(linearspace(0.0, 1.0, N_elem_i + 1))
    fensi_1, fesi_1 = L2blockx2D(xs_i1, ys_i1)

    xs_i2 = collect(linearspace(0.0, 1.0, N_elem_i + 1))
    ys_i2 = 0.5 .* ones(N_elem_i + 1)
    fensi_2, fesi_2 = L2blockx2D(xs_i2, ys_i2)

    fens_i = [fensi_1, fensi_2]
    fes_i = [fesi_1, fesi_2]

    u_i = Any[]
    F_lam = Vector{Float64}[]
    femm_i = Any[]
    geom_i = Any[]

    for j in 1:2
        if lam_order == 0
            push!(u_i, ElementalField(zeros(count(fes_i[j]), 1)))
        else
            push!(u_i, NodalField(zeros(size(fens_i[j].xyz, 1), 1)))
        end
        numberdofs!(u_i[j])

        femm = FEMMHeatDiff(IntegDomain(fes_i[j], GaussRule(1,4)), MatHeatDiff(reshape([1.0], 1, 1)))
        geom = NodalField(fens_i[j].xyz)
        push!(femm_i, femm)
        push!(geom_i, geom)

        push!(F_lam, zeros(nfreedofs(u_i[j])))
    end

    D = Vector{SparseMatrixCSC{Float64,Int}}()
    multipliers = [-1 -1; -1 1; 1 -1; 1 1]

    for i in 1:4
        D_rows = Vector{SparseMatrixCSC{Float64,Int}}()

        for j in 1:2
            fens_u, fes_u, _ = build_union_mesh(
                fens_i[j], fes_i[j],
                fens_all[i], edge_fes_all[i], 1;
                lam_order=lam_order, to_trim=true, dir=3-j
            )

            Dij, _, _ = build_D_matrix(
                fens_u, fes_u,
                fens_i[j], fes_i[j],
                fens_all[i], edge_fes_all[i];
                lam_order=lam_order
            )

            Dij = sparse(multipliers[i,j] * Dij)

            dnodes = dbc_nodes_all[i]
            if !isempty(dnodes)
                F_lam[j] .+= -Dij[:, dnodes] * gathersysvec(T_all[i], :d)
            end

            freecols = setdiff(1:count(fens_all[i]), dnodes)
            Dij = Dij[:, freecols]

            push!(D_rows, sparse(Dij))
        end

        push!(D, sparse(vcat(D_rows...)))
    end

    D_mat = sparse(hcat(D...))
    B = D_mat

    M1_lg = mass(femm_i[1], geom_i[1], u_i[1])
    M2_lg = mass(femm_i[2], geom_i[2], u_i[2])
    M_lg = blockdiag(M1_lg, M2_lg)

    h = 1.0 / N_elem_i
    M = h * M_lg

    Kfac = cholesky(Symmetric(K))

    # Form Schur complement on multiplier space:
    # S = B * K^{-1} * B'
    X = Kfac \ Matrix(B')
    S = Symmetric(B * X)

    # Since M is diagonal for P0, use cheap scaling
    mdiag = diag(M)
    invsqrtM = Diagonal(1.0 ./ sqrt.(mdiag))

    G = Symmetric(invsqrtM * S * invsqrtM)

    vals = eigvals(G)
    tol = 1e-12
    vals_pos = vals[vals .> tol]
    beta = isempty(vals_pos) ? 0.0 : sqrt(minimum(vals_pos))

    println("beta_h = $beta")

    return beta, h
end

# ============================================================
# parameters
# ============================================================

lam_order = 1
mults = 1:7

betavals = Float64[]
hvals = Float64[]

# ============================================================
# refinement loop
# ============================================================

for mult in mults
    beta, h = compute_beta_for_refinement(mult; lam_order=lam_order)
    push!(betavals, beta)
    push!(hvals, h)
end

# ============================================================
# plot
# ============================================================

using LaTeXStrings
default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)

plot(
    [0:6], betavals,
    yscale = :log10,
    marker = :circle,
    xlabel = L"Refinement factor $r$",
    ylabel = L"$\beta_h$",
    title  = "LBB verification: P1 Lagrange multipliers",
    grid   = true,
    legend = false
)

savefig("4quad_P1_stability.pdf")