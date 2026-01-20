using SparseArrays
using LinearAlgebra
using KrylovKit
using Plots

# ============================================================
# utilities
# ============================================================

function drop_zero_cols(D::SparseMatrixCSC)
    colnnz = diff(D.colptr)
    keep = findall(>(0), colnnz)
    return D[:, keep], keep
end

function smallest_positive_eig(T; howmany=10, reltol=1e-12)
    evals, _, _ = eigsolve(T, howmany, :SR)
    ev = sort(real.(evals))
    tol = reltol * maximum(abs.(ev))
    evp = filter(x -> x > tol, ev)
    return isempty(evp) ? NaN : evp[1]
end

# ============================================================
# P1 mass matrix on 1D interface
# ============================================================

function p1_mass_matrix_1d(fens, fes)
    nn = size(fens.xyz, 1)
    M = spzeros(nn, nn)

    for conn in fes.conn
        # conn is a Tuple{Int,Int}
        n1, n2 = conn

        y1 = fens.xyz[n1, 2]
        y2 = fens.xyz[n2, 2]
        h  = abs(y2 - y1)

        # local P1 mass matrix on a segment
        Me = (h / 6.0) * [2.0 1.0;
                          1.0 2.0]

        # assemble
        M[n1, n1] += Me[1,1]
        M[n1, n2] += Me[1,2]
        M[n2, n1] += Me[2,1]
        M[n2, n2] += Me[2,2]
    end

    return M
end


# ============================================================
# parameters
# ============================================================

lam_order = 1
mults = 0:5

βvals = Float64[]
hvals = Float64[]

# ============================================================
# refinement loop
# ============================================================
# ------------------------------------------------------------
# refinement loop
# ------------------------------------------------------------

for mult in mults
    println("\n===== mult = $mult =====")


using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
using KrylovKit
using Arpack
using SparseArrays
include("utilities.jl")

N_elem1 = 2 * 2^mult
N_elem2 = 3 * 2^mult
N_elem_i = min(N_elem1, N_elem2)
left_m = "q"
right_m = "t"
skew = 0.
lam_order = 1

kappa = [1.0 0; 0 1.0] 
material = MatHeatDiff(kappa)

#########################################################################################
width1 = 0.5
height1 = 1.0
if left_m == "t"
    fens1, fes1 = T3block(width1, height1, floor(Int, N_elem1/2), N_elem1)
    Rule1 = TriRule(1)
else
    fens1, fes1 = Q4block(width1, height1, floor(Int, N_elem1/2), N_elem1)
    Rule1 = GaussRule(2,2)
end

# edge_nodes1 = selectnode(fens1; box=[width1,width1, 0.0,height1], inflate=1e-8)
boundaryfes1 = meshboundary(fes1)
edge_fes1 = subset(boundaryfes1, selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1], inflate=1e-8))

fens1.xyz[:, 1] .+= skew * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5)


geom1 = NodalField(fens1.xyz)
T1 = NodalField(zeros(size(fens1.xyz, 1), 1)) # displacement field

box1 = [0.0,0.0,0.0,0.0]
dbc_nodes1 = selectnode(fens1; box=box1, inflate=1e-8)
for i in dbc_nodes1
    setebc!(T1, [i], 1, -1.0)
end

applyebc!(T1)
numberdofs!(T1)
femm1 = FEMMHeatDiff(IntegDomain(fes1, Rule1), material)
K1 = conductivity(femm1, geom1, T1)
K1_ff = matrix_blocked(K1, nfreedofs(T1), nfreedofs(T1))[:ff]

l1 = selectelem(fens1, meshboundary(fes1), box = [0.0,0.0, 0.0,height1], inflate=1e-8)
el1femm = FEMMBase(IntegDomain(subset(meshboundary(fes1), l1), GaussRule(1,2)))
fi1 = ForceIntensity(Float64[-1.0])
F1 = distribloads(el1femm, geom1, T1, fi1, 2)
F1_ff = vector_blocked(F1, nfreedofs(T1))[:f]


#########################################################################################


width2 = 0.5
height2 = 1.0
if right_m == "t"
    fens2, fes2 = T3block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    Rule2 = TriRule(1)
else
    fens2, fes2 = Q4block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    Rule2 = GaussRule(2,2)
end
# shift the second mesh to the right by 1.0
fens2.xyz[:, 1] .+= 0.5

boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2, box=[0.5,0.5, 0.0,height2], inflate=1e-8))

fens2.xyz[:, 1] .+= skew * (1.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 0.5)

geom2 = NodalField(fens2.xyz)
T2 = NodalField(zeros(size(fens2.xyz, 1), 1)) # displacement field

box2 = [1.0,1.0,0.0,0.0]
dbc_nodes2 = selectnode(fens2; box=box2, inflate=1e-8)
for i in dbc_nodes2
    setebc!(T2, [i], 1, 0.0)
end

applyebc!(T2)

numberdofs!(T2)
femm2 = FEMMHeatDiff(IntegDomain(fes2, Rule2), material)
K2 = conductivity(femm2, geom2, T2)
K2_ff = matrix_blocked(K2, nfreedofs(T2), nfreedofs(T2))[:ff]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]

l2 = selectelem(fens2, meshboundary(fes2), box = [1.0,1.0, 0.0,height2], inflate=1e-8)
el2femm = FEMMBase(IntegDomain(subset(meshboundary(fes2), l2), GaussRule(1,2)))
fi2 = ForceIntensity(Float64[1.0])
F2 = distribloads(el2femm, geom2, T2, fi2, 2)
F2_ff = vector_blocked(F2, nfreedofs(T2))[:f]


##########################################################################################

xs_i = 0.5*ones(N_elem_i+1)
ys_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
fens_i, fes_i = L2blockx2D(xs_i, ys_i)
fens_i.xyz[:, 1] .+= skew * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5)

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 1)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1)) # Lagrange multipliers field
end
numberdofs!(u_i)
femm_i = FEMMHeatDiff(IntegDomain(fes_i, GaussRule(1,2)), material)
D1,_,_ = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
D2,_,_ = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)

D1 = D1[:, setdiff(1:count(fens1), dbc_nodes1)]
D2 = D2[:, setdiff(1:count(fens2), dbc_nodes2)]

A = [K1_ff          zeros(size(K1_ff,1), size(K2_ff,2))    D1';
     zeros(size(K2_ff,1), size(K1_ff,2))     K2_ff          -D2';
     D1               -D2               zeros(size(D1,1), size(D1,1))]


    D1s = sparse(D1)
    D2s = sparse(D2)

    D1r, keep1 = drop_zero_cols(D1s)
    D2r, keep2 = drop_zero_cols(D2s)

    B = [D1s  -D2s]   # constraint operator

    # reduced primal stiffness
    K1s = sparse(K1_ff)
    K2s = sparse(K2_ff)

    # K1r = K1s[keep1, keep1]
    # K2r = K2s[keep2, keep2]

    K = blockdiag(K1s, K2s)
    Kchol = cholesky(Symmetric(K))

    # --------------------------------------------------------
    # multiplier mass matrix (P1 on 1D mesh)
    # --------------------------------------------------------

    Mλ = mass(femm_i, geom_i, u_i)
    Mλchol = cholesky(Symmetric(Mλ))

    # --------------------------------------------------------
    # inf-sup operator:  Mλ^{-1} B K^{-1} B'
    # --------------------------------------------------------

    function apply_T(y)
        tmp = B' * y
        tmp = Kchol \ tmp
        tmp = B * tmp
        return Mλchol \ (Mλchol' \ tmp)
    end

nλ = size(B, 1)

λs, _, _ = eigsolve(apply_T, nλ, 12, :SR)
λs = sort(real.(λs))

tol = 1e-12 * maximum(abs.(λs))
λpos = filter(x -> x > tol, λs)

βh = sqrt(isempty(λpos) ? NaN : λpos[1])


    println("β_h ≈ $βh")

    push!(βvals, βh)
    push!(hvals, 2.0^(-mult))
end

# ============================================================
# plot
# ============================================================

plot(
    hvals, βvals,
    xscale = :log10,
    yscale = :log10,
    marker = :circle,
    xlabel = "h",
    ylabel = "β_h (discrete inf-sup constant)",
    title  = "LBB verification (lam_order = 1, P1 mortar)",
    grid   = true,
    legend = false
)