using SparseArrays
using LinearAlgebra
using KrylovKit
using Plots


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

function p0_mass_diag_1d(fens, fes)
    ne = count(fes)
    m = Vector{Float64}(undef, ne)

    @inbounds for e in 1:ne
        conn = fes.conn[e]
        n1, n2 = conn
        y1 = fens.xyz[n1, 2]
        y2 = fens.xyz[n2, 2]
        m[e] = abs(y2 - y1)
    end

    return m
end

# ============================================================
# parameters
# ============================================================

lam_order = 0
mults =1:7

betavalss = Float64[]
hvals = Float64[]

# make dir to save results
if !isdir("gpt_results")
    mkpath("gpt_results")
end


# ============================================================
# refinement loop
# ============================================================

for mult in mults
    println("\n===== mult = $mult =====")




N_elem1 = 2 * ceil(Int, 2^mult)
N_elem2 = 2 * ceil(Int, 2^mult)
N_elem_i = min(N_elem1, N_elem2)
left_m = "q"
right_m = "t"
skew = 0.
# lam_order = 1

kappa = [1.0 0; 0 1.0] 
material = MatHeatDiff(kappa)

#########################################################################################
width1 = 0.5
height1 = 1.0
xs1 = collect(linearspace(0.0, width1, ceil(Int, N_elem1/2)))
# ys1 = log10.(logspace(0.0, height1, N_elem1))
ys1 = gradedspace(0.0, height1, N_elem1+1, 1.5)
if left_m == "t"
    fens1, fes1 = T3blockx(xs1, ys1)

    # fens1, fes1 = T3block(width1, height1, floor(Int, N_elem1/2), N_elem1)
    Rule1 = TriRule(1)
else
    fens1, fes1 = Q4blockx(xs1, ys1)
    # fens1, fes1 = Q4block(width1, height1, floor(Int, N_elem1/2), N_elem1)
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


xs2 = collect(linearspace(0.5, 1.0, ceil(Int, N_elem2/2)))
ys2 = gradedspace(0.0, 1.0, N_elem2+1, 0.8)
width2 = 0.5
height2 = 1.0
if right_m == "t"
    # fens2, fes2 = T3block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    fens2, fes2 = T3blockx(xs2, ys2)
    Rule2 = TriRule(1)
else
    # fens2, fes2 = Q4block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    fens2, fes2 = Q4blockx(xs2, ys2)    
    Rule2 = GaussRule(2,2)
end


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

xs_i = 0.5*ones(N_elem_i÷2+1)
ys_i = collect(linearspace(0.0, 1.0, N_elem_i÷2+1))

# ys_i = ys2
# xs_i = 0.5*ones(length(ys_i))

fens_i, fes_i = L2blockx2D(xs_i, ys_i)
fens_i.xyz[:, 1] .+= skew * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5)

fn1 = "gpt_results/sd1_$(N_elem1).vtk"
vtkexportmesh(fn1, fens1, fes1)

fn2 = "gpt_results/sd2_$(N_elem2).vtk"
vtkexportmesh(fn2, fens2, fes2)

fn_i = "gpt_results/sd_i_$(N_elem_i).vtk"
vtkexportmesh(fn_i, fens_i, fes_i)

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 1)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 1)) # Lagrange multipliers field
end
numberdofs!(u_i)
femm_i = FEMMHeatDiff(IntegDomain(fes_i, GaussRule(1,2)), material)
D1, _, _ = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order, tol=1e-8)
D2, _, _ = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order, tol=1e-8)

free1 = setdiff(1:count(fens1), dbc_nodes1)
free2 = setdiff(1:count(fens2), dbc_nodes2)

D1 = D1[:, free1]
D2 = D2[:, free2]

B = hcat(sparse(D1), -sparse(D2))

M1 = mass(femm1, geom1, T1)
M1_ff = matrix_blocked(M1, nfreedofs(T1), nfreedofs(T1))[:ff]

M2 = mass(femm2, geom2, T2)
M2_ff = matrix_blocked(M2, nfreedofs(T2), nfreedofs(T2))[:ff]

K = blockdiag(K1_ff, K2_ff)
M_sd = blockdiag(M1_ff, M2_ff)
KM = K + M_sd

F_KM = cholesky(Symmetric(KM))

h = 1.0 / N_elem_i
mlg_diag = p0_mass_diag_1d(fens_i, fes_i)

dinv = 1.0 ./ sqrt.(h .* mlg_diag)
Bscaled = Diagonal(dinv) * B

X = F_KM \ Matrix(Bscaled')
S = Bscaled * X

evals = eigvals(Symmetric(Matrix(S)))
tol = 1e-14
epos = filter(>(tol), evals)

beta = isempty(epos) ? 0.0 : sqrt(minimum(epos))

push!(betavalss, beta)
push!(hvals, h)

println("beta_h = $beta")
end

# ============================================================
# plot
# ============================================================
using LaTeXStrings
default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)

plot(
    [1:7], [betavalss betavalsr betavalsl],
    label = ["Skeleton" "Left trace" "Right trace"],
    # xscale = :log10,
    yscale = :log10,
    marker = :circle,
    xlabel = L"Refinement factor $r$",
    ylabel = L"$\sigma_{min}$: Lowest Singular Value of $\mathbf{G}$",
    title  = "LBB verification: P0 Lagrange multipliers",
    grid   = true,
    legend = :right,
    markershape=[:cross :circle :diamond],
    line = [:solid :dash :dot],
    # xflip = true
)
savefig("LBB_all.pdf")