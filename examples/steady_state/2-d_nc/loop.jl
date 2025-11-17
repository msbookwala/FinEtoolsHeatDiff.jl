using Plots
using LaTeXStrings

default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)
# Plots.scalefontsizes(1.2)
# L2_errors0 = []
# lagrange_errors0 = []
# n = 5
# for a in 0:n
#     global r = a
#     global lam_order = 0
#     include("linear_quadratic.jl")
#     println("Refinement level: $r, Total L2 error: $tot_l2", ", Total Lagrange error: $tot_lag_err")
#     push!(L2_errors0, tot_l2)
#     push!(lagrange_errors0, tot_lag_err)
# end

# L2_errors1 = []
# lagrange_errors1 = []
# n = 5
# for a in 0:n
#     global r = a
#     global lam_order = 1
#     include("linear_quadratic.jl")
#     println("Refinement level: $r, Total L2 error: $tot_l2", ", Total Lagrange error: $tot_lag_err")
#     push!(L2_errors1, tot_l2)
#     push!(lagrange_errors1, tot_lag_err)
# end

# plot(0:n, L2_errors0, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = "Temperature", line = (5, :dash), dpi=1000, palette=:lighttest)
# plot!(0:n, lagrange_errors0, label = "P0 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# # plot!(0:n, L2_errors1, label = "Temperature (Quad)", line = (5, :dash), dpi=1000)
# plot!(0:n, lagrange_errors1, label = "P1 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# savefig("Refinement.png")

L2_errors1 = []
lagrange_errors1 = []
n = 5
for a in 0:n
    l2 = []
    lag = []
    for b in [0,2]
        global r = a
        global lam_order = b
        include("curved_quadratic.jl")
        println("Refinement level: $r, Total L2 error: $tot_l2", ", Total Lagrange error: $tot_lag_err")
        push!(l2, tot_l2)
        push!(lag, tot_lag_err)
    end
    push!(L2_errors1, l2)
    push!(lagrange_errors1, lag)
end
L2_errors1 = hcat(L2_errors1...)'
lagrange_errors1 = hcat(lagrange_errors1...)'

plot(0:n, L2_errors1, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = ["P0 Temperature" "P2 Temperature"], line = (5, :dash), dpi=1000, palette=:lighttest,legend=:bottomleft)
# plot!(0:n, lagrange_errors0, label = "P0 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# plot!(0:n, L2_errors1, label = "Temperature (Quad)", line = (5, :dash), dpi=1000)
plot!(0:n, lagrange_errors1, label = ["P0 LM" "P2 LM"], line = (5, :solid), dpi=1000, palette=:lighttest)
savefig("Refinement.png")