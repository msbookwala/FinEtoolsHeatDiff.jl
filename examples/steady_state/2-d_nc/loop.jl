using Plots
using LaTeXStrings

default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)
# Plots.scalefontsizes(1/1.2)
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

# plot(0:n, L2_errors0, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = "P0 Temperature", line = (3, :dash), dpi=1000, palette=:lighttest, markershape=[:rect])
# plot!(0:n, lagrange_errors0, label = "P0 LM", line = (3, :solid), dpi=1000, palette=:lighttest, markershape=[:circle])
# plot!(0:n, L2_errors1, label = "P1 Temperature", line = (3, :dash), dpi=1000, palette=:lighttest, markershape=[:cross])
# plot!(0:n, lagrange_errors1, label = "P1 LM", line = (3, :solid), dpi=1000, palette=:lighttest, markershape=[:diamond])
# savefig("p1quadraticRefinement.pdf")

# L2_errors1 = []
# lagrange_errors1 = []
# n = 5
# for a in 0:n
#     l2 = []
#     lag = []
#     for b in [0]
#         global r = a
#         global lam_order = b
#         include("curved_quadratic.jl")
#         println("Refinement level: $r, Total L2 error: $tot_l2", ", Total Lagrange error: $tot_lag_err")
#         push!(l2, tot_l2)
#         push!(lag, tot_lag_err)
#     end
#     push!(L2_errors1, l2)
#     push!(lagrange_errors1, lag)
# end
# L2_errors1 = hcat(L2_errors1...)'
# lagrange_errors1 = hcat(lagrange_errors1...)'

# plot(0:n, L2_errors1, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = ["P0 Temperature" "P2 Temperature"], line = (3, :dash), dpi=1000, palette=:lighttest,legend=:bottomleft, markershape=[:rect :cross])
# # plot!(0:n, lagrange_errors0, label = "P0 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# # plot!(0:n, L2_errors1, label = "Temperature (Quad)", line = (5, :dash), dpi=1000)
# plot!(0:n, lagrange_errors1, label = ["P0 LM" "P2 LM"], line = (3, :solid), dpi=1000, palette=:lighttest, markershape=[:circle :diamond])
# savefig("Refinement.pdf")


# L2_errors =[]
# lagrange_errors = []
# levels=6
# for b in 0:1
#         l2 = []
#         lag = []
# for a in 0:levels
#         global r = a
#         global lam_order = b
#         println("Refinement level: $r")
#         include("4_quads.jl")
#         push!(l2, T_l2)
#         push!(lag, l_l2)
# end
#         push!(L2_errors, l2)
#         push!(lagrange_errors, lag)
# end


# plot(0:levels, L2_errors, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = ["P0 Temperature" "P1 Temperature"], line = (3, :dash), dpi=1000, palette=:lighttest,legend=:bottomleft, markershape=[:rect :cross])
# # plot!(0:n, lagrange_errors0, label = "P0 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# # plot!(0:n, L2_errors1, label = "Temperature (Quad)", line = (5, :dash), dpi=1000)
# plot!(0:levels, lagrange_errors, label = ["P0 LM" "P1 LM"], line = (3, :solid), dpi=1000, palette=:lighttest, markershape=[:circle :diamond])
# savefig("4QuadRefinement.pdf")


L2_errors =[]
lagrange_errors = []
levels=5
for b in 0:1
        l2 = []
        lag = []
for a in 0:levels
        global r = a
        global lam_order = b
        println("Refinement level: $r")
        include("wohlmuth3.jl")
        push!(l2, T_l2)
        push!(lag, l_l2)
end
        push!(L2_errors, l2)
        push!(lagrange_errors, lag)
end


plot(0:levels, L2_errors, title=L"$L^2$ Error vs Refinement Level ", xlabel="Refinement Level", ylabel="Error", yscale=:log10, label = ["P0 Temperature" "P1 Temperature"], line = (3, :dash), dpi=1000, palette=:lighttest,legend=:bottomleft, markershape=[:rect :cross])
# plot!(0:n, lagrange_errors0, label = "P0 LM", line = (5, :solid), dpi=1000, palette=:lighttest)
# plot!(0:n, L2_errors1, label = "Temperature (Quad)", line = (5, :dash), dpi=1000)
plot!(0:levels, lagrange_errors, label = ["P0 LM" "P1 LM"], line = (3, :solid), dpi=1000, palette=:lighttest, markershape=[:circle :diamond])
savefig("Wohlmuth.pdf")
