using Plots

L2_errors = []
lagrange_errors = []
n = 5
for a in 0:n
    global r = a
    include("linear_quadratic.jl")
    println("Refinement level: $r, Total L2 error: $tot_l2", ", Total Lagrange error: $tot_lag_err")
    push!(L2_errors, tot_l2)
    push!(lagrange_errors, tot_lag_err)
end
plot(0:n, L2_errors, title="L2 Error vs Refinement Level", xlabel="Refinement Level", ylabel="Total L2 Error", yscale=:log10, label = "L2 Error", line = (5, :dash))
plot!(0:n, lagrange_errors, title="L2 Error vs Refinement Level", xlabel="Refinement Level", ylabel="Total L2 Error", yscale=:log10, label = "Lagrange Error", line = (5, :solid))