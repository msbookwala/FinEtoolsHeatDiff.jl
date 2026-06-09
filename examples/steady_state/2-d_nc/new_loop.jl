using Plots
using LaTeXStrings
using Printf

default(fontfamily="Computer Modern", linewidth=2, framestyle=:box)

# ============================================================
# Choose which example to run
# ============================================================
# Set only this block when you want to switch examples.
# The solver file must set the error variables named below

# example_file = "sinxsiny.jl"
# output_file  = "SinXSinYRefinement.pdf"

# example_file = "linear_quadratic.jl"
# output_file  = "QuadRefinement.pdf"

# example_file = "p2_quadratic.jl"
# output_file  = "P0P2_QuadRefinement.pdf"


# example_file = "curved_quadratic.jl"
# output_file  = "curvedRefinement.pdf"

# example_file = "4_quads.jl"
# output_file  = "4QuadRefinement.pdf"

example_file = "wohlmuth3.jl"
output_file  = "WohlmuthRefinement.pdf"

n = 5

lm_orders = [0, 1]  # LM orders to test. Set to [0] to test only P0 LM, etc.
# lm_orders = [0]  # LM orders to test. Set to [0] to test only P0 LM, etc.
# lm_orders = [0,2]  # LM orders to test. Set to [0] to test only P0 LM, etc.

# Error variable names produced by the included file.
# For most two-subdomain examples:
# temperature_error_symbol = :tot_l2
# lagrange_error_symbol    = :tot_lag_err

# For examples such as 4_quads.jl / wohlmuth*.jl, use:
temperature_error_symbol = :T_l2
lagrange_error_symbol    = :l_l2


# ============================================================
# Helpers
# ============================================================
function get_global_value(sym::Symbol)
    return getfield(Main, sym)
end

function run_convergence(example_file, n, lm_orders;
                         temperature_error_symbol=:tot_l2,
                         lagrange_error_symbol=:tot_lag_err)

    L2_errors = Vector{Vector{Float64}}()
    lagrange_errors = Vector{Vector{Float64}}()

    for lo in lm_orders
        l2 = Float64[]
        lag = Float64[]

        for a in 0:n
            global r = a
            global lam_order = lo

            include(example_file)

            tot_l2_current = Float64(get_global_value(temperature_error_symbol))
            tot_lag_current = Float64(get_global_value(lagrange_error_symbol))

            println("LM order: $lo, Refinement level: $r, Total L2 error: $tot_l2_current, Total Lagrange error: $tot_lag_current")

            push!(l2, tot_l2_current)
            push!(lag, tot_lag_current)
        end

        push!(L2_errors, l2)
        push!(lagrange_errors, lag)
    end

    return L2_errors, lagrange_errors
end

function measured_slope(y, idx)
    # Since h halves when r increases by 1:
    # error ≈ C h^p, h_r = h0 * 2^(-r)
    # p = -log(e_{r+1}/e_r)/log(2)
    return -log(y[idx+1] / y[idx]) / log(2)
end

function slope_string(y, idx; digits=2)
    return @sprintf("%.2f", measured_slope(y, idx))
end

function plot_convergence_with_slopes_in_legend(L2_errors, lagrange_errors, lm_orders, n, output_file)
    x = collect(0:n)

    # For n=5, idx=5 uses the segment from r=4 to r=5.
    # If the last point is noisy, change this to idx=4.
    slope_idx = length(x) - 1

    temp0_slope = slope_string(L2_errors[1], slope_idx)
    lag0_slope  = slope_string(lagrange_errors[1], slope_idx)

    p = plot(
        x, L2_errors[1],
        title=L"$L^2$ Error vs Refinement Level ",
        xlabel=L"Refinement Level $r$",
        ylabel="Error",
        yscale=:log10,
        label="Temperature (P$(lm_orders[1]) LM), slope = $(temp0_slope)",
        line=(3, :dash),
        dpi=1000,
        palette=:lighttest,
        legend=:bottomleft,
        markershape=:rect
    )

    if length(lm_orders) >= 2
        temp1_slope = slope_string(L2_errors[2], slope_idx)
        plot!(
            x, L2_errors[2],
            label="Temperature (P$(lm_orders[2]) LM), slope = $(temp1_slope)",
            line=(3, :dash),
            dpi=1000,
            markershape=:cross
        )
    end

    plot!(
        x, lagrange_errors[1],
        label="P$(lm_orders[1]) LM, slope = $(lag0_slope)",
        line=(3, :solid),
        dpi=1000,
        palette=:lighttest,
        markershape=:circle
    )

    if length(lm_orders) >= 2
        lag1_slope = slope_string(lagrange_errors[2], slope_idx)
        plot!(
            x, lagrange_errors[2],
            label="P$(lm_orders[2]) LM, slope = $(lag1_slope)",
            line=(3, :solid),
            dpi=1000,
            palette=:lighttest,
            markershape=:diamond
        )
    end

    savefig(output_file)
    return p
end


# ============================================================
# Run selected example
# ============================================================
L2_errors, lagrange_errors = run_convergence(
    example_file,
    n,
    lm_orders;
    temperature_error_symbol=temperature_error_symbol,
    lagrange_error_symbol=lagrange_error_symbol
)

plot_convergence_with_slopes_in_legend(L2_errors, lagrange_errors, lm_orders, n, output_file)
