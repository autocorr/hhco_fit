include("hhco_util.jl")

using Distributions


function get_priors()
    n_samples = 100
    u = LinRange(0, 1, n_samples)

    x_n = LinRange(  4.5,   7.5, n_samples)  # log(cm-3)
    x_N = LinRange( 13.0,  14.5, n_samples)  # log(cm-2)
    x_T = LinRange( 10.0, 300.0, n_samples)  # K
    x_σ = LinRange(  0.5,   8.5, n_samples)  # km/s
    x_v = LinRange(-10.0,  10.0, n_samples)  # km/s
    x_f = LinRange(  0.5,   1.0, n_samples)  # unitless

    y_n = pdf.(Beta(2.0, 3.5), u)
    y_N = pdf.(Beta(4.0, 4.0), u)
    y_T = pdf.(Beta(1.6, 2.5), u)
    y_σ = pdf.(Beta(2.0, 2.0), u)
    y_v = pdf.(Beta(5.0, 5.0), u)
    y_f = pdf.(     Uniform(), u)
    return (
            (x_n, y_n, LABEL_n),
            (x_N, y_N, LABEL_N),
            (x_T, y_T, LABEL_T),
            (x_σ, y_σ, LABEL_σ),
            (x_v, y_v, LABEL_v),
            (x_f, y_f, LABEL_f),
           )
end


function plot_priors(priors)
    fig = new_figure(7, 5)
    grid = fig[1,1] = GridLayout()
    for (i, (x, y, xlabel)) in enumerate(priors)
        i_row = (i - 1) ÷ 2
        i_col = (i - 1) % 2
        ylabel = i == 5 ? "PDF" : ""
        ax = Axis(grid[i_row,i_col],
                  xminorticksvisible=true, yminorticksvisible=true,
                  xminorticks=IntervalsBetween(4), yminorticks=IntervalsBetween(2),
                  yticklabelsvisible=false, yticksvisible=true,
                  xgridvisible=false, ygridvisible=false,
                  xlabel=xlabel, ylabel=ylabel)
        lines!(x, y, color=:orangered)
        band!(x, fill(0, length(x)), y; color=(:orange, 0.25))
        ymin, ymax = extrema(y)
        ylims!(nothing, 1.15 * ymax)
    end
    tight_layout!(grid)
    save_figure("priors")
end

