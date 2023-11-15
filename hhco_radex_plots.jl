"""
Analyze and plot ratios of the para-H2CO 218 GHz lines for varying physical
conditions.
"""

include("hhco_util.jl")

const DATADIR = expanduser("~/lustre/temp/moldata")
const MOL = Specie("ph2co-h2", datadir=DATADIR)
const BG = galactic_isrf(MOL)


function get_df(n, T, N, δ)
    rdf = RunDef(MOL, density=Dict("p-h2"=>exp10(n)), tkin=T, cdmol=exp10(N),
                 deltav=δ, escprob=βsphere, bg=BG)
    df = get_results(rdf)
    # Select lines for transitions:
    #    3: 3(0,3) - 2(0,2)  218.222 GHz
    #   10: 3(2,2) - 2(2,1)  218.476
    #   13: 3(2,1) - 2(2,0)  218.760
    df[[3,10,13], :]
end


function calc_ratios(df)
    I1, I2, I3 = df[!, :t_rad]
    R21 = I2 / I1
    R31 = I3 / I1
    R21, R31
end


function collect_ratios(n, T, N, δ)
    ratios = calc_ratios.(get_df.(n, T, N, δ))
    hcat(collect.(ratios)...)'
end


function get_ratio_nT_grid(n_series, T_series, N, δ)
    ratios = zeros(length(n_series), length(T_series), 2)
    for (j, T) in enumerate(T_series)
        for (i, n) in enumerate(n_series)
            ratios[i,j,:] .= calc_ratios(get_df(n, T, N, δ))
        end
    end
    ratios
end


function get_intens_nT_grid(n_series, T_series, N, δ)
    intens = zeros(length(n_series), length(T_series), 3)
    for (j, T) in enumerate(T_series)
        for (i, n) in enumerate(n_series)
            df = get_df(n, T, N, δ)
            intens[i,j,:] .= df[!, :t_rad]
        end
    end
    intens
end


function get_axis_kwargs(;nx=4, ny=4)
    return (;
            xminorticksvisible=true, yminorticksvisible=true,
            xminorticks=IntervalsBetween(nx), yminorticks=IntervalsBetween(ny),
           )
end


function add_series_lines!(ax, series, params...)
    ratios = collect_ratios(params...)
    lines!(ax, series, ratios[:,1]; color=:black,
           label=L"\frac{T_\mathrm{R}(%$LABEL_1r)}{T_\mathrm{R}(%$LABEL_2r)}")
    lines!(ax, series, ratios[:,2]; color=:red,
           label=L"\frac{T_\mathrm{R}(%$LABEL_1r)}{T_\mathrm{R}(%$LABEL_3r)}")
end


function plot_ratio_1d_multipanel()
    n =  5.0
    T = 30.0
    N = 13.0
    δ =  0.7
    samples = 50
    n_series = LinRange( 4,     5,  samples)
    T_series = LinRange(10,   100,  samples)
    N_series = LinRange(12,   14.5, samples)
    δ_series = LinRange( 0.3,  2.0, samples)
    # plot
    fig = new_figure(8, 6)
    axis_kwargs = get_axis_kwargs()
    grid = fig[1,1] = GridLayout()
    axes = [Axis(grid[i,j]; axis_kwargs...) for i=1:2, j=1:2]
    # volume density
    ax = axes[1,1]
    ax.xlabel = LABEL_n
    ax.ylabel = LABEL_R
    add_series_lines!(ax, n_series, n_series, T, N, δ)
    axislegend(ax, position=:rb, framevisible=false)
    # temperature
    ax = axes[1,2]
    ax.xlabel = LABEL_T
    add_series_lines!(ax, T_series, n, T_series, N, δ)
    # column density
    ax = axes[2,1]
    ax.xlabel = LABEL_N
    ax.ylabel = LABEL_R
    add_series_lines!(ax, N_series, n, T, N_series, δ)
    # velocity dispersion
    ax = axes[2,2]
    ax.xlabel = LABEL_δ
    add_series_lines!(ax, δ_series, n, T, N, δ_series)
    # save figure
    tight_layout!(grid)
    save_figure("ratio_1d")
end


function plot_ratio_pcolor()
    samples = 60
    δ = 0.7
    N = 13.0
    n_series = LinRange( 4.5,   8, samples)
    T_series = LinRange(10.0, 300, samples)
    ratio_levels  = LinRange(0,  0.6, 13)
    intens_levels = LinRange(0,  6.0, 13)
    # plot
    fig = new_figure(8, 6.5)
    axis_kwargs = get_axis_kwargs(nx=2, ny=5)
    grid = fig[1,1] = GridLayout()
    axes = [Axis(grid[i,j]; axis_kwargs...) for i=1:2, j=1:2]
    # column density high
    ax = axes[1,1]
    ax.ylabel = LABEL_T
    intens = get_intens_nT_grid(n_series, T_series, N, δ)
    contourf!(ax, n_series, T_series, intens[:,:,1]; levels=intens_levels,
              extendhigh=:auto, colormap=:magma)
    tightlimits!(ax)
    text!(ax, 6.6, 265, text=L"T_\mathrm{R}(%$LABEL_1r)", color=:white,
          textsize=12)
    # column density, high
    ax = axes[1,2]
    hm = contourf!(ax, n_series, T_series, intens[:,:,2]; levels=intens_levels,
                   extendhigh=:auto, colormap=:magma)
    tightlimits!(ax)
    text!(ax, 6.6, 265, text=L"T_\mathrm{R}(%$LABEL_2r)", color=:white,
          textsize=12)
    cb = Colorbar(grid[1,3], hm, label=LABEL_I)
    # column density low, R21
    ax = axes[2,1]
    ax.xlabel = LABEL_n
    ax.ylabel = LABEL_T
    ratios = get_ratio_nT_grid(n_series, T_series, N, δ)
    contourf!(ax, n_series, T_series, ratios[:,:,1]; levels=ratio_levels,
              extendhigh=:auto, colormap=:magma)
    tightlimits!(ax)
    text!(ax, 6.6, 235,
          text=L"\frac{T_\mathrm{R}(%$LABEL_1r)}{T_\mathrm{R}(%$LABEL_2r)}",
          color=:black, textsize=12)
    # column density low, R31
    ax = axes[2,2]
    ax.xlabel = LABEL_n
    hm = contourf!(ax, n_series, T_series, ratios[:,:,2]; levels=ratio_levels,
                   extendhigh=:auto, colormap=:magma)
    tightlimits!(ax)
    text!(ax, 6.6, 235,
          text=L"\frac{T_\mathrm{R}(%$LABEL_1r)}{T_\mathrm{R}(%$LABEL_3r)}",
          color=:black, textsize=12)
    cb = Colorbar(grid[2,3], hm, label=LABEL_R)
    # save figure
    tight_layout!(grid)
    save_figure("ratio_2d")
end

