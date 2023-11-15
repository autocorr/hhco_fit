"""
Analyze and plot ratios of the para-H2CO 218 GHz lines for varying physical
conditions.
"""

include("hhco_util.jl")

using Glob
using Logging

using CSV
using FITSIO
using KissMCMC
using Distributions
using LoopVectorization

import Plots
import PairPlots

using PyCall
pyemcee = pyimport("emcee")
corner  = pyimport("corner")
getdist = pyimport("getdist")
gdplots = pyimport("getdist.plots")

const H    = 6.62607015e-27  # erg s, Planck's constant
const KB   = 1.380649e-16    # erg/K, Boltzmann's constant
const TCMB = 2.72548         # K, T(CMB); Fixsen (2009) ApJ 707 916F
const FWHM = 2√(2log(2))     # convert σ to full-width at half-maximum
const F0   = 218.5e9         # Hz, approximate center freq
const T0   = H * F0 / KB
const TBG  = T0 / expm1(T0 / TCMB)
const MAD_TO_STD = 1.4826

const DATADIR = expanduser("~/lustre/temp/moldata")


SYS_VELOS = Dict(
    "G0.068-0.075"   => 50.0,  # km/s
    "G0.106-0.082"   => 55.0,
    "G0.145-0.086"   => 55.0,
    "G0.212-0.001"   => 45.0,
    "G0.316-0.201"   => 20.0,
    "G0.326-0.085"   => 10.0,
    "G0.340+0.055"   => 15.0,
    "G0.380+0.050"   => 40.0,
    "G0.393-0.034"   => 90.0,
    "G0.412+0.052"   => 20.0,
    "G0.483+0.010"   => 25.0,
    "G0.489+0.010"   => 25.0,
    "G1.651-0.050"   => 55.0,
    "G1.683-0.089"   => NaN,
    "G359.137+0.031" => 0.0,
    "G359.615-0.243" => 20.0,
    "G359.734+0.002" => NaN,
    "G359.865+0.022" => 55.0,
    "G359.889-0.093" => NaN,
)


mad(x) = begin μ = median(x); median(abs.(x .- μ)) end
read_tsv(name) = CSV.File(open("data/$name.tsv"), header=false, skipto=7)
get_spec(f::FITS, i, j) = read(f[1], i, j, :)

struct Spec{F <: AbstractFloat, V <: AbstractArray}
    v::V
    I::V
    rms::F
    vsys::F
    function Spec(v::V, I::V; vsys=NaN) where {V <: AbstractArray}
        @assert ndims(v) == 1
        @assert size(v) == size(I)
        v = copy(v)
        I = copy(I)
        I[isnan.(I)] .= 0
        rms = MAD_TO_STD * mad(I)
        # If the velocity is in descending order, reverse.
        if v[2] < v[1]
            reverse!(v)
            reverse!(I)
        end
        new{eltype(I), V}(v, I, rms, vsys)
    end
end
function Spec(h::ImageHDU, i, j)
    v = get_velocity_axis(h)
    T = read(h, i, j, :)
    Spec(v, T)
end
Spec(f::FITS, i, j) = Spec(f[1], i, j)

function jy_to_k(freq, bmaj, bmin)
    θ_maj = bmaj * 60^2  # deg to arcsec
    θ_min = bmin * 60^2  # deg to arcsec
    1.222e3 * 1e3 / (freq^2 * θ_maj * θ_min)
end

function read_field_spectra(path)
    # convert Jy/bm to K
    hdul = FITS(path)
    table = hdul[2]
    columns = FITSIO.colnames(table)
    field = split(basename(path), "_")[1]
    sources = unique([split(f, "_")[1] for f in columns])
    bmaj = read_header(table)["BMAJ"]
    bmin = read_header(table)["BMIN"]
    vsys = SYS_VELOS[field]
    @info "Sources in $field: $(length(sources))"
    spectra_by_source = Dict()
    for source in sources
        source_spectra = []
        for freq in (218.2, 218.5, 218.8)
            v = convert.(Float64, read(table, "$(source)_velocity_h2co_$(freq)GHz"))
            I = convert.(Float64, read(table, "$(source)_intensity_h2co_$(freq)GHz"))
            T = jy_to_k(freq, bmaj, bmin) .* I
            spec = Spec(v, T; vsys=vsys)
            push!(source_spectra, spec)
        end
        spectra_by_source[source] = source_spectra
    end
    field, spectra_by_source
end

function get_velocity_axis(hdu)
    @assert read_key(hdu, "CUNIT3")[1] == "km/s"
    n, _ = read_key(hdu, "NAXIS3")
    δ, _ = read_key(hdu, "CDELT3")
    c, _ = read_key(hdu, "CRVAL3")
    range(c, c+(n-1)*δ, length=n) |> collect
end

function read_spectra_tsv(name="core1")
    f1 = read_tsv("$(name)_218.2GHz")
    f2 = read_tsv("$(name)_218.5GHz")
    f3 = read_tsv("$(name)_218.8GHz")
    (
     Spec(f1.Column1, f1.Column2),
     Spec(f2.Column1, f2.Column2),
     Spec(f3.Column1, f3.Column2),
    )
end

function read_spectra(i, j; field="G0.106-0.082")
    (
     Spec(FITS("data/$field.H2CO.218.2GHz.fits"), i, j),
     Spec(FITS("data/$field.H2CO.218.5GHz.fits"), i, j),
     Spec(FITS("data/$field.H2CO.218.8GHz.fits"), i, j),
    )
end


function predict!(I, v, vcen, σ, τmain, tex)
    @assert size(I) == size(v)
    @turbo for i in eachindex(I)
        τ = τmain * exp(-0.5 * (v[i] - vcen)^2 / σ^2)
        I[i] = (
                (T0 / (exp(T0 / tex) - 1) - TBG)
                * (1 - exp(-τ))
               )
    end
    I
end
predict!(s::Spec, ŝ::Spec, args...) = predict!(s.I, ŝ.v, args...)


function calc_likelihood(y, ŷ, rms)
    @assert size(y) == size(ŷ)
    lnL = zero(eltype(y))
    @turbo for i in eachindex(y)
        lnL += -0.5 * (y[i] - ŷ[i])^2 / rms^2
    end
    lnL
end
calc_likelihood(s::Spec, ŝ::Spec) = calc_likelihood(s.I, ŝ.I, ŝ.rms)


function get_logprob_interp(dataspec; N=20)
    # create work-array spectra from observed
    workspec = deepcopy(dataspec)
    # parameter axes
    log_density = LinRange( 4.5,   7.5, 5N)
    log_cdmol   = LinRange(13.0,  14.5, N)
    tkin        = LinRange(10.0, 300.0, N)
    deltav      = LinRange(0.5FWHM, 8.5FWHM, N)  # ~1 to 20 km/s FWHM
    cdmol = exp10.(log_cdmol)
    density = exp10.(log_density)
    vsys = dataspec[1].vsys
    # Select lines for transitions:
    #    3: 3(0,3) - 2(0,2)  218.222 GHz
    #   10: 3(2,2) - 2(2,1)  218.476
    #   13: 3(2,1) - 2(2,0)  218.760
    trans = [3, 10, 13]
    inlog = [true, false, true, false]
    log_axes = (log_density, tkin, log_cdmol, deltav)
    lin_axes = (density, tkin, cdmol, deltav)
    # compute cube grid
    @info "Calculating grid"
    mol = Specie("ph2co-h2", datadir=DATADIR)
    bg  = galactic_isrf(mol)
    iterpars = Jadex.Solver.IterationParams(5, 500, 4, 4)
    # ignore warnings for temperature interpolation
    logger = SimpleLogger(Logging.Error)
    τcube, xcube = with_logger(logger) do
        rungrid(mol, (lin_axes..., trans); collision_partner="p-h2",
                iterpars=iterpars, escprob=βsphere, bg=bg)
    end
    # compute cube
    @info "Calculating interpolation"
    τinterp = get_interp(τcube, (log_axes..., trans))
    xinterp = get_interp(xcube, (log_axes..., trans))
    function logprob(params)
        log_n, log_N, T, σ, v, f_b = params
        # Calculate prior probabilities
        # Prior ranges:
        #   log_n :     4.5     7.5             # log(cm-3)
        #   log_N :    13.0    15.0             # log(cm-2)
        #       T :    10.0   300.0             # K
        #       σ :     0.5     8.0             # km/s
        #       v :   -20.0    20.0             # km/s
        #     f_b :     0.5     1.0             # unitless
        lnL = zero(log_n)                       # log-likelihood
        lnL += logpdf(  4.5  +   3.0 * Beta(2.0, 3.5), log_n)
        lnL += logpdf( 13.0  +   1.5 * Beta(4.0, 4.0), log_N)
        lnL += logpdf( 10.0  + 290.0 * Beta(1.6, 2.5), T)
        lnL += logpdf(  0.5  +   8.0 * Beta(2.0, 2.0), σ)
        lnL += logpdf(-20.0  +  40.0 * Beta(5.0, 5.0), v)
        lnL += logpdf(  0.5  +   0.5 * Uniform(), f_b)
        if lnL == -Inf
            return lnL
        end
        # Derived parameters
        δv = σ * FWHM               # FWHM linewidth
        vcen = v + vsys             # LSRK velocity
        # Predict spectra
        for (i, (s, ŝ)) in enumerate(zip(workspec, dataspec))
            tau = τinterp(log_n, T, log_N, δv, i)
            tex = xinterp(log_n, T, log_N, δv, i)
            predict!(s, ŝ, vcen, σ, tau, tex)
            @turbo s.I .*= f_b
            lnL += calc_likelihood(s, ŝ)
        end
        lnL
    end
    logprob
end


function get_logprob_direct(dataspec)
    # create work-array spectra from observed
    workspec = deepcopy(dataspec)
    # set-up Jadex workspace
    iterpars = Jadex.Solver.IterationParams(5, 500, 4, 4)
    mol = Specie("ph2co-h2", datadir=DATADIR)
    sol = Solution(mol, iterpars)
    bg  = galactic_isrf(mol)
    trans = [3, 10, 13]
    vsys =  dataspec[1].vsys
    function logprob(params)
        log_n, log_N, T, σ, v, f_b = params
        # Calculate prior probabilities
        # Prior ranges:
        #   log_n :     4.5     7.5             # log(cm-3)
        #   log_N :    13.0    15.0             # log(cm-2)
        #       T :    10.0   300.0             # K
        #       σ :     0.5     8.0             # km/s
        #       v :   -20.0    20.0             # km/s
        #     f_b :     0.5     1.0             # unitless
        lnL = zero(log_n)                       # log-likelihood
        lnL += logpdf(  4.5  +   3.0 * Beta(2.0, 3.5), log_n)
        lnL += logpdf( 13.0  +   1.5 * Beta(4.0, 4.0), log_N)
        lnL += logpdf( 10.0  + 290.0 * Beta(1.6, 2.5), T)
        lnL += logpdf(  0.5  +   8.0 * Beta(2.0, 2.0), σ)
        lnL += logpdf(-20.0  +  40.0 * Beta(5.0, 5.0), v)
        lnL += logpdf(  0.5  +   0.5 * Uniform(), f_b)
        if lnL == -Inf
            return lnL
        end
        # Derived parameters
        δv = σ * FWHM                  # FWHM linewidth
        vcen = v + vsys                # LSRK velocity
        # Predict spectra
        rdf = RunDef(mol, density=Dict("p-h2"=>exp10(log_n)), tkin=T,
                     cdmol=exp10(log_N), deltav=δv, escprob=βsphere, bg=bg)
        solve!(sol, rdf)
        for (i, (s, ŝ)) in enumerate(zip(workspec, dataspec))
            tau = sol.τl[trans[i]]
            tex = sol.tex[trans[i]]
            predict!(s, ŝ, vcen, σ, tau, tex)
            @turbo s.I .*= f_b
            lnL += calc_likelihood(s, ŝ)
        end
        reset!(sol)
        lnL
    end
    logprob
end


function run_mcmc(logprob; niter=5_000)
    initial_params = [
             5.0  # log_n
            14.0  # log_N
           100.0  # T
             3.0  # σ
             1.0  # v
             0.8  # f_b
    ]
    nwalkers = 30
    radii = 0.1 .* initial_params
    θ₀ = make_theta0s(initial_params, radii, logprob, nwalkers)
    chain, accept = emcee(logprob, θ₀, niter=niter, nthin=5, use_progress_meter=true)
    squash_chain, squash_accept = squash_walkers(chain, accept)
    @info "Acceptance ratio: $squash_accept"
    stack(squash_chain, dims=1), squash_accept
end


function run_pymcmc(logprob; niter=500)
    initial_params = [
             5.0  # log_n
            14.0  # log_N
           100.0  # T
             3.0  # σ
             1.0  # v
             0.8  # f_b
    ]
    nwalkers = 30
    ndim = 6
    nburn = niter ÷ 5
    radii = 0.1 .* initial_params
    θ₀ = make_theta0s(initial_params, radii, logprob, nwalkers)
    sampler = pyemcee.EnsembleSampler(nwalkers, ndim, logprob)
    state = sampler.run_mcmc(θ₀, nburn)
    sampler.reset()
    sampler.run_mcmc(nothing, niter)
    chain = sampler.get_chain(flat=true)
    chain, sampler
end


function get_predicted_profile(params, dataspec)
    workspec = deepcopy(dataspec)
    mol = Specie("ph2co-h2", datadir=DATADIR)
    bg = galactic_isrf(mol)
    log_n, log_N, T, σ, v, f_b = params
    δv = σ * FWHM
    vcen = v + dataspec[1].vsys
    rdf = RunDef(mol, density=Dict("p-h2"=>exp10(log_n)), tkin=T,
                 cdmol=exp10(log_N), deltav=δv, escprob=βsphere, bg=bg)
    _, sol = solve(rdf)
    trans = [3, 10, 13]
    for (i, (s, ŝ)) in enumerate(zip(workspec, dataspec))
        tau = sol.τl[trans[i]]
        tex = sol.tex[trans[i]]
        predict!(s, ŝ, vcen, σ, tau, tex)
        @turbo s.I .*= f_b
    end
    workspec
end


function plot_spectra(spectra; params=nothing, xlim=(-50, 150),
        ylim=(-1, 3), outname="test")
    fig = new_figure(4, 5)
    grid = fig[1,1] = GridLayout()
    labels = [LABEL_1, LABEL_2, LABEL_3]
    if ~isnothing(params)
        predicted_spectra = get_predicted_profile(params, spectra)
    end
    for (i, spec) in enumerate(spectra)
        xlabel = i == 3 ? L"v_\mathrm{lsr} \; [ \mathrm{km\,s^{-1}} ]" : ""
        ylabel = i == 3 ? L"T_\mathrm{b} \; [ \mathrm{K} ]" : ""
        ax = Axis(grid[i,1],
                  xminorticksvisible=true, xminorticks=IntervalsBetween(5),
                  xgridvisible=false, ygridvisible=false, xlabel=xlabel,
                  ylabel=ylabel, xlabelsize=12, ylabelsize=12,
                  xticklabelsize=10, yticklabelsize=10)
        hlines!([0], color=:gray, linestyle=:dash, linewidth=0.7)
        vlines!([spec.vsys-20, spec.vsys+20], color=:gray, linestyle=:dot, linewidth=0.7)
        band!(spec.v, zero(spec.I), spec.I, color=(:orange, 0.25))
        stairs!(spec.v, spec.I, color=:orangered, step=:center)
        label_x = xlim[1] + 0.05xlim[2]
        label_y = ylim[1] + 0.80ylim[2]
        text!(label_x, label_y, text=labels[i], fontsize=12)
        xlims!(xlim...)
        ylims!(ylim...)
        if ~isnothing(params)
            pspec = predicted_spectra[i]
            lines!(pspec.v, pspec.I, color=:red)
        end
    end
    tight_layout!(grid)
    save_figure("$(outname)_spectra")
end


function plot_test_spectra()
    coords = [(110, 80), (100, 63), (194, 116)]
    for (i, j) in coords
        @info "plotting $i, $j"
        spectra = read_spectra(i, j)
        plot_spectra(spectra; outname="$(i)_$(j)")
    end
end


function astable(chain; nthin=1)
    (;
     logn=chain[1:nthin:end,1],
     logN=chain[1:nthin:end,2],
        T=chain[1:nthin:end,3],
        s=chain[1:nthin:end,4],
        v=chain[1:nthin:end,5],
       fb=chain[1:nthin:end,6],
    )
end


function plot_corner(chain; nthin=1, outname="corner")
    Plots.gr()
    table = astable(chain; nthin=nthin)
    PairPlots.corner(table; plotcontours=false, filterscatter=true, dpi=300)
    Plots.savefig("$outname.pdf")
    Plots.savefig("$outname.png")
end


function plot_pycorner(chain; outname="corner_py")
    labels = [LABEL_n, LABEL_N, LABEL_T, LABEL_σ, LABEL_v, LABEL_f]
    corner.corner(chain, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=true)
    corner.core.pl.savefig("$outname.pdf", dpi=300)
    corner.core.pl.savefig("$outname.png", dpi=300)
    corner.core.pl.close("all")
end


function plot_corner_gd(chain; outname="corner_gd")
    names = ["logn", "logN", "T", "s", "v", "fb"]
    samples = getdist.MCSamples(samples=chain, names=names)
    g = gdplots.get_subplot_plotter()
    g.triangle_plot(samples, filled=true)
    gdplots.plt.savefig("$outname.pdf", dpi=300)
    gdplots.plt.savefig("$outname.png", dpi=300)
    gdplots.plt.close("all")
end


function plot_tau_profile()
    v = LinRange(-5, 5, 1_000)
    σ = 1.0
    vcen = 0.0
    τmain = 10.0
    τ = @. τmain * exp(-0.5 * (v - vcen)^2 / σ^2)
    y = @. 1 - exp(-τ)
    fig = new_figure(4, 3)
    ax = Axis(fig[1,1], xlabel="velocity", ylabel="optical depth")
    hlines!([0], color=:gray, linestyle=:dash, linewidth=0.7)
    lines!(v, y, color=:dodgerblue)
    save_figure("tau_profile")
end


function plot_all_integrated_fields(data_dir="data/mina_data")
    for path in glob(joinpath(data_dir, "*.fits"))
        field, spectra_by_source = read_field_spectra(path)
        n_spectra = length(spectra_by_source)
        workspec = deepcopy(spectra_by_source["$(field)a"])
        for w in workspec
            w.I .= 0
        end
        for (source, spectra) in spectra_by_source
            for (w, s) in zip(workspec, spectra)
                w.I .+= s.I
            end
        end
        for w in workspec
            w.I ./= n_spectra
        end
        plot_spectra(workspec; ylim=(-0.1, 1.0), outname="$(field)_average")
    end
end


function run_all_mina_fits(data_dir="data/mina_data"; niter=50_000)
    for path in glob(joinpath(data_dir, "*.fits"))
        field, spectra_by_source = read_field_spectra(path)
        for (source, spectra) in spectra_by_source
            if isnan(spectra[1].vsys)
                continue
            end
            @info "Fitting source: $source"
            logprob = get_logprob_direct(spectra)
            chain, _ = run_mcmc(logprob; niter=niter)
            plot_pycorner(chain, outname="$(source)_corner_py")
            med_params = median(chain, dims=1)
            xlim = (spectra[1].vsys - 30, spectra[1].vsys + 30)
            plot_spectra(spectra; params=med_params, xlim=xlim, outname="$(source)")
            CSV.write("$(source)_chain.csv", astable(chain))
        end
    end
end


function run_test_fits()
    coords = [(110, 80), (100, 63), (194, 116)]
    for (i, j) in coords
        @info "fitting $i, $j"
        spectra = read_spectra(i, j)
        logprob = get_logprob_direct(spectra)
        chain, _ = run_mcmc(logprob; niter=500_000)
        plot_pycorner(chain, outname="$(i)_$(j)_corner_py")
        CSV.write("$(i)_$(j)_chain.csv", astable(chain))
    end
end

