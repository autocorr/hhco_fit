push!(LOAD_PATH, expanduser("~/lustre/temp/jadex/Jadex.jl"))

using Jadex

using CairoMakie
using Makie.FileIO
using LaTeXStrings

const LABEL_R = L"\mathrm{Ratio}"
const LABEL_I = L"T_\mathrm{R}\; [\mathrm{K}]"
const LABEL_n = L"\log_{10}(n\; [\mathrm{cm^{-3}}])"
const LABEL_T = L"T_\mathrm{k}\; [\mathrm{K}]"
const LABEL_N = L"\log_{10}(N(\mathrm{para\!\cdot\!H_2CO})\; [\mathrm{cm^{-2}}])"
const LABEL_δ = L"\Delta v\; [\mathrm{km\,s^{-1}}])"
const LABEL_σ = L"\sigma_\mathrm{v}\; [\mathrm{km\,s^{-1}}]"
const LABEL_v = L"v_\mathrm{lsr}\; [\mathrm{km\,s^{-1}}]"
const LABEL_f = L"f_\mathrm{b}"
const LABEL_1 = L"3_{0,3} \rightarrow 2_{0,2}"
const LABEL_2 = L"3_{2,2} \rightarrow 2_{2,1}"
const LABEL_3 = L"3_{2,1} \rightarrow 2_{2,0}"
const LABEL_1r = LaTeXString(raw"3_{0,3} \rightarrow 2_{0,2}")
const LABEL_2r = LaTeXString(raw"3_{2,2} \rightarrow 2_{2,1}")
const LABEL_3r = LaTeXString(raw"3_{2,1} \rightarrow 2_{2,0}")

const FWHM = 2√(2log(2))


function new_figure(size)
    size_pt = 72 .* size
    Figure(resolution=size_pt, fonts=(; regular="CMU Serif"))
end
new_figure(x, y) = new_figure((x, y))

function save_figure(name, fig)
    FileIO.save("$name.pdf", fig, pt_per_unit=1)
    FileIO.save("$name.png", fig, px_per_unit=10)
end
save_figure(name) = save_figure(name, current_figure())

function tight_layout!(grid, n=10)
    colgap!(grid, n)
    rowgap!(grid, n)
end

