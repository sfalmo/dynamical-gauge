using Plots, JLD2, LaTeXStrings

include("utils.jl")

function plot_histograms(filename="data/current.jld2")
    titles = Dict(
        "1" => "",
        "A" => raw"\hat{A}",
        "vs[i]" => raw"v_i",
        "vs[i]²" => raw"v_i^2",
        "Fints[i]" => raw"f_{\mathrm{int},i}",
        "Fexts[i]" => raw"f_{\mathrm{int},i}",
        "∇A[i]" => raw"\nabla_i \hat{A}",
    )
    for X in ["x", "v"]
        for k in ["", "kin", "int", "ext"]
            for A in ["", "A"]
                Â = A == "" ? "" : raw"\hat{A}(t)"
                X̂ = raw"\hat{" * X * raw"}"
                kk = k == "" ? "" : raw",\mathrm{" * k * raw"}"
                titles[A*"{"*X*"s[i],H"*k*"0}"] = Â * raw"\{" * X̂ * raw"_i(t), H_{0" * kk * raw"}\}"
                titles[A*"vs[i]{"*X*"[i],H"*k*"0}"] = Â * raw"\hat{v}_i(t) \{" * X̂ * raw"_i(t), H_{0" * kk * raw"}\}"
                titles[A*"{"*X*"s[i],H"*k*"0}_finite_diff"] = Â * raw"\{" * X̂ * raw"_i(t), H_{0" * kk * raw"}\}"
                titles[A*"vs[i]{"*X*"s[i],H"*k*"0}_finite_diff"] = Â * raw"\hat{v}_i(t) \{" * X̂ * raw"_i(t), H_{0" * kk * raw"}\}"
            end
        end
    end
    mkpath("plot_all_histograms")
    jldopen(filename, "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        dt = result_histograms.dt
        rs = result_histograms.rs
        dr = result_histograms.dr
        calc_computed!(result_histograms)
        xlims = (ts[begin] - dt / 2, ts[end] + dt / 2)
        ylims = (rs[begin] - dr / 2, rs[end] + dr / 2)
        size = (800, 400)
        xlabel = L"$t / \tau$"
        ylabel = L"$x / \sigma$"
        for (key, hist) in result_histograms.hists
            title = key
            if key in keys(titles)
                title = LaTeXString(raw"$\langle \sum_i \delta(x - \hat{x}_i(t)) " * titles[key] * raw"\rangle$")
            end
            cmin, cmax = extrema(@view hist[:, 2:end])
            cabsmax = max(abs(cmin), abs(cmax))
            cabsmax = max(cabsmax, 0.0001)
            clims = (-cabsmax, cabsmax)
            c = cgrad(:RdBu, rev=true)
            p = heatmap(ts, rs, hist; c, xlims, ylims, clims, size, xlabel, ylabel, title, top_margin=10Plots.px)
            savefig(p, "plot_all_histograms/$(key).pdf")
        end
    end
end

plot_histograms("data/current.jld2")
