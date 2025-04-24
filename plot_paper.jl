using PyCall, PyPlot, LaTeXStrings, JLD2, DelimitedFiles
ioff()

include("utils.jl")

rcParams = PyDict(matplotlib["rcParams"])
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["font.size"] = 9
rcParams["text.usetex"] = true
rcParams["figure.dpi"] = 300


function analytical()
    fig, ax = subplots(3, 3, figsize=(7, 5), sharex=true, sharey=true, layout="constrained")
    cmap = "RdBu_r"

    cbar_ρ = nothing
    cbar_Cext = nothing
    cbar_Ckin = nothing
    for (i, case) in enumerate(["ideal_harmonic_eq", "ideal_harmonic_switch", "ideal_harmonic_inverted"])
        xs = vec(readdlm("data/analytical/$(case)/xs.dat"))
        ts = vec(readdlm("data/analytical/$(case)/ts.dat"))
        ts ./= ts[end]
        ρ = readdlm("data/analytical/$(case)/rho.dat")
        Cext = readdlm("data/analytical/$(case)/Cext.dat")
        Ckin = readdlm("data/analytical/$(case)/Ckin.dat")
        # ρ_absmax = maximum(abs, ρ)
        ρ_absmax = 0.4
        im = ax[1, i].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        if i == 3
            cbar_ρ = fig.colorbar(im, ax=ax[1, i], pad=0.01)
            cbar_ρ.ax.set_ylim(0, ρ_absmax)
        end
        Ckin_absmax = 0.4
        im = ax[2, i].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        if i == 3
            cbar_Ckin = fig.colorbar(im, ax=ax[2, i], pad=0.01)
        end
        Cext_absmax = 0.4
        im = ax[3, i].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        if i == 3
            cbar_Cext = fig.colorbar(im, ax=ax[3, i], pad=0.01)
        end
        ax[1, i].contour(ts, xs, ρ, levels=range(0, ρ_absmax, 10), linewidths=0.5, cmap="Greys")
        ax[2, i].contour(ts, xs, Ckin, levels=range(-Ckin_absmax, Ckin_absmax, 11), linewidths=0.5, cmap="Greys")
        ax[3, i].contour(ts, xs, Cext, levels=range(-Cext_absmax, Cext_absmax, 11), linewidths=0.5, cmap="Greys")
    end
    cbar_ρ.set_label(L"\rho a_0")
    cbar_Cext.set_label(L"\langle \hat{C}_\mathrm{ext} \hat{R} \rangle a_0")
    cbar_Ckin.set_label(L"\langle \hat{C}_\mathrm{id} \hat{R} \rangle a_0")

    for a in ax[:, begin]
        a.set_ylabel(L"x / a_0")
    end
    for a in ax[end, :]
        a.set_xlabel(L"t \omega / \pi")
    end

    ax[begin, 1].set_title(L"k = k_0")
    ax[begin, 2].set_title(L"k = k_0 / 4")
    ax[begin, 3].set_title(L"k = -k_0")

    mkpath("plot_paper")
    fig.savefig("plot_paper/analytical.pdf")
end


function compare()
    fig, ax = subplots(3, 3, figsize=(7, 5), sharex=true, sharey=true, layout="constrained")
    cmap = "RdBu_r"

    case = "ideal_harmonic_switch"
    xs = vec(readdlm("data/analytical/$(case)/xs.dat"))
    ts = vec(readdlm("data/analytical/$(case)/ts.dat"))
    ts ./= ts[end]
    ρ_analytical = readdlm("data/analytical/$(case)/rho.dat")
    Cext_analytical = readdlm("data/analytical/$(case)/Cext.dat")
    Ckin_analytical = readdlm("data/analytical/$(case)/Ckin.dat")
    ρ_absmax = 0.4
    ax[1, 1].pcolormesh(ts, xs, ρ_analytical, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
    Ckin_absmax = 0.4
    ax[2, 1].pcolormesh(ts, xs, Ckin_analytical, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
    Cext_absmax = 0.4
    ax[3, 1].pcolormesh(ts, xs, Cext_analytical, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)

    jldopen("data/ideal_harmonic_switch_finite_diff.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        ts ./= ts[end]
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"] / 5
        Cext = result_histograms.hists["ACext"]
        Ckin = result_histograms.hists["ACkin"]
        ax[1, 2].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        ax[2, 2].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        ax[3, 2].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        Ckin_finite_diff = result_histograms.hists["ACkin_finite_diff"]
        Cext_finite_diff = result_histograms.hists["ACext_finite_diff"]
        im = ax[1, 3].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        cbar_ρ = fig.colorbar(im, ax=ax[1, 3])
        cbar_ρ.ax.set_ylim(0, ρ_absmax)
        im = ax[2, 3].pcolormesh(ts, xs, Ckin_finite_diff, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        cbar_Ckin = fig.colorbar(im, ax=ax[2, 3])
        im = ax[3, 3].pcolormesh(ts, xs, Cext_finite_diff, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        cbar_Cext = fig.colorbar(im, ax=ax[3, 3])
        cbar_ρ.set_label(L"\rho a_0")
        cbar_Cext.set_label(L"\langle \hat{C}_\mathrm{ext} \hat{R} \rangle a_0")
        cbar_Ckin.set_label(L"\langle \hat{C}_\mathrm{id} \hat{R} \rangle a_0")
    end

    for a in ax[:, begin]
        a.set_ylabel(L"x / a_0")
    end
    for a in ax[end, :]
        a.set_xlabel(L"t \omega / \pi")
        a.set_ylim(-5, 5)
    end

    ax[begin, 1].set_title("analytical")
    ax[begin, 2].set_title("autodiff")
    ax[begin, 3].set_title("finite difference")

    mkpath("plot_paper")
    fig.savefig("plot_paper/compare.pdf")
end


function doublewell()
    fig, ax = subplots(4, 2, figsize=(3.5, 6), sharex=true, sharey=true, layout="constrained")
    cmap = "RdBu_r"
    Ckin_absmax = 0.6
    Cext_absmax = 0.6
    Cint_absmax = 0.6
    ρ_absmax = 0.6

    jldopen("data/ideal_doublewell_eq.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        Ckin = result_histograms.hists["ACkin"]
        Cext = result_histograms.hists["ACext"]
        Cint = result_histograms.hists["ACint"]
        im = ax[1, 1].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        im = ax[2, 1].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        im = ax[3, 1].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        im = ax[4, 1].pcolormesh(ts, xs, Cint, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
    end

    jldopen("data/WCA_doublewell_eq.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        Ckin = result_histograms.hists["ACkin"]
        Cext = result_histograms.hists["ACext"]
        Cint = result_histograms.hists["ACint"]
        im = ax[1, 2].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        cbar_ρ = fig.colorbar(im, ax=ax[1, 2])
        cbar_ρ.ax.set_ylim(0, ρ_absmax)
        im = ax[2, 2].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        cbar_Ckin = fig.colorbar(im, ax=ax[2, 2])
        im = ax[3, 2].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        cbar_Cext = fig.colorbar(im, ax=ax[3, 2])
        im = ax[4, 2].pcolormesh(ts, xs, Cint, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
        cbar_Cint = fig.colorbar(im, ax=ax[4, 2])
        cbar_ρ.set_label(L"\rho a")
        cbar_Ckin.set_label(L"\langle \hat{C}_\mathrm{id} \hat{R} \rangle a")
        cbar_Cext.set_label(L"\langle \hat{C}_\mathrm{ext} \hat{R} \rangle a")
        cbar_Cint.set_label(L"\langle \hat{C}_\mathrm{int} \hat{R} \rangle a")
    end

    for a in ax[:, begin]
        a.set_ylabel(L"x / a")
        a.set_ylim(-5, 5)
    end
    for a in ax[end, :]
        a.set_xticks([0, 1, 2, 3])
        a.set_xlabel(L"t / t_\mathrm{MD}")
    end

    ax[begin, 1].set_title("ideal gas")
    ax[begin, 2].set_title("WCA")
    ax[begin, 1].set_xlim(0, 3.15)
    ax[begin, 2].set_xlim(0, 3.15)

    mkpath("plot_paper")
    fig.savefig("plot_paper/doublewell.pdf")
end


function interacting()
    fig, ax = subplots(4, 4, figsize=(7, 6), sharex=true, sharey=true, layout="constrained")
    cmap = "RdBu_r"
    Ckin_absmax = 0.4
    Cext_absmax = 0.4
    Cint_absmax = 0.4
    ρ_absmax = 0.8

    jldopen("data/3WCA_harmonic_stronger_finite_diff.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        Ckin = result_histograms.hists["ACkin"]
        Cext = result_histograms.hists["ACext"]
        Cint = result_histograms.hists["ACint"]
        Ckin_finite_diff = result_histograms.hists["ACkin_finite_diff"]
        Cext_finite_diff = result_histograms.hists["ACext_finite_diff"]
        Cint_finite_diff = result_histograms.hists["ACint_finite_diff"]
        ax[1, 1].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        ax[2, 1].pcolormesh(ts, xs, Ckin_finite_diff, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        ax[3, 1].pcolormesh(ts, xs, Cext_finite_diff, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        ax[4, 1].pcolormesh(ts, xs, Cint_finite_diff, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
        ax[1, 2].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        ax[2, 2].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        ax[3, 2].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        ax[4, 2].pcolormesh(ts, xs, Cint, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
    end

    jldopen("data/5WCA_harmonic_stronger_finite_diff.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        Ckin = result_histograms.hists["ACkin"]
        Cext = result_histograms.hists["ACext"]
        Cint = result_histograms.hists["ACint"]
        im = ax[1, 3].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        im = ax[2, 3].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        im = ax[3, 3].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        im = ax[4, 3].pcolormesh(ts, xs, Cint, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
    end

    jldopen("data/5WCA_harmonic-doublewell_finite_diff.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        Ckin = result_histograms.hists["ACkin"]
        Cext = result_histograms.hists["ACext"]
        Cint = result_histograms.hists["ACint"]
        im = ax[1, end].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        cbar_ρ = fig.colorbar(im, ax=ax[1, end])
        cbar_ρ.ax.set_ylim(0, ρ_absmax)
        im = ax[2, end].pcolormesh(ts, xs, Ckin, cmap=cmap, vmin=-Ckin_absmax, vmax=Ckin_absmax, rasterized=true)
        cbar_Ckin = fig.colorbar(im, ax=ax[2, end])
        im = ax[3, end].pcolormesh(ts, xs, Cext, cmap=cmap, vmin=-Cext_absmax, vmax=Cext_absmax, rasterized=true)
        cbar_Cext = fig.colorbar(im, ax=ax[3, end])
        im = ax[4, end].pcolormesh(ts, xs, Cint, cmap=cmap, vmin=-Cint_absmax, vmax=Cint_absmax, rasterized=true)
        cbar_Cint = fig.colorbar(im, ax=ax[4, end])
        cbar_ρ.set_label(L"\rho a")
        cbar_Ckin.set_label(L"\langle \hat{C}_\mathrm{id} \hat{R} \rangle a")
        cbar_Cext.set_label(L"\langle \hat{C}_\mathrm{ext} \hat{R} \rangle a")
        cbar_Cint.set_label(L"\langle \hat{C}_\mathrm{int} \hat{R} \rangle a")
    end

    for a in ax[:, begin]
        a.set_ylabel(L"x / a")
        a.set_ylim(-5, 5)
    end
    for a in ax[end, :]
        a.set_xticks([0, 1, 2, 3])
        a.set_xlabel(L"t / t_\mathrm{MD}")
        a.set_xlim(0, 3.15)
    end

    ax[begin, 1].set_title("finite difference")
    ax[begin, 2].set_title("autodiff")

    mkpath("plot_paper")
    fig.savefig("plot_paper/interacting.pdf");
end

function acceleration()
    fig, ax = subplots(4, 2, figsize=(3.5, 6), sharex=true, sharey=true, layout="constrained")
    cmap = "RdBu_r"
    ρ_absmax = 0.8
    ∇ρ_absmax = 0.8
    τc_absmax = ρ_absmax
    Cacc_absmax = ∇ρ_absmax

    jldopen("data/5WCA_harmonic_stronger_finite_diff.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        ∇ρ = result_histograms.hists["∇ρ"]
        τc = result_histograms.hists["vs[i]{xs[i],H0}"]
        Cacc = result_histograms.hists["{vs[i],H0}"]
        im = ax[1, 1].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        im = ax[2, 1].pcolormesh(ts, xs, ∇ρ, cmap=cmap, vmin=-∇ρ_absmax, vmax=∇ρ_absmax, rasterized=true)
        im = ax[3, 1].pcolormesh(ts, xs, τc, cmap=cmap, vmin=-τc_absmax, vmax=τc_absmax, rasterized=true)
        im = ax[4, 1].pcolormesh(ts, xs, Cacc, cmap=cmap, vmin=-Cacc_absmax, vmax=Cacc_absmax, rasterized=true)
    end

    jldopen("data/5WCA_harmonic_stronger_short1e5.jld2", "r") do file
        result_histograms = file["result_histograms"]
        ts = result_histograms.ts
        xs = result_histograms.rs
        calc_computed!(result_histograms)
        ρ = result_histograms.hists["1"]
        ∇ρ = result_histograms.hists["∇ρ"]
        τc = result_histograms.hists["vs[i]{xs[i],H0}"]
        Cacc = result_histograms.hists["{vs[i],H0}"]
        im = ax[1, 2].pcolormesh(ts, xs, ρ, cmap=cmap, vmin=-ρ_absmax, vmax=ρ_absmax, rasterized=true)
        cbar_ρ = fig.colorbar(im, ax=ax[1, 2])
        cbar_ρ.ax.set_ylim(0, ρ_absmax)
        cbar_ρ.set_label(L"\rho a")
        im = ax[2, 2].pcolormesh(ts, xs, ∇ρ, cmap=cmap, vmin=-∇ρ_absmax, vmax=∇ρ_absmax, rasterized=true)
        cbar_∇ρ = fig.colorbar(im, ax=ax[2, 2])
        cbar_∇ρ.set_label(L"\nabla \rho a^2")
        im = ax[3, 2].pcolormesh(ts, xs, τc, cmap=cmap, vmin=-τc_absmax, vmax=τc_absmax, rasterized=true)
        cbar_τc = fig.colorbar(im, ax=ax[3, 2])
        cbar_τc.ax.set_ylim(0, τc_absmax)
        cbar_τc.set_label(L"-\tau_c a")
        im = ax[4, 2].pcolormesh(ts, xs, Cacc, cmap=cmap, vmin=-Cacc_absmax, vmax=Cacc_absmax, rasterized=true)
        cbar_Cacc = fig.colorbar(im, ax=ax[4, 2])
        cbar_Cacc.set_label(L"C_\mathrm{acc} a^2")
    end

    for a in ax[:, begin]
        a.set_ylim(-5, 5)
        a.set_ylabel(L"x / a")
    end
    for a in ax[end, :]
        a.set_xticks([0, 1, 2, 3])
        a.set_xlabel(L"t / t_\mathrm{MD}")
    end

    ax[begin, 1].set_title(L"10^7" * " trajectories")
    ax[begin, 2].set_title(L"10^5" * " trajectories")

    mkpath("plot_paper")
    fig.savefig("plot_paper/acceleration_splitting.pdf")
end


analytical()
compare()
interacting()
doublewell()
acceleration()
