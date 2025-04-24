using PlotlyJS, JLD2

include("simulation.jl")

function EWCA(r)
    if abs(r) > 2^(1 / 6)
        return zero(r)
    end
    rinv = 1 / r
    4 * (rinv^12 - rinv^6) + 1
end

function FWCA(r)
    if abs(r) > 2^(1 / 6)
        return zero(r)
    end
    rinv = 1 / r
    24 * (2 * rinv^12 - rinv^6) * rinv
end

Eharmonic(x; k=1.0) = 0.5 * k * x^2
Fharmonic(x; k=1.0) = -k * x

Edoublewell(x; ϵ=1.0, Δ=2.5) = ϵ * (x^2 − Δ^2)^2 / Δ^4
Fdoublewell(x; ϵ=1.0, Δ=2.5) = -4 * ϵ * x * (x^2 - Δ^2) / Δ^4

function norm_histograms!(result_histograms, acc_histograms, acc)
    norm = 1 / acc / result_histograms.dr
    for key in keys(acc_histograms.hists)
        result_histograms.hists[key] .= acc_histograms.hists[key] .* norm
    end
end

function plot_histograms(result_histograms::Histograms; key="1", p=nothing)
    if isnothing(p)
        p = PlotlyJS.plot(PlotlyJS.heatmap(; x=result_histograms.ts, y=result_histograms.rs, z=result_histograms.hists[key], zmid=0, colorscale="RdBu"), Layout(xaxis_title="t", yaxis_title="x"))
        display(p)
    else
        restyle!(p, 1, z=[result_histograms.hists[key]])
    end
    p
end

function plot_histograms(filename="data/current.jld2"; key="1", p=nothing)
    jldopen(filename, "r+") do file
        result_histograms = file["result_histograms"]
        return plot_histograms(result_histograms; key, p)
    end
end

function run_simulations()
    N = 5
    L = 20.0
    T = 1.0

    Vext0(x) = Eharmonic(x; k=0.25)
    Fext0(x) = Fharmonic(x; k=0.25)

    Vext1(x) = Eharmonic(x; k=1)
    Fext1(x) = Fharmonic(x; k=1)

    ϕ(r) = EWCA(r)
    Fint(r) = FWCA(r)

    # A_funcs = (xs -> 1, xs -> zero(xs))
    A_funcs = (xs -> sum(xs) / length(xs), xs -> ones(length(xs)) ./ length(xs))

    filename = "data/current.jld2"

    n_Xs = 100000
    production_time = 3.15
    dt = 5e-3
    finite_diff = true

    Xs_inits = mc(N, L, T, ϕ, Fint, Vext0, Fext0; production_sweeps=n_Xs)

    histograms, success = simulate(L, Xs_inits[1], ϕ, Fint, Vext0, Fext0, Vext1, Fext1; production_time, dt, A_funcs, finite_diff)
    for key in keys(histograms.hists)
        histograms.hists[key] .= 0.0
    end

    result_histograms = deepcopy(histograms)
    acc_histograms = deepcopy(histograms)
    accs = 0
    p = nothing
    for runid in eachindex(Xs_inits)
        Xs_init = Xs_inits[runid]
        println("Doing run $(runid)")
        histograms, success = simulate(L, Xs_init, ϕ, Fint, Vext0, Fext0, Vext1, Fext1; production_time, dt, A_funcs, finite_diff)
        m = 0.0
        if success
            for hist in values(histograms.hists)
                m = max(m, maximum(abs, hist))
            end
            if m > 1e6
                println("Detected large histogram entry ($(m)), discarding this run")
                continue
            end
            accs += 1
            for (key, hist) in histograms.hists
                acc_histograms.hists[key] .+= hist
            end
        end
        if mod(runid, 100) == 1
            norm_histograms!(result_histograms, acc_histograms, accs)
            jldsave(filename; result_histograms=result_histograms, acc_histograms=acc_histograms, acc=accs)
            p = plot_histograms(filename; key="1", p=p)
        end
    end
end

run_simulations()
