using PlotlyJS
using ForwardDiff

struct System{T<:Real,F1<:Function,F2<:Function,F3<:Function,F4<:Function}
    N::Int
    ϕ::F1
    Fint::F2
    Vext::F3
    Fext::F4
    Xs::Vector{T}
    Fints::Vector{T}
    Fexts::Vector{T}
    Fs::Vector{T}
    function System(Xs::AbstractVector, ϕ::Function, Fint::Function, Vext::Function, Fext::Function)
        N = length(Xs) ÷ 2
        Fints = similar(Xs[1:N])
        Fexts = similar(Xs[1:N])
        Fs = similar(Xs[1:N])
        new{eltype(Xs),typeof(ϕ),typeof(Fint),typeof(Vext),typeof(Fext)}(N, ϕ, Fint, Vext, Fext, Xs, Fints, Fexts, Fs)
    end
end

function make_dual(system::System)
    Xs = copy(system.Xs)
    jaccfg = ForwardDiff.JacobianConfig(nothing, Xs, ForwardDiff.Chunk{2*system.N}())
    Xsdual = jaccfg.duals
    ForwardDiff.seed!(Xsdual, Xs, jaccfg.seeds)
    System(Xsdual, system.ϕ, system.Fint, system.Vext, system.Fext)
end

struct Histograms
    rs::Vector{Float64}
    ts::Vector{Float64}
    dr::Float64
    dt::Float64
    hists::Dict{String,Array{Float64,2}}
    function Histograms(rs, ts, ∂X0s, with_finite_diff=false)
        dr = rs[2] - rs[1]
        dt = ts[2] - ts[1]
        n_rs = length(rs)
        n_ts = length(ts)
        hists = Dict(
            "1" => zeros(n_rs, n_ts),
            "A" => zeros(n_rs, n_ts),
            "vs[i]" => zeros(n_rs, n_ts),
            "vs[i]²" => zeros(n_rs, n_ts),
            "Fints[i]" => zeros(n_rs, n_ts),
            "Fexts[i]" => zeros(n_rs, n_ts),
            "∇A[i]" => zeros(n_rs, n_ts),
        )
        for X in keys(∂X0s)
            for k in ["kin", "int", "ext"]
                for A in ["", "A"]
                    hists[A*"{"*X*"[i],H"*k*"0}"] = zeros(n_rs, n_ts)
                    hists[A*"vs[i]{"*X*"[i],H"*k*"0}"] = zeros(n_rs, n_ts)
                    if with_finite_diff
                        hists[A*"{"*X*"[i],H"*k*"0}_finite_diff"] = zeros(n_rs, n_ts)
                        hists[A*"vs[i]{"*X*"[i],H"*k*"0}_finite_diff"] = zeros(n_rs, n_ts)
                    end
                end
            end
        end
        new(rs, ts, dr, dt, hists)
    end
end

function calc_forces!(system, t)
    xs = @view system.Xs[begin:system.N]
    system.Fints .= 0
    for i in 1:system.N
        system.Fexts[i] = system.Fext(xs[i], t)
        for j in i+1:system.N
            Fint = system.Fint(xs[i] - xs[j])
            system.Fints[i] += Fint
            system.Fints[j] -= Fint
        end
    end
end

function step!(system, t; dt=0.001, kin=true, int=true, ext=true, p=nothing)
    if !isnothing(p)
        restyle!(p, 1, x=[ForwardDiff.value.(system.Xs[begin:system.N])])
        sleep(dt)
    end
    xs = @view system.Xs[begin:system.N]
    vs = @view system.Xs[system.N+1:end]
    dthalf = 0.5 * dt
    system.Fs .= 0
    if int
        system.Fs .+= system.Fints
    end
    if ext
        system.Fs .+= system.Fexts
    end
    if kin
        xs .+= vs .* dt
    end
    xs .+= system.Fs .* dthalf .* dt
    vs .+= system.Fs .* dthalf
    calc_forces!(system, t + dt)
    system.Fs .= 0
    if int
        system.Fs .+= system.Fints
    end
    if ext
        system.Fs .+= system.Fexts
    end
    vs .+= system.Fs .* dthalf
end

function sample!(histograms, system, t, ∂X0s, X0s, A_funcs=(xs -> 1, xs -> zero(xs)), systems_finite_diff=nothing, Δt_diff=nothing)
    xs = @view system.Xs[begin:system.N]
    vs = @view system.Xs[system.N+1:end]
    t_bin = round(Int, t / histograms.dt) + 1
    A = A_funcs[1](ForwardDiff.value.(xs))
    ∇A = A_funcs[2](ForwardDiff.value.(xs))
    ForwardDiff.extract_jacobian!(Nothing, ∂X0s["xs"], xs, 2 * system.N)
    ForwardDiff.extract_jacobian!(Nothing, ∂X0s["vs"], vs, 2 * system.N)
    # ForwardDiff.extract_jacobian!(Nothing, ∂X0s["Fints"], system.Fints, 2 * system.N)
    # ForwardDiff.extract_jacobian!(Nothing, ∂X0s["Fexts"], system.Fexts, 2 * system.N)
    poisson = Dict()
    if !isnothing(systems_finite_diff)
        Δ = Dict()
        for k in ["kin", "int", "ext"]
            ΔXs = systems_finite_diff[k].Xs .- ForwardDiff.value.(system.Xs)
            Δ["xs,"*k] = ΔXs[begin:system.N]
            Δ["vs,"*k] = ΔXs[system.N+1:end]
            Δ["Fints,"*k] = systems_finite_diff[k].Fints .- ForwardDiff.value.(system.Fints)
            Δ["Fexts,"*k] = systems_finite_diff[k].Fexts .- ForwardDiff.value.(system.Fexts)
        end
    end
    for X in keys(∂X0s)
        for k in ["kin", "int", "ext"]
            poisson[X*","*k] = ∂X0s[X] * X0s[k]
            if !isnothing(systems_finite_diff)
                poisson[X*","*k*"_finite_diff"] = Δ[X*","*k] / Δt_diff
            end
        end
    end
    for i in 1:system.N
        r_bin = ceil(Int, (ForwardDiff.value(xs[i]) + histograms.rs[end] + histograms.dr / 2) / histograms.dr)
        if r_bin < 1 || r_bin > length(histograms.rs)
            continue
        end
        histograms.hists["1"][r_bin, t_bin] += 1
        histograms.hists["A"][r_bin, t_bin] += A
        histograms.hists["∇A[i]"][r_bin, t_bin] += ∇A[i]
        histograms.hists["vs[i]"][r_bin, t_bin] += ForwardDiff.value(vs[i])
        histograms.hists["vs[i]²"][r_bin, t_bin] += ForwardDiff.value(vs[i])^2
        histograms.hists["Fints[i]"][r_bin, t_bin] += ForwardDiff.value(system.Fints[i])
        histograms.hists["Fexts[i]"][r_bin, t_bin] += ForwardDiff.value(system.Fexts[i])
        for X in keys(∂X0s)
            for k in ["kin", "int", "ext"]
                histograms.hists["{"*X*"[i],H"*k*"0}"][r_bin, t_bin] += poisson[X*","*k][i]
                histograms.hists["vs[i]{"*X*"[i],H"*k*"0}"][r_bin, t_bin] += ForwardDiff.value(vs[i]) * poisson[X*","*k][i]
                histograms.hists["A{"*X*"[i],H"*k*"0}"][r_bin, t_bin] += A * poisson[X*","*k][i]
                histograms.hists["Avs[i]{"*X*"[i],H"*k*"0}"][r_bin, t_bin] += A * ForwardDiff.value(vs[i]) * poisson[X*","*k][i]
                if !isnothing(systems_finite_diff)
                    histograms.hists["{"*X*"[i],H"*k*"0}_finite_diff"][r_bin, t_bin] += poisson[X*","*k*"_finite_diff"][i]
                    histograms.hists["vs[i]{"*X*"[i],H"*k*"0}_finite_diff"][r_bin, t_bin] += ForwardDiff.value(vs[i]) * poisson[X*","*k*"_finite_diff"][i]
                    histograms.hists["A{"*X*"[i],H"*k*"0}_finite_diff"][r_bin, t_bin] += A * poisson[X*","*k*"_finite_diff"][i]
                    histograms.hists["Avs[i]{"*X*"[i],H"*k*"0}_finite_diff"][r_bin, t_bin] += A * ForwardDiff.value(vs[i]) * poisson[X*","*k*"_finite_diff"][i]
                end
            end
        end
    end
end

function calc_energy(system, i)
    xi = system.Xs[i]
    vi = system.Xs[system.N+i]
    E = system.Vext(xi)
    E += vi^2 / 2
    for j in 1:system.N
        if i == j
            continue
        end
        xj = system.Xs[j]
        E += system.ϕ(xi - xj)
        if isinf(E)
            break
        end
    end
    E
end

function trial_move!(system::System, T::Number; Δxmax=0.5, Δvmax=0.5)
    i = rand(1:system.N)
    x_before = system.Xs[i]
    v_before = system.Xs[system.N+i]
    Ebefore = calc_energy(system, i)
    system.Xs[i] += Δxmax * (2 * rand() - 1)
    system.Xs[system.N+i] += Δvmax * (2 * rand() - 1)
    Eafter = calc_energy(system, i)
    if rand() > exp(-(Eafter - Ebefore) / T)
        system.Xs[i] = x_before
        system.Xs[system.N+i] = v_before
    end
end

function mc(N::Int, L::Number, T::Number, ϕ::Function, Fint::Function, Vext0::Function, Fext0::Function; equilibration_sweeps=100000, production_sweeps=1000000, sweep_transitions=100)
    Xs = [collect(range(-L/2 + 1.0, L/2 - 1.0, length=N))..., [0 for _ in 1:N]...]
    system = System(Xs, ϕ, Fint, Vext0, Fext0)
    println("MC: generating initial configurations")
    println("Equilibration")
    for i in 1:equilibration_sweeps
        for _ in 1:sweep_transitions
            trial_move!(system, T)
        end
        mod(i, 10000) == 0 && println(i)
    end
    Xs = [zeros(Float64, 2 * N) for _ in 1:production_sweeps]
    println("Production")
    for i in 1:production_sweeps
        for _ in 1:sweep_transitions
            trial_move!(system, T)
        end
        Xs[i] .= system.Xs
        mod(i, 10000) == 0 && println(i)
    end
    Xs
end

function simulate(L::Number, Xs_init::AbstractArray, ϕ::Function, Fint::Function, Vext0::Function, Fext0::Function, Vext1::Function, Fext1::Function; production_time=1.0, dt=0.001, dr=0.05, A_funcs=(xs -> 1, xs -> zero(xs)), run_with_dual=true, finite_diff=false, plot=false)
    Vext(x, t) = t <= 0 ? Vext0(x) : Vext1(x)
    Fext(x, t) = t <= 0 ? Fext0(x) : Fext1(x)
    system = System(Xs_init, ϕ, Fint, Vext, Fext)
    calc_forces!(system, 0)
    rs = (-L/2+dr/2):dr:L/2
    p = nothing
    if plot
        p = PlotlyJS.plot(PlotlyJS.scatter(; x=system.Xs[begin:system.N], y=zeros(system.N), mode="markers", marker_size=20), Layout(xaxis_range=[-L/2, L/2]))
        display(p)
    end
    ts = 0:dt:production_time
    X0s = Dict(
        "kin" => [system.Xs[system.N+1:end]..., zeros(system.N)...],
        "int" => [zeros(system.N)..., system.Fints...],
        "ext" => [zeros(system.N)..., system.Fexts...],
    )
    if finite_diff
        Δt_diff = 0.0001
        systems_finite_diff = Dict(
            "kin" => deepcopy(system),
            "int" => deepcopy(system),
            "ext" => deepcopy(system),
        )
        step!(systems_finite_diff["kin"], -Δt_diff; dt=Δt_diff, kin=true, int=false, ext=false)
        step!(systems_finite_diff["int"], -Δt_diff; dt=Δt_diff, kin=false, int=true, ext=false)
        step!(systems_finite_diff["ext"], -Δt_diff; dt=Δt_diff, kin=false, int=false, ext=true)
    end
    if run_with_dual
        system = make_dual(system)
        calc_forces!(system, 0)
    end
    xs = system.Xs[begin:system.N]
    vs = system.Xs[system.N+1:end]
    ∂X0s = Dict(
        "xs" => similar(xs, ForwardDiff.valtype(eltype(xs)), system.N, 2 * system.N),
        "vs" => similar(vs, ForwardDiff.valtype(eltype(vs)), system.N, 2 * system.N),
        # "Fints" => similar(system.Fints, ForwardDiff.valtype(eltype(system.Fints)), system.N, 2 * system.N),
        # "Fexts" => similar(system.Fexts, ForwardDiff.valtype(eltype(system.Fexts)), system.N, 2 * system.N),
    )
    histograms = Histograms(rs, ts, ∂X0s, finite_diff)
    for t in ts
        if finite_diff
            sample!(histograms, system, t, ∂X0s, X0s, A_funcs, systems_finite_diff, Δt_diff)
        else
            sample!(histograms, system, t, ∂X0s, X0s, A_funcs)
        end
        step!(system, t; dt, p)
        if finite_diff
            for s in values(systems_finite_diff)
                step!(s, t; dt)
            end
        end
    end
    success = true
    for (key, hist) in histograms.hists
        if any(isnan, hist)
            println("$(key) contains NaN")
            success = false
        end
    end
    histograms, success
end
