using Statistics

function finite_diff(a; dx)
    dv = diff(a; dims=1) / 2
    v = vcat(dv[1:1, :], dv)
    v .+= vcat(dv, dv[end:end, :])
    v / dx
end

function calc_computed!(histogram)
    histogram.hists["∇ρ"] = finite_diff(histogram.hists["1"]; dx=histogram.dr)
    histogram.hists["∇A"] = finite_diff(histogram.hists["A"]; dx=histogram.dr)
    τ = histogram.hists["vs[i]²"]
    ∇τ = finite_diff(τ; dx=histogram.dr)
    histogram.hists["∇τ"] = ∇τ
    histogram.hists["∇A-∇A[i]"] = histogram.hists["∇A"] - histogram.hists["∇A[i]"]
    histogram.hists["Fs[i]"] = histogram.hists["Fints[i]"] + histogram.hists["Fexts[i]"]
    for A in ["", "A"]
        histogram.hists[A*"{xs[i],H0}"] = histogram.hists[A*"{xs[i],Hkin0}"] + histogram.hists[A*"{xs[i],Hint0}"] + histogram.hists[A*"{xs[i],Hext0}"]
        histogram.hists[A*"vs[i]{xs[i],H0}"] = histogram.hists[A*"vs[i]{xs[i],Hkin0}"] + histogram.hists[A*"vs[i]{xs[i],Hint0}"] + histogram.hists[A*"vs[i]{xs[i],Hext0}"]
        histogram.hists[A*"{vs[i],H0}"] = histogram.hists[A*"{vs[i],Hkin0}"] + histogram.hists[A*"{vs[i],Hint0}"] + histogram.hists[A*"{vs[i],Hext0}"]
        histogram.hists[A*"Ckin"] = histogram.hists[A*"{vs[i],Hkin0}"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hkin0}"]; dx=histogram.dr)
        histogram.hists[A*"Cint"] = histogram.hists[A*"{vs[i],Hint0}"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hint0}"]; dx=histogram.dr)
        histogram.hists[A*"Cext"] = histogram.hists[A*"{vs[i],Hext0}"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hext0}"]; dx=histogram.dr)
        histogram.hists[A*"C"] = histogram.hists[A*"Ckin"] + histogram.hists[A*"Cint"] + histogram.hists[A*"Cext"]
        if "{vs[i],Hkin0}_finite_diff" in keys(histogram.hists)
            histogram.hists[A*"Ckin_finite_diff"] = histogram.hists[A*"{vs[i],Hkin0}_finite_diff"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hkin0}_finite_diff"]; dx=histogram.dr)
            histogram.hists[A*"Cint_finite_diff"] = histogram.hists[A*"{vs[i],Hint0}_finite_diff"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hint0}_finite_diff"]; dx=histogram.dr)
            histogram.hists[A*"Cext_finite_diff"] = histogram.hists[A*"{vs[i],Hext0}_finite_diff"] - finite_diff(histogram.hists[A*"vs[i]{xs[i],Hext0}_finite_diff"]; dx=histogram.dr)
        end
    end
end

function coarsen_t(ts, A; factor=50)
    A_coarse = hcat([mean(A[:, 1+bin-factor÷2:bin-factor÷2+factor]; dims=2) for bin in factor÷2:factor:size(A)[2]]...)
    return ts[1+factor÷2:factor:end], A_coarse
end
