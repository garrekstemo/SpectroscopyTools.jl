"""
Curve fitting routines for spectroscopy data.

Model functions come from CurveFitModels.jl — this file contains
automatic fitting routines that use those models.
"""

# Internal IRF functions

function _erfc(x)
    if x < 0
        return 2.0 - _erfc(-x)
    end
    t = 1.0 / (1.0 + 0.3275911 * x)
    tau = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
          t * (-1.453152027 + t * 1.061405429))))
    return tau * exp(-x^2)
end

function _exp_decay_irf_conv(t, A, tau, t0, sigma)
    t_shifted = t - t0
    arg_exp = sigma^2 / (2 * tau^2) - t_shifted / tau
    arg_erfc = (sigma^2 / tau - t_shifted) / (sigma * sqrt(2))
    return (A / 2) * exp(arg_exp) * _erfc(arg_erfc)
end

# Internal helpers for fitting

# From pre-computed residual sum of squares (e.g., rss(sol) from CurveFit)
function _rsquared(y_data, ss_res::Real)
    ss_tot = sum((y_data .- mean(y_data)).^2)
    return ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0
end

# From fitted values (e.g., for per-trace R² in global fitting)
_rsquared(y_data, y_fit::AbstractVector) = _rsquared(y_data, sum((y_data .- y_fit).^2))

function _detect_signal_type(signal)
    max_val = maximum(signal)
    min_val = minimum(signal)
    if abs(max_val) >= abs(min_val)
        return :esa, max_val, argmax(signal)
    else
        return :gsb, min_val, argmin(signal)
    end
end

"""
    fit_decay_irf(t, signal; sigma_init=5.0) -> ExpDecayFit

Fit an exponential decay convolved with a Gaussian IRF to pump-probe data.
"""
function fit_decay_irf(t::AbstractVector{<:Real}, signal::AbstractVector{<:Real};
                       sigma_init::Real=5.0)
    signal_type, peak_val, peak_idx = _detect_signal_type(signal)

    t0_init = t[peak_idx]
    n_edge = max(5, length(signal) ÷ 20)
    offset_init = mean(signal[1:n_edge])
    A_init = peak_val - offset_init

    half_val = (peak_val + offset_init) / 2
    half_idx = peak_idx
    for i in peak_idx:length(signal)
        if signal_type == :esa && signal[i] < half_val
            half_idx = i
            break
        elseif signal_type == :gsb && signal[i] > half_val
            half_idx = i
            break
        end
    end
    tau_init = max(abs(t[half_idx] - t[peak_idx]) / log(2), 1.0)

    function model(p, t_vec)
        A, tau, t0, sigma, offset = p
        tau = abs(tau)
        sigma = abs(sigma)
        return [_exp_decay_irf_conv(ti, A, tau, t0, sigma) + offset for ti in t_vec]
    end

    p0 = [A_init, tau_init, t0_init, sigma_init, offset_init]

    prob = NonlinearCurveFitProblem(model, p0, t, signal)
    sol = solve(prob)

    A, tau, t0, sigma, offset = coef(sol)
    tau = abs(tau)
    sigma = abs(sigma)

    return ExpDecayFit(A, tau, t0, sigma, offset, signal_type,
                       residuals(sol), _rsquared(signal, rss(sol)))
end

# Pulse width estimation

const FWHM_FACTOR = 2 * sqrt(2 * log(2))  # ≈ 2.355

"""
    irf_fwhm(sigma)

Convert IRF standard deviation σ to full-width at half-maximum (FWHM).
"""
irf_fwhm(sigma) = FWHM_FACTOR * sigma

"""
    pulse_fwhm(sigma_irf)

Estimate individual pulse FWHM from fitted IRF width, assuming identical
Gaussian pump and probe pulses.
"""
pulse_fwhm(sigma_irf) = FWHM_FACTOR * sigma_irf / sqrt(2)


# =============================================================================
# Multi-exponential IRF convolution
# =============================================================================

function _multiexp_irf_conv(t, taus, amplitudes, t0, sigma, offset)
    result = offset
    for i in eachindex(taus)
        result += _exp_decay_irf_conv(t, amplitudes[i], taus[i], t0, sigma)
    end
    return result
end

# =============================================================================
# Unified TA API: fit_exp_decay
# =============================================================================

"""
    fit_exp_decay(trace::TATrace; n_exp=1, irf=false, irf_width=0.15, t_start=0.0, t_range=nothing)

Fit exponential decay to a transient absorption trace.

# Arguments
- `trace`: TATrace
- `n_exp`: Number of exponential components (default 1)
- `irf`: Include IRF convolution (default false)
- `irf_width`: Initial guess for IRF σ in ps (default 0.15)
- `t_start`: Start time for fitting when irf=false (default 0.0)
- `t_range`: Optional (t_min, t_max) to restrict fit region

# Returns
- `n_exp=1`: `ExpDecayFit`
- `n_exp>1`: `MultiexpDecayFit`
"""
function fit_exp_decay(trace::TATrace; n_exp::Int=1, irf::Bool=false, irf_width::Float64=0.15,
                       t_start::Float64=0.0, t_range=nothing)
    @assert n_exp >= 1 "n_exp must be at least 1"

    if n_exp > 1
        return _fit_multiexp_decay(trace; n_exp=n_exp, irf=irf, irf_width=irf_width,
                                   t_start=t_start, t_range=t_range)
    end

    t = trace.time
    signal = trace.signal

    if !isnothing(t_range)
        t_min, t_max = t_range
        mask = (t .>= t_min) .& (t .<= t_max)
        t = t[mask]
        signal = signal[mask]
    end

    if irf
        return fit_decay_irf(t, signal; sigma_init=irf_width)
    else
        mask = t .>= t_start
        t_fit = t[mask]
        signal_fit = signal[mask]

        signal_type = first(_detect_signal_type(signal))
        peak_val = signal_type == :esa ? maximum(signal_fit) : minimum(signal_fit)

        n_end = max(1, min(10, length(signal_fit) ÷ 4))
        offset0 = mean(signal_fit[end-n_end+1:end])

        A0 = peak_val - offset0
        tau0 = (t_fit[end] - t_fit[1]) / 3.0

        function model(p, t_vec)
            A, tau, offset = p
            return @. A * exp(-(t_vec - t_start) / abs(tau)) + offset
        end

        prob = NonlinearCurveFitProblem(model, [A0, tau0, offset0], t_fit, signal_fit)
        sol = solve(prob)

        A, tau, offset = coef(sol)
        tau = abs(tau)

        return ExpDecayFit(
            A, tau, t_start, NaN, offset,
            signal_type, residuals(sol), _rsquared(signal_fit, rss(sol))
        )
    end
end

# =============================================================================
# Multi-exponential fitting (internal)
# =============================================================================

function _fit_multiexp_decay(trace::TATrace; n_exp::Int, irf::Bool, irf_width::Float64,
                             t_start::Float64, t_range)
    t = trace.time
    signal = trace.signal

    if !isnothing(t_range)
        t_min, t_max = t_range
        mask = (t .>= t_min) .& (t .<= t_max)
        t = t[mask]
        signal = signal[mask]
    end

    if irf
        t_fit = t
        signal_fit = signal
    else
        mask = t .>= t_start
        t_fit = t[mask]
        signal_fit = signal[mask]
    end

    signal_type = first(_detect_signal_type(signal))
    peak_val = signal_type == :esa ? maximum(signal_fit) : minimum(signal_fit)

    n_end = max(1, min(10, length(signal_fit) ÷ 4))
    offset0 = mean(signal_fit[end-n_end+1:end])

    half_val = (peak_val + offset0) / 2
    half_idx = findfirst(i -> begin
        if signal_type == :esa
            signal_fit[i] <= half_val
        else
            signal_fit[i] >= half_val
        end
    end, eachindex(signal_fit))
    tau_est = isnothing(half_idx) ? (t_fit[end] - t_fit[1]) / 3 : t_fit[half_idx] - t_fit[1]
    tau_est = max(tau_est, 0.1)

    taus_init = Float64[]
    if n_exp == 1
        push!(taus_init, tau_est)
    elseif n_exp == 2
        push!(taus_init, tau_est / 3.0)
        push!(taus_init, tau_est * 2.0)
    else
        log_min = log10(tau_est / 10)
        log_max = log10(tau_est * 10)
        for i in 1:n_exp
            push!(taus_init, 10^(log_min + (log_max - log_min) * (i - 1) / (n_exp - 1)))
        end
    end

    total_amp = peak_val - offset0
    amps_init = fill(total_amp / n_exp, n_exp)

    if irf
        p0 = vcat(taus_init, amps_init, [0.0, irf_width, offset0])

        function multiexp_irf_model(p, t_vec)
            taus = abs.(p[1:n_exp])
            amps = p[n_exp+1:2*n_exp]
            t0 = p[2*n_exp+1]
            sigma = abs(p[2*n_exp+2])
            offset = p[2*n_exp+3]
            return [_multiexp_irf_conv(ti, taus, amps, t0, sigma, offset) for ti in t_vec]
        end

        prob = NonlinearCurveFitProblem(multiexp_irf_model, p0, t_fit, signal_fit)
        sol = solve(prob)
        p_opt = coef(sol)

        taus_fit = abs.(p_opt[1:n_exp])
        amps_fit = p_opt[n_exp+1:2*n_exp]
        t0 = p_opt[2*n_exp+1]
        sigma = abs(p_opt[2*n_exp+2])
        offset = p_opt[2*n_exp+3]
    else
        model = n_exponentials(n_exp)

        p0 = Float64[]
        for i in 1:n_exp
            push!(p0, amps_init[i])
            push!(p0, taus_init[i])
        end
        push!(p0, offset0)

        t_shifted = t_fit .- t_fit[1]

        prob = NonlinearCurveFitProblem(model, p0, t_shifted, signal_fit)
        sol = solve(prob)
        p_opt = coef(sol)

        amps_fit = Float64[]
        taus_fit = Float64[]
        for i in 1:n_exp
            push!(amps_fit, p_opt[2*i - 1])
            push!(taus_fit, abs(p_opt[2*i]))
        end
        offset = p_opt[end]
        t0 = t_fit[1]
        sigma = NaN
    end

    sort_idx = sortperm(taus_fit)
    taus_sorted = taus_fit[sort_idx]
    amps_sorted = amps_fit[sort_idx]

    rsquared = _rsquared(signal_fit, rss(sol))

    max_reasonable_tau = 10 * (t_fit[end] - t_fit[1])
    if any(τ -> τ > max_reasonable_tau, taus_sorted)
        @warn "Multi-exponential fit may have failed — time constants unreasonably large. Consider fewer components."
    end

    return MultiexpDecayFit(
        taus_sorted, amps_sorted,
        t0, sigma, offset,
        signal_type, residuals(sol), rsquared
    )
end

# =============================================================================
# Global fitting: shared time constant across multiple traces
# =============================================================================

"""
    fit_global(traces::Vector{TATrace}; irf_width=0.15, labels=nothing) -> GlobalFitResult

Fit multiple traces simultaneously with a shared time constant τ.
"""
function fit_global(traces::Vector{TATrace}; irf_width::Float64=0.15, labels=nothing)
    n_traces = length(traces)
    @assert n_traces >= 2 "Need at least 2 traces for global fitting"

    if isnothing(labels)
        labels = ["Trace $i" for i in 1:n_traces]
    end

    individual_fits = [fit_exp_decay(tr; irf=true, irf_width=irf_width) for tr in traces]

    tau_init = mean([f.tau for f in individual_fits])
    sigma_init = mean([f.sigma for f in individual_fits])
    t0_init = mean([f.t0 for f in individual_fits])

    p0 = Float64[tau_init, sigma_init, t0_init]
    for f in individual_fits
        push!(p0, f.amplitude)
        push!(p0, f.offset)
    end

    total_len = sum(length(tr.time) for tr in traces)

    function global_model(p, dummy_x)
        tau, sigma, t0 = abs(p[1]), abs(p[2]), p[3]

        y_pred = similar(p, total_len)
        idx = 1
        for i in 1:n_traces
            A = p[3 + 2*(i-1) + 1]
            offset = p[3 + 2*(i-1) + 2]
            t_vec = traces[i].time

            for t in t_vec
                y_pred[idx] = _exp_decay_irf_conv(t, A, tau, t0, sigma) + offset
                idx += 1
            end
        end
        return y_pred
    end

    x_all = Float64[]
    y_all = Float64[]
    for tr in traces
        append!(x_all, tr.time)
        append!(y_all, tr.signal)
    end

    prob = NonlinearCurveFitProblem(global_model, p0, x_all, y_all)
    sol = solve(prob)
    p_opt = coef(sol)

    tau = abs(p_opt[1])
    sigma = abs(p_opt[2])
    t0 = p_opt[3]

    amplitudes = [p_opt[3 + 2*(i-1) + 1] for i in 1:n_traces]
    offsets = [p_opt[3 + 2*(i-1) + 2] for i in 1:n_traces]

    residuals_vec = Vector{Vector{Float64}}(undef, n_traces)
    rsquared_individual = zeros(n_traces)

    resid_all = residuals(sol)
    y_pred_all = fitted(sol)
    idx = 1
    for i in 1:n_traces
        n_pts = length(traces[i].time)
        residuals_vec[i] = resid_all[idx:idx+n_pts-1]
        rsquared_individual[i] = _rsquared(traces[i].signal, y_pred_all[idx:idx+n_pts-1])
        idx += n_pts
    end

    rsquared_global = _rsquared(y_all, rss(sol))

    return GlobalFitResult(
        tau, sigma, t0,
        amplitudes, offsets,
        labels,
        rsquared_global, rsquared_individual,
        residuals_vec
    )
end

# =============================================================================
# Predict functions for fit results
# =============================================================================

function predict(fit::ExpDecayFit, time::AbstractVector)
    if isnan(fit.sigma)
        return [t >= fit.t0 ? fit.amplitude * exp(-(t - fit.t0) / fit.tau) + fit.offset : fit.offset
                for t in time]
    else
        return [_exp_decay_irf_conv(t, fit.amplitude, fit.tau, fit.t0, fit.sigma) + fit.offset
                for t in time]
    end
end

predict(fit::ExpDecayFit, trace::TATrace) = predict(fit, trace.time)

function predict(fit::GlobalFitResult, traces::Vector{TATrace})
    n = length(traces)
    curves = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        A = fit.amplitudes[i]
        offset = fit.offsets[i]
        curves[i] = [_exp_decay_irf_conv(t, A, fit.tau, fit.t0, fit.sigma) + offset
                     for t in traces[i].time]
    end
    return curves
end


function predict(fit::MultiexpDecayFit, time::AbstractVector)
    if isnan(fit.sigma)
        return [begin
            val = fit.offset
            for i in eachindex(fit.taus)
                if t >= fit.t0
                    val += fit.amplitudes[i] * exp(-(t - fit.t0) / fit.taus[i])
                end
            end
            val
        end for t in time]
    else
        return [_multiexp_irf_conv(t, fit.taus, fit.amplitudes, fit.t0, fit.sigma, fit.offset)
                for t in time]
    end
end

predict(fit::MultiexpDecayFit, trace::TATrace) = predict(fit, trace.time)

# =============================================================================
# TA Spectrum Fitting (generalized N-peak model)
# =============================================================================

const _PEAK_SIGNS = Dict(:esa => 1, :gsb => -1, :se => -1, :positive => 1, :negative => -1)

function _build_ta_model(fns, signs, npps, fit_offset)
    function model(p, x)
        y = similar(p, length(x))
        fill!(y, zero(eltype(p)))

        p_idx = 1
        for i in eachindex(fns)
            npp = npps[i]
            peak_p = p[p_idx:p_idx+npp-1]
            y = y .+ signs[i] .* fns[i](peak_p, x)
            p_idx += npp
        end

        if fit_offset
            y = y .+ p[p_idx]
        end

        return y
    end
    return model
end

function _ta_initial_guesses(ν, y, peak_specs, signs, npps, fit_offset)
    p0 = Float64[]

    n_pos = count(s -> s > 0, signs)
    n_neg = count(s -> s < 0, signs)

    # Detect positive peaks (ESA-type) using find_peaks on the raw signal
    pos_peaks = PeakInfo[]
    if n_pos > 0
        detected = find_peaks(ν, y; min_prominence=0.01)
        if length(detected) < n_pos
            detected = _synthesize_peak_guesses(ν, y, n_pos, detected)
        elseif length(detected) > n_pos
            sort!(detected, by=p -> p.prominence, rev=true)
            detected = detected[1:n_pos]
            sort!(detected, by=p -> p.position)
        end
        pos_peaks = detected
    end

    # Detect negative peaks (GSB/SE-type) by inverting the signal
    neg_peaks = PeakInfo[]
    if n_neg > 0
        detected = find_peaks(ν, -y; min_prominence=0.01)
        if length(detected) < n_neg
            detected = _synthesize_peak_guesses(ν, -y, n_neg, detected)
        elseif length(detected) > n_neg
            sort!(detected, by=p -> p.prominence, rev=true)
            detected = detected[1:n_neg]
            sort!(detected, by=p -> p.position)
        end
        neg_peaks = detected
    end

    pos_idx = 0
    neg_idx = 0

    for (i, (label, fn)) in enumerate(peak_specs)
        if signs[i] > 0
            pos_idx += 1
            pk = pos_peaks[pos_idx]
        else
            neg_idx += 1
            pk = neg_peaks[neg_idx]
        end

        push!(p0, pk.prominence)
        push!(p0, pk.position)
        push!(p0, _width_guess(pk.width, fn))

        if npps[i] >= 4 && fn === pseudo_voigt
            push!(p0, 0.5)
        end
    end

    if fit_offset
        push!(p0, 0.0)
    end

    return p0
end

"""
    fit_ta_spectrum(spec::TASpectrum; kwargs...) -> TASpectrumFit

Fit a transient absorption spectrum with N peaks of arbitrary lineshape.

Uses `find_peaks` to automatically detect initial peak positions from the data,
so multiple well-separated peaks of the same type (e.g., three GSB peaks for
W(CO)₆) are initialized correctly.

# Keywords
- `peaks=[:esa, :gsb]` — Peak types. Each element is either a `Symbol`
  (`:esa`, `:gsb`, `:se`, `:positive`, `:negative`) or a `(Symbol, Function)`
  tuple specifying label and lineshape model.
- `model=gaussian` — Default lineshape for peaks specified as symbols only.
- `region=nothing` — Optional `(x_min, x_max)` fitting region.
- `fit_offset=false` — Whether to fit a constant offset.
- `p0=nothing` — Manual initial parameter vector. Overrides automatic detection.

# Peak signs
- `:esa`, `:positive` → +1 (positive ΔA)
- `:gsb`, `:se`, `:negative` → -1 (negative ΔA)

# Examples
```julia
# Default: 1 Gaussian ESA + 1 Gaussian GSB
result = fit_ta_spectrum(spec)

# Three GSB peaks (e.g., W(CO)₆ carbonyl stretches)
result = fit_ta_spectrum(spec; peaks=[:esa, :esa, :esa, :gsb, :gsb, :gsb])

# Per-peak lineshapes
result = fit_ta_spectrum(spec; peaks=[(:esa, lorentzian), (:gsb, gaussian)])

# Access results
result[:esa].center      # first ESA peak
result[2].center         # second peak by index
anharmonicity(result)    # GSB - ESA center (only if exactly 1 of each)
predict(result, ν)       # full fitted curve
predict_peak(result, 1)  # single peak contribution
```
"""
function fit_ta_spectrum(spec::TASpectrum;
                         peaks::AbstractVector=[:esa, :gsb],
                         model::Function=gaussian,
                         region=nothing,
                         fit_offset::Bool=false,
                         p0::Union{Nothing, AbstractVector}=nothing)

    if isnothing(region)
        ν = collect(Float64, spec.wavenumber)
        y = collect(Float64, spec.signal)
    else
        mask = (spec.wavenumber .>= region[1]) .& (spec.wavenumber .<= region[2])
        ν = collect(Float64, spec.wavenumber[mask])
        y = collect(Float64, spec.signal[mask])
    end

    peak_specs = if eltype(peaks) <: Symbol
        [(label, model) for label in peaks]
    else
        [(label, fn) for (label, fn) in peaks]
    end

    n_peaks = length(peak_specs)

    signs = Int[]
    for (label, _) in peak_specs
        haskey(_PEAK_SIGNS, label) || error("Unknown peak type :$label. Use :esa, :gsb, :se, :positive, or :negative.")
        push!(signs, _PEAK_SIGNS[label])
    end

    fns = [fn for (_, fn) in peak_specs]
    npps = [_n_peak_params(fn) for fn in fns]

    p0_use = if isnothing(p0)
        _ta_initial_guesses(ν, y, peak_specs, signs, npps, fit_offset)
    else
        collect(Float64, p0)
    end

    composite = _build_ta_model(fns, signs, npps, fit_offset)
    prob = NonlinearCurveFitProblem(composite, p0_use, ν, y)
    sol = solve(prob)
    p_opt = coef(sol)

    ta_peaks = TAPeak[]
    idx = 1
    for i in 1:n_peaks
        label, fn = peak_specs[i]
        npp = npps[i]
        amp = abs(p_opt[idx])
        center = p_opt[idx + 1]
        width = abs(p_opt[idx + 2])
        idx += npp

        push!(ta_peaks, TAPeak(label, _model_name(fn), center, width, amp))
    end

    offset_val = fit_offset ? p_opt[end] : 0.0

    return TASpectrumFit(
        ta_peaks, offset_val,
        _rsquared(y, rss(sol)), residuals(sol),
        collect(p_opt), fns, signs, npps, fit_offset, ν
    )
end

function predict(fit::TASpectrumFit)
    return predict(fit, fit._x)
end

function predict(fit::TASpectrumFit, x::AbstractVector)
    x_f = collect(Float64, x)
    y = zeros(length(x_f))
    idx = 1
    for i in eachindex(fit._peak_fns)
        npp = fit._peak_npp[i]
        peak_p = fit._coef[idx:idx+npp-1]
        y .+= fit._peak_signs[i] .* fit._peak_fns[i](peak_p, x_f)
        idx += npp
    end
    if fit._fit_offset
        y .+= fit._coef[end]
    end
    return y
end

predict(fit::TASpectrumFit, spec::TASpectrum) = predict(fit, spec.wavenumber)

function predict_peak(fit::TASpectrumFit, i::Int)
    return predict_peak(fit, i, fit._x)
end

function predict_peak(fit::TASpectrumFit, i::Int, x::AbstractVector)
    1 <= i <= length(fit.peaks) || throw(BoundsError(fit.peaks, i))
    x_f = collect(Float64, x)
    idx = 1
    for j in 1:i-1
        idx += fit._peak_npp[j]
    end
    npp = fit._peak_npp[i]
    peak_p = fit._coef[idx:idx+npp-1]
    return fit._peak_signs[i] .* fit._peak_fns[i](peak_p, x_f)
end
