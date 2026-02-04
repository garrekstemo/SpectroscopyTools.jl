"""
Curve fitting routines for spectroscopy data.

Model functions come from CurveFitModels.jl — this file contains
automatic fitting routines that use those models.
"""

using CurveFit
using CurveFitModels

"""
    fit_decay(signal::Vector{Float64}; t0_threshold::Float64=0.1) -> ExpDecayFit

Fit an exponential decay to a pump-probe signal.

Automatically detects time zero and signal type (ESA vs GSB).
"""
function fit_decay(signal::Vector{Float64}; t0_threshold::Float64=0.1)
    max_val = maximum(signal)
    min_val = minimum(signal)

    if abs(max_val) > abs(min_val)
        signal_type = :esa
        peak_val = max_val
        peak_idx = argmax(signal)
    else
        signal_type = :gsb
        peak_val = min_val
        peak_idx = argmin(signal)
    end

    threshold = t0_threshold * peak_val
    t0 = 1
    for i in 1:peak_idx
        if signal_type == :esa && signal[i] > threshold
            t0 = i
            break
        elseif signal_type == :gsb && signal[i] < threshold
            t0 = i
            break
        end
    end

    decay_signal = signal[peak_idx:end]
    t_decay = Float64.(collect(0:length(decay_signal)-1))

    n_end = min(10, length(decay_signal) ÷ 4)
    offset0 = mean(decay_signal[end-n_end+1:end])

    A0 = peak_val - offset0
    tau0 = length(decay_signal) / 3.0

    prob = NonlinearCurveFitProblem(single_exponential, [A0, tau0, offset0], t_decay, decay_signal)
    sol = solve(prob)

    A, tau, offset = coef(sol)

    fitted_vals = single_exponential([A, tau, offset], t_decay)
    resid = decay_signal .- fitted_vals
    ss_res = sum(resid.^2)
    ss_tot = sum((decay_signal .- mean(decay_signal)).^2)
    rsquared = 1 - ss_res / ss_tot

    return ExpDecayFit(A, tau, offset, t0, signal_type, resid, rsquared)
end

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

"""
    fit_decay_irf(t::Vector{Float64}, signal::Vector{Float64};
                  sigma_init::Float64=5.0) -> ExpDecayIRFFit

Fit an exponential decay convolved with a Gaussian IRF to pump-probe data.
"""
function fit_decay_irf(t::Vector{Float64}, signal::Vector{Float64};
                       sigma_init::Float64=5.0)
    max_val = maximum(signal)
    min_val = minimum(signal)

    if abs(max_val) > abs(min_val)
        signal_type = :esa
        peak_val = max_val
        peak_idx = argmax(signal)
    else
        signal_type = :gsb
        peak_val = min_val
        peak_idx = argmin(signal)
    end

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

    fitted_vals = model([A, tau, t0, sigma, offset], t)
    resid = signal .- fitted_vals
    ss_res = sum(resid.^2)
    ss_tot = sum((signal .- mean(signal)).^2)
    rsquared = 1 - ss_res / ss_tot

    return ExpDecayIRFFit(A, tau, t0, sigma, offset, signal_type, resid, rsquared)
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
    fit_exp_decay(trace::TATrace; n_exp=1, irf=true, irf_width=0.15, t_start=0.0, t_range=nothing)

Fit exponential decay to a transient absorption trace.

# Arguments
- `trace`: TATrace
- `n_exp`: Number of exponential components (default 1)
- `irf`: Include IRF convolution (default true)
- `irf_width`: Initial guess for IRF σ in ps (default 0.15)
- `t_start`: Start time for fitting when irf=false (default 0.0)
- `t_range`: Optional (t_min, t_max) to restrict fit region

# Returns
- `n_exp=1`: `ExpDecayIRFFit`
- `n_exp>1`: `MultiexpDecayFit`
"""
function fit_exp_decay(trace::TATrace; n_exp::Int=1, irf::Bool=true, irf_width::Float64=0.15,
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

        if abs(maximum(signal)) >= abs(minimum(signal))
            signal_type = :esa
            peak_val = maximum(signal_fit)
        else
            signal_type = :gsb
            peak_val = minimum(signal_fit)
        end

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

        fitted_vals = model([A, tau, offset], t_fit)
        resid = signal_fit .- fitted_vals
        ss_res = sum(resid.^2)
        ss_tot = sum((signal_fit .- mean(signal_fit)).^2)
        rsquared = 1 - ss_res / ss_tot

        return ExpDecayIRFFit(
            A, tau, t_start, NaN, offset,
            signal_type, resid, rsquared
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

    if abs(maximum(signal)) >= abs(minimum(signal))
        signal_type = :esa
        peak_val = maximum(signal_fit)
    else
        signal_type = :gsb
        peak_val = minimum(signal_fit)
    end

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

        fitted_vals = multiexp_irf_model(p_opt, t_fit)
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

        fitted_vals = model(p_opt, t_shifted)
    end

    sort_idx = sortperm(taus_fit)
    taus_sorted = taus_fit[sort_idx]
    amps_sorted = amps_fit[sort_idx]

    resid = signal_fit .- fitted_vals
    ss_res = sum(resid.^2)
    ss_tot = sum((signal_fit .- mean(signal_fit)).^2)
    rsquared = 1 - ss_res / ss_tot

    return MultiexpDecayFit(
        taus_sorted, amps_sorted,
        t0, sigma, offset,
        signal_type, resid, rsquared
    )
end

# =============================================================================
# Biexponential fitting
# =============================================================================

function _biexp_irf_conv(t, A1, tau1, A2, tau2, t0, sigma, offset)
    return _exp_decay_irf_conv(t, A1, tau1, t0, sigma) +
           _exp_decay_irf_conv(t, A2, tau2, t0, sigma) + offset
end

"""
    fit_biexp_decay(trace::TATrace; irf=true, irf_width=0.15, t_range=nothing) -> BiexpDecayFit

Fit biexponential decay to a transient absorption trace.
"""
function fit_biexp_decay(trace::TATrace; irf::Bool=true, irf_width::Float64=0.15,
                         t_start::Float64=0.0, t_range=nothing)
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

    if abs(maximum(signal)) >= abs(minimum(signal))
        signal_type = :esa
        peak_val = maximum(signal_fit)
    else
        signal_type = :gsb
        peak_val = minimum(signal_fit)
    end

    n_end = max(1, min(10, length(signal_fit) ÷ 4))
    offset0 = mean(signal_fit[end-n_end+1:end])

    t_shifted = t_fit .- t_start

    half_val = (peak_val + offset0) / 2
    half_idx = findfirst(i -> begin
        if signal_type == :esa
            signal_fit[i] <= half_val
        else
            signal_fit[i] >= half_val
        end
    end, eachindex(signal_fit))
    tau_est = isnothing(half_idx) ? t_shifted[end] / 3 : t_shifted[half_idx]
    tau_est = max(tau_est, 0.1)

    total_amp = peak_val - offset0
    A1_0 = 0.3 * total_amp
    A2_0 = 0.7 * total_amp
    tau1_0 = tau_est / 3.0
    tau2_0 = tau_est * 2.0

    if irf
        function biexp_irf_model(p, t_vec)
            A1, tau1, A2, tau2, t0, sigma, offset = p
            return [_biexp_irf_conv(t, A1, abs(tau1), A2, abs(tau2), t0, abs(sigma), offset)
                    for t in t_vec]
        end

        p0 = [A1_0, tau1_0, A2_0, tau2_0, 0.0, irf_width, offset0]
        prob = NonlinearCurveFitProblem(biexp_irf_model, p0, t_fit, signal_fit)
        sol = solve(prob)

        A1, tau1, A2, tau2, t0, sigma, offset = coef(sol)
        tau1, tau2, sigma = abs(tau1), abs(tau2), abs(sigma)

        fitted_vals = biexp_irf_model([A1, tau1, A2, tau2, t0, sigma, offset], t_fit)
    else
        biexp_model = n_exponentials(2)

        t_shifted = t_fit .- t_fit[1]

        p0 = [A1_0, tau1_0, A2_0, tau2_0, offset0]
        prob = NonlinearCurveFitProblem(biexp_model, p0, t_shifted, signal_fit)
        sol = solve(prob)

        A1, tau1, A2, tau2, offset = coef(sol)
        tau1, tau2 = abs(tau1), abs(tau2)
        t0 = t_fit[1]
        sigma = NaN

        fitted_vals = biexp_model([A1, tau1, A2, tau2, offset], t_shifted)
    end

    max_reasonable_tau = 10 * (t_fit[end] - t_fit[1])
    if tau1 > max_reasonable_tau || tau2 > max_reasonable_tau
        @warn "Biexponential fit may have failed - time constants unreasonably large. Consider using single exponential fit."
    end

    if tau1 > tau2
        tau1, tau2 = tau2, tau1
        A1, A2 = A2, A1
    end

    resid = signal_fit .- fitted_vals
    ss_res = sum(resid.^2)
    ss_tot = sum((signal_fit .- mean(signal_fit)).^2)
    rsquared = 1 - ss_res / ss_tot

    return BiexpDecayFit(
        tau1, tau2, A1, A2,
        t0, sigma, offset,
        signal_type, resid, rsquared
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

    individual_fits = [fit_exp_decay(tr; irf_width=irf_width) for tr in traces]

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

    y_pred_all = global_model(p_opt, x_all)
    idx = 1
    for i in 1:n_traces
        n_pts = length(traces[i].time)
        y_pred_i = y_pred_all[idx:idx+n_pts-1]
        y_data_i = traces[i].signal

        residuals_vec[i] = y_data_i .- y_pred_i

        ss_res = sum(residuals_vec[i].^2)
        ss_tot = sum((y_data_i .- mean(y_data_i)).^2)
        rsquared_individual[i] = 1 - ss_res / ss_tot

        idx += n_pts
    end

    ss_res_global = sum(sum.(r -> r.^2, residuals_vec))
    ss_tot_global = sum((y_all .- mean(y_all)).^2)
    rsquared_global = 1 - ss_res_global / ss_tot_global

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

function predict(fit::ExpDecayIRFFit, time::AbstractVector)
    if isnan(fit.sigma)
        return [t >= fit.t0 ? fit.amplitude * exp(-(t - fit.t0) / fit.tau) + fit.offset : fit.offset
                for t in time]
    else
        return [_exp_decay_irf_conv(t, fit.amplitude, fit.tau, fit.t0, fit.sigma) + fit.offset
                for t in time]
    end
end

predict(fit::ExpDecayIRFFit, trace::TATrace) = predict(fit, trace.time)

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

function predict(fit::BiexpDecayFit, time::AbstractVector)
    if isnan(fit.sigma)
        return [t >= fit.t0 ?
                fit.amplitude1 * exp(-(t - fit.t0) / fit.tau1) +
                fit.amplitude2 * exp(-(t - fit.t0) / fit.tau2) + fit.offset :
                fit.offset
                for t in time]
    else
        return [_biexp_irf_conv(t, fit.amplitude1, fit.tau1, fit.amplitude2, fit.tau2,
                                fit.t0, fit.sigma, fit.offset) for t in time]
    end
end

predict(fit::BiexpDecayFit, trace::TATrace) = predict(fit, trace.time)

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
# TA Spectrum Fitting (ESA + GSB peaks)
# =============================================================================

_ta_model_full(p, ν) = @. p[1] * exp(-4 * log(2) * ((ν - p[2]) / p[3])^2) -
                         p[4] * exp(-4 * log(2) * ((ν - p[5]) / p[6])^2) + p[7]

_ta_model_no_offset(p, ν) = @. p[1] * exp(-4 * log(2) * ((ν - p[2]) / p[3])^2) -
                               p[4] * exp(-4 * log(2) * ((ν - p[5]) / p[6])^2)

_ta_model_gsb_fwhm_fixed(p, ν, Γ_gsb) = @. p[1] * exp(-4 * log(2) * ((ν - p[2]) / p[3])^2) -
                                           p[4] * exp(-4 * log(2) * ((ν - p[5]) / Γ_gsb)^2) + p[6]

_ta_model_both_fixed(p, ν, Γ_gsb) = @. p[1] * exp(-4 * log(2) * ((ν - p[2]) / p[3])^2) -
                                        p[4] * exp(-4 * log(2) * ((ν - p[5]) / Γ_gsb)^2)

"""
    fit_ta_spectrum(spec::TASpectrum; region=nothing, gsb_fwhm=nothing, fit_offset=false) -> TASpectrumFit

Fit a pump-probe spectrum with ESA and GSB Gaussian peaks.
"""
function fit_ta_spectrum(spec::TASpectrum; region=nothing, gsb_fwhm=nothing, fit_offset::Bool=false)
    if isnothing(region)
        ν = spec.wavenumber
        y = spec.signal
    else
        mask = (spec.wavenumber .>= region[1]) .& (spec.wavenumber .<= region[2])
        ν = spec.wavenumber[mask]
        y = spec.signal[mask]
    end

    max_val = maximum(y)
    min_val = minimum(y)
    esa_idx = argmax(y)
    gsb_idx = argmin(y)

    ν_esa_init = ν[esa_idx]
    ν_gsb_init = ν[gsb_idx]
    A_esa_init = max_val
    A_gsb_init = abs(min_val)
    fwhm_init = 20.0

    fix_gsb_fwhm = !isnothing(gsb_fwhm)
    Γ_gsb_fixed = fix_gsb_fwhm ? gsb_fwhm : 0.0

    if fit_offset && !fix_gsb_fwhm
        p0 = [A_esa_init, ν_esa_init, fwhm_init, A_gsb_init, ν_gsb_init, fwhm_init, 0.0]
        fit = solve(NonlinearCurveFitProblem(_ta_model_full, p0, ν, y))
        p = coef(fit)
        esa_amp, esa_center, esa_fwhm = p[1], p[2], abs(p[3])
        gsb_amp, gsb_center, gsb_fwhm_fit = p[4], p[5], abs(p[6])
        offset = p[7]

    elseif fit_offset && fix_gsb_fwhm
        model_gsb_fixed(p, ν) = _ta_model_gsb_fwhm_fixed(p, ν, Γ_gsb_fixed)
        p0 = [A_esa_init, ν_esa_init, fwhm_init, A_gsb_init, ν_gsb_init, 0.0]
        fit = solve(NonlinearCurveFitProblem(model_gsb_fixed, p0, ν, y))
        p = coef(fit)
        esa_amp, esa_center, esa_fwhm = p[1], p[2], abs(p[3])
        gsb_amp, gsb_center = p[4], p[5]
        gsb_fwhm_fit = gsb_fwhm
        offset = p[6]

    elseif !fit_offset && !fix_gsb_fwhm
        p0 = [A_esa_init, ν_esa_init, fwhm_init, A_gsb_init, ν_gsb_init, fwhm_init]
        fit = solve(NonlinearCurveFitProblem(_ta_model_no_offset, p0, ν, y))
        p = coef(fit)
        esa_amp, esa_center, esa_fwhm = p[1], p[2], abs(p[3])
        gsb_amp, gsb_center, gsb_fwhm_fit = p[4], p[5], abs(p[6])
        offset = 0.0

    else  # !fit_offset && fix_gsb_fwhm
        model_both_fixed(p, ν) = _ta_model_both_fixed(p, ν, Γ_gsb_fixed)
        p0 = [A_esa_init, ν_esa_init, fwhm_init, A_gsb_init, ν_gsb_init]
        fit = solve(NonlinearCurveFitProblem(model_both_fixed, p0, ν, y))
        p = coef(fit)
        esa_amp, esa_center, esa_fwhm = p[1], p[2], abs(p[3])
        gsb_amp, gsb_center = p[4], p[5]
        gsb_fwhm_fit = gsb_fwhm
        offset = 0.0
    end

    y_fit = predict_ta_spectrum_internal(ν, esa_amp, esa_center, esa_fwhm,
                                          gsb_amp, gsb_center, gsb_fwhm_fit, offset)
    residuals = y .- y_fit
    ss_res = sum(residuals.^2)
    ss_tot = sum((y .- mean(y)).^2)
    rsquared = 1 - ss_res / ss_tot

    anharmonicity = gsb_center - esa_center

    return TASpectrumFit(esa_center, esa_fwhm, esa_amp,
                         gsb_center, gsb_fwhm_fit, gsb_amp,
                         offset, anharmonicity, rsquared, residuals)
end

function predict_ta_spectrum_internal(ν, A_esa, ν_esa, Γ_esa, A_gsb, ν_gsb, Γ_gsb, offset)
    return @. A_esa * exp(-4 * log(2) * ((ν - ν_esa) / Γ_esa)^2) -
              A_gsb * exp(-4 * log(2) * ((ν - ν_gsb) / Γ_gsb)^2) + offset
end

function predict(fit::TASpectrumFit, wavenumber::AbstractVector)
    return predict_ta_spectrum_internal(wavenumber,
                                        fit.esa_amplitude, fit.esa_center, fit.esa_fwhm,
                                        fit.gsb_amplitude, fit.gsb_center, fit.gsb_fwhm,
                                        fit.offset)
end

predict(fit::TASpectrumFit, spec::TASpectrum) = predict(fit, spec.wavenumber)
