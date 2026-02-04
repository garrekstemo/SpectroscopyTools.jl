# Core spectroscopy analysis functions

using SavitzkyGolay: savitzky_golay as _sg_filter

# ============================================================================
# BASIC UTILITY FUNCTIONS
# ============================================================================

"""
    normalize(x)

Normalize array by maximum absolute value. Returns zeros if max is zero.
"""
normalize(x) = (m = maximum(abs.(x)); m == 0 ? zeros(length(x)) : x ./ m)

"""
    time_index(times, t_target)

Find index closest to target time value.
"""
time_index(times, t_target) = argmin(abs.(times .- t_target))

# ============================================================================
# TRANSMITTANCE <-> ABSORBANCE CONVERSIONS
# ============================================================================

"""
    transmittance_to_absorbance(T; percent=false)

Convert transmittance to absorbance: `A = -log10(T)`.
Input `T` is fractional (0 to 1). Use `percent=true` for percent transmittance.
"""
function transmittance_to_absorbance(T::Real; percent::Bool=false)
    T_frac = percent ? T / 100 : T
    T_frac > 0 || throw(ArgumentError("Transmittance must be positive (got $T_frac)"))
    return -log10(T_frac)
end

function transmittance_to_absorbance(T::AbstractVector; percent::Bool=false)
    return transmittance_to_absorbance.(T; percent=percent)
end

"""
    absorbance_to_transmittance(A; percent=false)

Convert absorbance to transmittance: `T = 10^(-A)`.
"""
function absorbance_to_transmittance(A::Real; percent::Bool=false)
    T = 10.0^(-A)
    return percent ? T * 100 : T
end

function absorbance_to_transmittance(A::AbstractVector; percent::Bool=false)
    return absorbance_to_transmittance.(A; percent=percent)
end

# ============================================================================
# TRANSIENT ABSORPTION ANALYSIS
# ============================================================================

"""
    calc_ΔA(trace; mode=:transmission)

Calculate change in absorbance (ΔA) from transient absorption data.

Modes: `:difference`, `:transmission` (-ΔT/T), `:OD` (ΔOD = -log10(T_on/T_off))
"""
function calc_ΔA(trace; mode=:transmission)
    if mode == :OD
        return -log10.(trace.on ./ trace.off)
    elseif mode == :transmission
        return -trace.diff ./ trace.off
    else
        return trace.diff
    end
end

"""
    fit_decay_trace(time, signal; truncate_time=0.0, initial_tau, t0=nothing, fit_func=single_exponential)

Fit a decay trace after truncating early points.
Returns `(fit=sol, start_idx=idx)`.
"""
function fit_decay_trace(time, signal; truncate_time=0.0, initial_tau, t0=nothing, fit_func=single_exponential)
    t0 = something(t0, time[argmax(abs.(signal))])
    start_idx = searchsortedfirst(time, t0 + truncate_time)

    t_fit = @view time[start_idx:end]
    y_fit = @view signal[start_idx:end]
    mask = isfinite.(t_fit) .& isfinite.(y_fit)
    t_fit, y_fit = t_fit[mask], y_fit[mask]

    fit = solve(NonlinearCurveFitProblem(fit_func, [y_fit[1], initial_tau, minimum(y_fit)], t_fit, y_fit))

    return (fit=fit, start_idx=start_idx)
end

"""
    extract_tau(fit; idx=2, digits=2)

Extract time constant and standard error from exponential fit.
"""
extract_tau(fit; idx=2, digits=2) =
    (τ = round(coef(fit)[idx], digits=digits),
     στ = round(stderror(fit)[idx], digits=digits))

"""
    fit_global_decay(time, signals; start_idx, initial_tau)

Fit a global decay model with shared τ across multiple traces.
"""
function fit_global_decay(time, signals; start_idx, initial_tau)
    t_fit = @view time[start_idx:end]
    signal_fits = [@view s[start_idx:end] for s in signals]

    mask = isfinite.(t_fit)
    for s in signal_fits
        mask .&= isfinite.(s)
    end
    t_fit = t_fit[mask]
    signal_fits = [s[mask] for s in signal_fits]

    function global_model(p, t)
        τ = abs(p[1])
        n = length(signal_fits)
        models = [
            single_exponential((p[1 + i], τ, p[1 + n + i]), t)
            for i in 1:n
        ]
        return vcat(models...)
    end

    y_fit = vcat(signal_fits...)
    p0 = vcat(
        initial_tau,
        [signal_fits[i][1] for i in eachindex(signal_fits)]...,
        [minimum(signal_fits[i]) for i in eachindex(signal_fits)]...,
    )

    fit = solve(NonlinearCurveFitProblem(global_model, p0, t_fit, y_fit))
    return (fit=fit, start_idx=start_idx)
end

# ============================================================================
# SPECTRUM SUBTRACTION
# ============================================================================

"""
    subtract_spectrum(sample, reference; scale=1.0, interpolate=false)

Subtract a reference spectrum from a sample spectrum.
Expects objects with `.x` and `.y` fields.
"""
function subtract_spectrum(sample, reference; scale::Real=1.0, interpolate=false)
    if !interpolate
        if length(sample.x) != length(reference.x)
            error("""
                  Grid mismatch: sample has $(length(sample.x)) points, reference has $(length(reference.x)).

                  Use `interpolate=true` to interpolate the reference onto the sample grid:
                      subtract_spectrum(sample, reference, interpolate=true)
                  """)
        end

        max_diff = maximum(abs.(sample.x .- reference.x))
        if max_diff > 0.01
            error("""
                  Grid mismatch: x-values differ by up to $(round(max_diff, digits=3)) cm⁻¹.

                  Use `interpolate=true` to interpolate the reference onto the sample grid:
                      subtract_spectrum(sample, reference, interpolate=true)
                  """)
        end
    end

    if interpolate
        ref_interp = similar(sample.y)
        for i in eachindex(sample.x)
            x_target = sample.x[i]
            idx = searchsortedfirst(reference.x, x_target)
            if idx == 1
                ref_interp[i] = reference.y[1]
            elseif idx > length(reference.x)
                ref_interp[i] = reference.y[end]
            else
                x1, x2 = reference.x[idx-1], reference.x[idx]
                y1, y2 = reference.y[idx-1], reference.y[idx]
                t = (x_target - x1) / (x2 - x1)
                ref_interp[i] = y1 + t * (y2 - y1)
            end
        end
        y_subtracted = sample.y .- scale .* ref_interp
    else
        y_subtracted = sample.y .- scale .* reference.y
    end

    return (x=sample.x, y=y_subtracted)
end

# ============================================================================
# LINEAR BASELINE CORRECTION
# ============================================================================

"""
    linear_baseline_correction(x, y, anchor_points; window=5.0)

Apply linear baseline correction using anchor points.
"""
function linear_baseline_correction(x, y, anchor_points; window=5.0)
    lo, hi = anchor_points

    lo_mask = (lo - window) .<= x .<= (lo + window)
    hi_mask = (hi - window) .<= x .<= (hi + window)

    x_lo, y_lo = Statistics.mean(x[lo_mask]), Statistics.mean(y[lo_mask])
    x_hi, y_hi = Statistics.mean(x[hi_mask]), Statistics.mean(y[hi_mask])

    m = (y_hi - y_lo) / (x_hi - x_lo)
    b = y_lo - m * x_lo

    baseline = @. m * x + b

    return (x=x, y=y .- baseline, baseline=baseline)
end

# ============================================================================
# SMOOTHING AND PEAK ANALYSIS
# ============================================================================

"""
    smooth_data(y; window=3)

Apply moving average smoothing to data.
"""
function smooth_data(y; window=3)
    n = length(y)
    smoothed = similar(y)
    half_w = window ÷ 2

    for i in eachindex(y)
        left_extend = min(half_w, i - 1)
        right_extend = min(half_w, n - i)
        start_idx = i - left_extend
        end_idx = i + right_extend
        smoothed[i] = Statistics.mean(@view y[start_idx:end_idx])
    end
    return smoothed
end

"""
    savitzky_golay(y; window=5, order=2)

Apply Savitzky-Golay filter for smoothing while preserving peak shape.
"""
function savitzky_golay(y; window=5, order=2)
    _sg_filter(y, window, order).y
end

"""
    calc_fwhm(x, y; smooth_window=5)

Calculate full width at half maximum (FWHM) of the dominant positive peak.
"""
function calc_fwhm(x, y; smooth_window=5)
    y_smooth = smooth_window > 1 ? savitzky_golay(y; window=smooth_window) : y

    peak_idx = argmax(y_smooth)
    peak_val = y_smooth[peak_idx]
    half_max = peak_val / 2

    left_x = x[1]
    for i in (peak_idx-1):-1:1
        if y_smooth[i] <= half_max
            α = (half_max - y_smooth[i]) / (y_smooth[i+1] - y_smooth[i])
            left_x = x[i] + α * (x[i+1] - x[i])
            break
        end
    end

    right_x = x[end]
    for i in (peak_idx+1):length(y_smooth)
        if y_smooth[i] <= half_max
            α = (half_max - y_smooth[i-1]) / (y_smooth[i] - y_smooth[i-1])
            right_x = x[i-1] + α * (x[i] - x[i-1])
            break
        end
    end

    fwhm = abs(right_x - left_x)

    return (
        peak_position = x[peak_idx],
        peak_value = peak_val,
        fwhm = fwhm,
        bounds = (left_x, right_x)
    )
end
