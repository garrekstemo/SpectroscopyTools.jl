# Core spectroscopy analysis functions

# ============================================================================
# BASIC UTILITY FUNCTIONS
# ============================================================================

"""
    normalize(x)

Normalize array by maximum absolute value. Returns zeros if max is zero.
"""
normalize(x) = (m = maximum(abs.(x)); m == 0 ? zero(x) : x ./ m)

"""
    time_index(times, t_target)

Find index closest to target time value.
"""
time_index(times, t_target) = _find_nearest_idx(times, t_target)

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
# SPECTRUM SUBTRACTION
# ============================================================================

"""
    subtract_spectrum(sample, reference; scale=1.0, interpolate=false)

Subtract a reference spectrum from a sample spectrum.

Accepts `AbstractSpectroscopyData` types (uses `xdata`/`ydata` interface)
or any objects with `.x` and `.y` fields.

Returns `(x=..., y=...)` NamedTuple.
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

# Typed dispatch: AbstractSpectroscopyData → xdata/ydata interface
function subtract_spectrum(sample::AbstractSpectroscopyData,
                           reference::AbstractSpectroscopyData; kwargs...)
    subtract_spectrum((x=xdata(sample), y=ydata(sample)),
                      (x=xdata(reference), y=ydata(reference)); kwargs...)
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
        smoothed[i] = mean(@view y[start_idx:end_idx])
    end
    return smoothed
end

"""
    calc_fwhm(x, y; smooth_window=5)

Calculate full width at half maximum (FWHM) of the dominant positive peak.
"""
function calc_fwhm(x, y; smooth_window=5)
    y_smooth = smooth_window > 1 ? _sg_filter(y, smooth_window, 2).y : y

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
