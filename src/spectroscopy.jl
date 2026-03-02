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

# ============================================================================
# SPECTRAL MATH FUNCTIONS
# ============================================================================

"""
    savitzky_golay_smooth(y; window=11, order=3)

Apply Savitzky-Golay smoothing filter to data.

Uses polynomial fitting within a sliding window to smooth data while preserving
peak shape and width better than moving-average smoothing (Savitzky & Golay, 1964).

# Arguments
- `y::AbstractVector{<:Real}`: Input data vector.

# Keywords
- `window::Int=11`: Window size (must be odd and ≥ order + 2).
- `order::Int=3`: Polynomial order for the local fit.

# Returns
- `Vector{Float64}`: Smoothed data vector (same length as input).

# Examples
```julia
y_noisy = sin.(0:0.1:2π) .+ 0.1 * randn(63)
y_smooth = savitzky_golay_smooth(y_noisy; window=11, order=3)
```
"""
function savitzky_golay_smooth(y::AbstractVector{<:Real}; window::Int=11, order::Int=3)
    return _sg_filter(y, window, order).y
end

"""
    derivative(y; order=1, window=11, poly_order=3)
    derivative(x, y; order=1, window=11, poly_order=3)

Compute the derivative of a signal using Savitzky-Golay differentiation.

When `x` is provided, the derivative is correctly scaled by the x-spacing
(using the median point spacing as the rate parameter).

# Arguments
- `x::AbstractVector{<:Real}`: (optional) Independent variable (e.g., wavenumber, wavelength).
- `y::AbstractVector{<:Real}`: Signal to differentiate.

# Keywords
- `order::Int=1`: Derivative order (1 = first derivative, 2 = second, etc.).
- `window::Int=11`: Savitzky-Golay window size (must be odd, ≥ poly_order + 2).
- `poly_order::Int=3`: Polynomial order for the SG filter (must be ≥ `order`).

# Returns
- `Vector{Float64}`: Derivative of the input signal.

# Examples
```julia
x = 400.0:0.5:800.0
y = @. 100 * exp(-(x - 520)^2 / (2 * 20^2))
dy = derivative(x, y; order=1)          # First derivative (correctly scaled)
d2y = derivative(x, y; order=2)         # Second derivative
```
"""
function derivative(y::AbstractVector{<:Real}; order::Int=1, window::Int=11, poly_order::Int=3)
    return _sg_filter(y, window, poly_order, deriv=order).y
end

function derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                    order::Int=1, window::Int=11, poly_order::Int=3)
    dx = median(diff(collect(x)))
    rate = 1.0 / abs(dx)
    return _sg_filter(y, window, poly_order, deriv=order, rate=rate).y
end

"""
    band_area(x, y, x_min, x_max)

Compute the integrated area under a spectrum within a given range using
trapezoidal integration.

# Arguments
- `x::AbstractVector{<:Real}`: x-axis values (e.g., wavenumber, wavelength).
- `y::AbstractVector{<:Real}`: y-axis values (e.g., intensity).
- `x_min::Real`: Lower bound of integration range.
- `x_max::Real`: Upper bound of integration range.

# Returns
- `Float64`: Integrated area.

# Examples
```julia
x = 400.0:0.5:800.0
y = @. 100 * exp(-(x - 520)^2 / (2 * 10^2))
area = band_area(x, y, 480.0, 560.0)
```
"""
function band_area(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                   x_min::Real, x_max::Real)
    x_lo, x_hi = minmax(x_min, x_max)
    mask = findall(xi -> x_lo <= xi <= x_hi, x)
    length(mask) >= 2 || throw(ArgumentError(
        "Fewer than 2 points in range [$x_lo, $x_hi]"))
    xr = x[mask]
    yr = y[mask]
    area = 0.0
    for i in 2:length(xr)
        area += 0.5 * (yr[i] + yr[i-1]) * abs(xr[i] - xr[i-1])
    end
    return area
end

"""
    normalize_area(x, y)

Normalize a spectrum so its total integrated area equals 1.

# Arguments
- `x::AbstractVector{<:Real}`: x-axis values.
- `y::AbstractVector{<:Real}`: y-axis values.

# Returns
- `Vector{Float64}`: Area-normalized y-values.

# Examples
```julia
x = 1.0:0.1:10.0
y = ones(length(x))
y_norm = normalize_area(x, y)  # Total area ≈ 1.0
```
"""
function normalize_area(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    total = band_area(x, y, x[1], x[end])
    abs(total) < eps(Float64) && throw(ArgumentError(
        "Total area is zero; cannot normalize"))
    return Float64.(y ./ total)
end

"""
    normalize_to_peak(x, y, position; tolerance=5.0)

Normalize a spectrum by dividing by the intensity at a specified peak position.

Finds the data point nearest to `position` within `tolerance` and divides
the entire spectrum by that point's intensity.

# Arguments
- `x::AbstractVector{<:Real}`: x-axis values.
- `y::AbstractVector{<:Real}`: y-axis values.
- `position::Real`: Target x-position for normalization (e.g., 520 cm-1 for Si).

# Keywords
- `tolerance::Real=5.0`: Maximum allowed distance from `position` to nearest data point.

# Returns
- `Vector{Float64}`: Peak-normalized y-values.

# Examples
```julia
x = 400.0:1.0:800.0
y = @. 50 * exp(-(x - 520)^2 / (2 * 10^2))
y_norm = normalize_to_peak(x, y, 520.0)  # y at x≈520 becomes 1.0
```
"""
function normalize_to_peak(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                           position::Real; tolerance::Real=5.0)
    idx = _find_nearest_idx(collect(x), position)
    dist = abs(x[idx] - position)
    dist <= tolerance || throw(ArgumentError(
        "No data point within tolerance=$tolerance of position=$position " *
        "(nearest at x=$(x[idx]), distance=$dist)"))
    val = y[idx]
    abs(val) < eps(Float64) && throw(ArgumentError(
        "Intensity at position=$position is zero; cannot normalize"))
    return Float64.(y ./ val)
end

"""
    estimate_snr(y)

Estimate the signal-to-noise ratio using the DER-SNR method.

Uses the second-order finite difference of adjacent pixels to estimate noise
(Stoehr et al., 2008, "DER_SNR: A Simple & General Spectroscopic Signal-to-Noise
Measurement Algorithm").

# Arguments
- `y::AbstractVector{<:Real}`: Spectral intensity values.

# Returns
- `NamedTuple{(:snr, :signal, :noise)}`: Estimated SNR, signal level, and noise level.

# Examples
```julia
y_clean = 100 * ones(100)
y_noisy = y_clean .+ 2 * randn(100)
result = estimate_snr(y_noisy)
result.snr  # ≈ 50
```
"""
function estimate_snr(y::AbstractVector{<:Real})
    n = length(y)
    n >= 4 || throw(ArgumentError("Need at least 4 points to estimate SNR"))
    noise_arr = similar(y, n - 2)
    for i in 2:(n - 1)
        noise_arr[i - 1] = abs(2 * y[i] - y[i - 1] - y[i + 1])
    end
    noise = median(noise_arr) / 0.6744897501960817
    noise = noise / sqrt(6.0)
    signal = median(y)
    snr = noise > 0 ? signal / noise : Inf
    return (snr=snr, signal=signal, noise=noise)
end
