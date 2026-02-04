"""
Data types for spectroscopy analysis.
"""

# =============================================================================
# Abstract types for spectroscopy data
# =============================================================================

"""
    AbstractSpectroscopyData

Abstract base type for all spectroscopy data types.

Subtypes must implement the following interface:
- `xdata(d)` — Primary x-axis data (Vector{Float64})
- `ydata(d)` — Signal data (Vector{Float64}) or secondary axis for 2D
- `zdata(d)` — Matrix data for 2D types, `nothing` for 1D
- `xlabel(d)` — X-axis label string
- `ylabel(d)` — Y-axis or signal label string
- `is_matrix(d)` — Whether data is 2D (returns Bool)

Optional (have defaults):
- `source_file(d)` — Source filename (default: `""`)
- `npoints(d)` — Number of data points (default: `length(xdata(d))`, tuple for 2D)
- `title(d)` — Display title (default: `source_file(d)`)

This enables uniform handling in data viewers and plotting functions
while maintaining semantic field names in each concrete type.

# Example
```julia
# Works for any spectroscopy data type
function plot_data(data::AbstractSpectroscopyData)
    if is_matrix(data)
        heatmap(xdata(data), ydata(data), zdata(data)')
    else
        lines(xdata(data), ydata(data))
    end
end
```
"""
abstract type AbstractSpectroscopyData end

# Interface functions (default implementations raise errors)
xdata(::AbstractSpectroscopyData) = error("xdata not implemented for this type")
ydata(::AbstractSpectroscopyData) = error("ydata not implemented for this type")
zdata(::AbstractSpectroscopyData) = nothing  # Default: not a matrix
xlabel(::AbstractSpectroscopyData) = "X"
ylabel(::AbstractSpectroscopyData) = "Y"
zlabel(::AbstractSpectroscopyData) = "Signal"
is_matrix(::AbstractSpectroscopyData) = false  # Default: 1D data

"""
    source_file(d::AbstractSpectroscopyData) -> String

Return the source filename for the data.
"""
source_file(::AbstractSpectroscopyData) = ""

"""
    npoints(d::AbstractSpectroscopyData) -> Int
    npoints(d::TAMatrix) -> Tuple{Int,Int}

Return the number of data points. For 1D data, returns an `Int`.
For 2D data (`TAMatrix`), returns `(n_time, n_wavelength)`.
"""
npoints(d::AbstractSpectroscopyData) = length(xdata(d))

"""
    title(d::AbstractSpectroscopyData) -> String

Return a display title for the data. Defaults to `source_file(d)`.
"""
title(d::AbstractSpectroscopyData) = source_file(d)

# =============================================================================
# Pump-probe types
# =============================================================================

"""
    AxisType

Enum indicating what the x-axis of pump-probe data represents.

- `time_axis` — Kinetics measurement, x = time in picoseconds
- `wavelength_axis` — Spectral measurement, x = wavelength in nanometers
"""
@enum AxisType begin
    time_axis        # Kinetics: x = time (ps)
    wavelength_axis  # Spectrum: x = wavelength (nm)
end

"""
    PumpProbeData

Raw pump-probe measurement data.

# Fields
- `time::Vector{Float64}`: X-axis data (time in ps OR wavelength in nm, see `axis_type`)
- `on::Matrix{Float64}`: Pump-on signal (rows = x points, cols = channels 0-3)
- `off::Matrix{Float64}`: Pump-off signal
- `diff::Matrix{Float64}`: Lock-in difference (ON - OFF, computed by instrument)
- `timestamp::String`: Acquisition timestamp
- `axis_type::AxisType`: What the x-axis represents (time or wavelength)
"""
struct PumpProbeData
    time::Vector{Float64}
    on::Matrix{Float64}
    off::Matrix{Float64}
    diff::Matrix{Float64}
    timestamp::String
    axis_type::AxisType
end

"""
    xaxis(data::PumpProbeData) -> Vector{Float64}

Return the x-axis data (time or wavelength depending on measurement type).
"""
xaxis(d::PumpProbeData) = d.time

"""
    xaxis_label(data::PumpProbeData) -> String

Return the appropriate x-axis label based on measurement type.
"""
xaxis_label(d::PumpProbeData) = d.axis_type == time_axis ? "Time (ps)" : "Wavelength (nm)"

# =============================================================================
# Transient Absorption types (unified API)
# =============================================================================

"""
    TATrace <: AbstractSpectroscopyData

Single-wavelength transient absorption kinetic trace.

# Fields
- `time::Vector{Float64}`: Time axis (ps)
- `signal::Vector{Float64}`: ΔA signal
- `wavelength::Float64`: Probe wavelength, NaN if unknown
- `metadata::Dict{Symbol,Any}`: Additional info
"""
struct TATrace <: AbstractSpectroscopyData
    time::Vector{Float64}
    signal::Vector{Float64}
    wavelength::Float64
    metadata::Dict{Symbol,Any}
end

# Constructor with default wavelength
TATrace(time, signal; wavelength=NaN, metadata=Dict{Symbol,Any}()) =
    TATrace(time, signal, wavelength, metadata)

# AbstractSpectroscopyData interface
xdata(t::TATrace) = t.time
ydata(t::TATrace) = t.signal
xlabel(::TATrace) = "Time (ps)"
ylabel(::TATrace) = "ΔA"
source_file(t::TATrace) = get(t.metadata, :filename, "")

# Pretty printing
function Base.show(io::IO, t::TATrace)
    n = length(t.time)
    t_range = "$(round(minimum(t.time), digits=2)) to $(round(maximum(t.time), digits=2)) ps"
    λ_str = isnan(t.wavelength) ? "unknown" : "$(t.wavelength)"
    filename = get(t.metadata, :filename, "")

    print(io, "TATrace: $n points, $t_range")
    !isempty(filename) && print(io, " ($(filename))")
end

function Base.show(io::IO, ::MIME"text/plain", t::TATrace)
    println(io, "TATrace")
    println(io, "  Time points: $(length(t.time))")
    println(io, "  Time range:  $(round(minimum(t.time), digits=2)) to $(round(maximum(t.time), digits=2)) ps")
    println(io, "  Wavelength:  $(isnan(t.wavelength) ? "unknown" : t.wavelength)")
    if haskey(t.metadata, :filename)
        println(io, "  File:        $(t.metadata[:filename])")
    end
    if haskey(t.metadata, :mode)
        println(io, "  Mode:        $(t.metadata[:mode])")
    end
end

"""
    TASpectrum <: AbstractSpectroscopyData

Transient absorption spectrum at a fixed time delay.

# Fields
- `wavenumber::Vector{Float64}`: Wavenumber axis (cm⁻¹)
- `signal::Vector{Float64}`: ΔA signal
- `time_delay::Float64`: Time delay (ps), NaN if unknown
- `metadata::Dict{Symbol,Any}`: Additional info
"""
struct TASpectrum <: AbstractSpectroscopyData
    wavenumber::Vector{Float64}
    signal::Vector{Float64}
    time_delay::Float64
    metadata::Dict{Symbol,Any}
end

# Constructor with default time_delay
TASpectrum(wavenumber, signal; time_delay=NaN, metadata=Dict{Symbol,Any}()) =
    TASpectrum(wavenumber, signal, time_delay, metadata)

# AbstractSpectroscopyData interface
xdata(s::TASpectrum) = s.wavenumber
ydata(s::TASpectrum) = s.signal
xlabel(::TASpectrum) = "Wavenumber (cm⁻¹)"
ylabel(::TASpectrum) = "ΔA"
source_file(s::TASpectrum) = get(s.metadata, :filename, "")

# Pretty printing
function Base.show(io::IO, s::TASpectrum)
    n = length(s.wavenumber)
    ν_range = "$(round(Int, minimum(s.wavenumber))) - $(round(Int, maximum(s.wavenumber))) cm⁻¹"
    filename = get(s.metadata, :filename, "")

    print(io, "TASpectrum: $n points, $ν_range")
    !isempty(filename) && print(io, " ($(filename))")
end

function Base.show(io::IO, ::MIME"text/plain", s::TASpectrum)
    println(io, "TASpectrum")
    println(io, "  Points:      $(length(s.wavenumber))")
    println(io, "  Range:       $(round(Int, minimum(s.wavenumber))) - $(round(Int, maximum(s.wavenumber))) cm⁻¹")
    println(io, "  Time delay:  $(isnan(s.time_delay) ? "unknown" : "$(s.time_delay) ps")")
    if haskey(s.metadata, :filename)
        println(io, "  File:        $(s.metadata[:filename])")
    end
    if haskey(s.metadata, :mode)
        println(io, "  Mode:        $(s.metadata[:mode])")
    end
end

"""
    TASpectrumFit

Result of pump-probe spectrum fitting (ESA + GSB peaks).

# Fields
- `esa_center`, `esa_fwhm`, `esa_amplitude`: ESA peak parameters
- `gsb_center`, `gsb_fwhm`, `gsb_amplitude`: GSB peak parameters
- `offset`, `anharmonicity`, `rsquared`, `residuals`: Fit metadata
"""
struct TASpectrumFit
    esa_center::Float64
    esa_fwhm::Float64
    esa_amplitude::Float64
    gsb_center::Float64
    gsb_fwhm::Float64
    gsb_amplitude::Float64
    offset::Float64
    anharmonicity::Float64
    rsquared::Float64
    residuals::Vector{Float64}
end

function Base.show(io::IO, fit::TASpectrumFit)
    print(io, "TASpectrumFit: ESA $(round(Int, fit.esa_center)) cm⁻¹, GSB $(round(Int, fit.gsb_center)) cm⁻¹, Δν = $(round(fit.anharmonicity, digits=1)) cm⁻¹")
end

function Base.show(io::IO, ::MIME"text/plain", fit::TASpectrumFit)
    println(io, "TASpectrumFit")
    println(io)
    println(io, "  Peak        Center (cm⁻¹)    FWHM (cm⁻¹)    Amplitude")
    println(io, "  ─────────────────────────────────────────────────────────")
    println(io, "  ESA         $(lpad(round(fit.esa_center, digits=1), 10))    $(lpad(round(fit.esa_fwhm, digits=1), 10))    $(lpad(round(fit.esa_amplitude, sigdigits=3), 10))")
    println(io, "  GSB         $(lpad(round(fit.gsb_center, digits=1), 10))    $(lpad(round(fit.gsb_fwhm, digits=1), 10))    $(lpad(round(fit.gsb_amplitude, sigdigits=3), 10))")
    println(io)
    println(io, "  Anharmonicity: $(round(fit.anharmonicity, digits=1)) cm⁻¹")
    println(io, "  Offset:        $(round(fit.offset, sigdigits=3))")
    println(io)
    print(io, "  R² = $(round(fit.rsquared, digits=5))")
end

# =============================================================================
# Exponential decay fit types
# =============================================================================

"""
    ExpDecayFit

Result of exponential decay fitting (no IRF).

# Fields
- `amplitude`, `tau`, `offset`, `t0`, `signal_type`, `residuals`, `rsquared`
"""
struct ExpDecayFit
    amplitude::Float64
    tau::Float64
    offset::Float64
    t0::Int
    signal_type::Symbol
    residuals::Vector{Float64}
    rsquared::Float64
end

function Base.show(io::IO, fit::ExpDecayFit)
    signal_str = fit.signal_type == :esa ? "ESA" : "GSB"
    print(io, "ExpDecayFit: τ = $(round(fit.tau, digits=2)), R² = $(round(fit.rsquared, digits=4)) ($signal_str)")
end

function Base.show(io::IO, ::MIME"text/plain", fit::ExpDecayFit)
    signal_str = fit.signal_type == :esa ? "ESA (positive)" : "GSB (negative)"

    println(io, "ExpDecayFit ($signal_str)")
    println(io)
    println(io, "  Parameter       Value")
    println(io, "  ─────────────────────────────")
    println(io, "  τ          $(lpad(round(fit.tau, digits=3), 10))")
    println(io, "  Amplitude  $(lpad(round(fit.amplitude, sigdigits=4), 10))")
    println(io, "  Offset     $(lpad(round(fit.offset, sigdigits=4), 10))")
    println(io, "  t₀ index   $(lpad(fit.t0, 10))")
    println(io)
    print(io,   "  R² = $(round(fit.rsquared, digits=5))")
end

"""
    ExpDecayIRFFit

Result of exponential decay fitting with instrument response function convolution.

# Fields
- `amplitude`, `tau`, `t0`, `sigma`, `offset`, `signal_type`, `residuals`, `rsquared`
"""
struct ExpDecayIRFFit
    amplitude::Float64
    tau::Float64
    t0::Float64
    sigma::Float64
    offset::Float64
    signal_type::Symbol
    residuals::Vector{Float64}
    rsquared::Float64
end

function Base.show(io::IO, fit::ExpDecayIRFFit)
    signal_str = fit.signal_type == :esa ? "ESA" : "GSB"
    print(io, "ExpDecayIRFFit: τ = $(round(fit.tau, digits=2)) ps, R² = $(round(fit.rsquared, digits=4)) ($signal_str)")
end

function Base.show(io::IO, ::MIME"text/plain", fit::ExpDecayIRFFit)
    signal_str = fit.signal_type == :esa ? "ESA (positive)" : "GSB (negative)"
    has_irf = !isnan(fit.sigma)

    println(io, "ExpDecayIRFFit ($signal_str)")
    println(io)
    println(io, "  Parameter       Value")
    println(io, "  ─────────────────────────────")
    println(io, "  τ          $(lpad(round(fit.tau, digits=3), 10)) ps")
    if has_irf
        println(io, "  σ_IRF      $(lpad(round(fit.sigma, digits=3), 10)) ps")
    end
    println(io, "  t₀         $(lpad(round(fit.t0, digits=3), 10)) ps")
    println(io, "  Amplitude  $(lpad(round(fit.amplitude, sigdigits=4), 10))")
    println(io, "  Offset     $(lpad(round(fit.offset, sigdigits=4), 10))")
    println(io)
    print(io, "  R² = $(round(fit.rsquared, digits=5))")
end

"""
    BiexpDecayFit

Result of biexponential decay fitting.

Time constants are ordered so τ₁ < τ₂ (fast, slow).
"""
struct BiexpDecayFit
    tau1::Float64
    tau2::Float64
    amplitude1::Float64
    amplitude2::Float64
    t0::Float64
    sigma::Float64
    offset::Float64
    signal_type::Symbol
    residuals::Vector{Float64}
    rsquared::Float64
end

function Base.show(io::IO, fit::BiexpDecayFit)
    signal_str = fit.signal_type == :esa ? "ESA" : "GSB"
    print(io, "BiexpDecayFit: τ₁ = $(round(fit.tau1, digits=2)) ps, τ₂ = $(round(fit.tau2, digits=2)) ps, R² = $(round(fit.rsquared, digits=4)) ($signal_str)")
end

function Base.show(io::IO, ::MIME"text/plain", fit::BiexpDecayFit)
    signal_str = fit.signal_type == :esa ? "ESA (positive)" : "GSB (negative)"
    has_irf = !isnan(fit.sigma)

    total_amp = abs(fit.amplitude1) + abs(fit.amplitude2)
    frac1 = total_amp > 0 ? round(100 * abs(fit.amplitude1) / total_amp, digits=1) : 0.0
    frac2 = total_amp > 0 ? round(100 * abs(fit.amplitude2) / total_amp, digits=1) : 0.0

    println(io, "BiexpDecayFit ($signal_str)")
    println(io)
    println(io, "  Parameter       Value        Weight")
    println(io, "  ──────────────────────────────────────")
    println(io, "  τ₁         $(lpad(round(fit.tau1, digits=3), 10)) ps   $(lpad(frac1, 5))%")
    println(io, "  τ₂         $(lpad(round(fit.tau2, digits=3), 10)) ps   $(lpad(frac2, 5))%")
    println(io, "  A₁         $(lpad(round(fit.amplitude1, sigdigits=4), 10))")
    println(io, "  A₂         $(lpad(round(fit.amplitude2, sigdigits=4), 10))")
    if has_irf
        println(io, "  σ_IRF      $(lpad(round(fit.sigma, digits=3), 10)) ps")
    end
    println(io, "  t₀         $(lpad(round(fit.t0, digits=3), 10)) ps")
    println(io, "  Offset     $(lpad(round(fit.offset, sigdigits=4), 10))")
    println(io)
    print(io, "  R² = $(round(fit.rsquared, digits=5))")
end

"""
    MultiexpDecayFit

Result of multi-exponential decay fitting (n >= 1 components).

# Fields
- `taus::Vector{Float64}`: Time constants (sorted fast->slow)
- `amplitudes::Vector{Float64}`: Corresponding amplitudes
- `t0`, `sigma`, `offset`, `signal_type`, `residuals`, `rsquared`

# Derived properties
- `n_exp(fit)`: Number of exponential components
- `weights(fit)`: Relative amplitude weights (normalized to 100%)
"""
struct MultiexpDecayFit
    taus::Vector{Float64}
    amplitudes::Vector{Float64}
    t0::Float64
    sigma::Float64
    offset::Float64
    signal_type::Symbol
    residuals::Vector{Float64}
    rsquared::Float64
end

n_exp(fit::MultiexpDecayFit) = length(fit.taus)
function weights(fit::MultiexpDecayFit)
    total = sum(abs, fit.amplitudes)
    return total > 0 ? abs.(fit.amplitudes) ./ total : fill(1.0 / length(fit.amplitudes), length(fit.amplitudes))
end

function Base.show(io::IO, fit::MultiexpDecayFit)
    n = length(fit.taus)
    signal_str = fit.signal_type == :esa ? "ESA" : "GSB"
    tau_str = join(["$(round(τ, digits=2))" for τ in fit.taus], ", ")
    print(io, "MultiexpDecayFit: n=$n, τ = [$tau_str] ps, R² = $(round(fit.rsquared, digits=4)) ($signal_str)")
end

function Base.show(io::IO, ::MIME"text/plain", fit::MultiexpDecayFit)
    n = length(fit.taus)
    signal_str = fit.signal_type == :esa ? "ESA (positive)" : "GSB (negative)"
    has_irf = !isnan(fit.sigma)
    w = weights(fit)

    println(io, "MultiexpDecayFit ($signal_str, $n components)")
    println(io)
    println(io, "  Component     τ (ps)        Amplitude      Weight")
    println(io, "  ────────────────────────────────────────────────────")
    for i in 1:n
        println(io, "  $(lpad(i, 5))       $(lpad(round(fit.taus[i], digits=3), 8))    $(lpad(round(fit.amplitudes[i], sigdigits=4), 12))    $(lpad(round(100*w[i], digits=1), 5))%")
    end
    println(io)
    if has_irf
        println(io, "  σ_IRF      $(lpad(round(fit.sigma, digits=3), 10)) ps")
    end
    println(io, "  t₀         $(lpad(round(fit.t0, digits=3), 10)) ps")
    println(io, "  Offset     $(lpad(round(fit.offset, sigdigits=4), 10))")
    println(io)
    print(io, "  R² = $(round(fit.rsquared, digits=5))")
end

"""
    PumpProbeResult

Complete result from pump-probe analysis.
"""
struct PumpProbeResult
    t::Vector{Float64}
    signal::Vector{Float64}
    fit::ExpDecayIRFFit
    fit_curve::Vector{Float64}
    filename::String
end

function Base.show(io::IO, r::PumpProbeResult)
    println(io, "PumpProbeResult: $(r.filename)")
    println(io, "  Signal: $(r.fit.signal_type == :esa ? "ESA" : "GSB")")
    println(io, "  τ = $(round(r.fit.tau, digits=2))")
    println(io, "  σ_IRF = $(round(r.fit.sigma, digits=2))")
    println(io, "  R² = $(round(r.fit.rsquared, digits=4))")
end

# =============================================================================
# Peak fitting result types
# =============================================================================

"""
    PeakFitResult

Result of peak fitting with any lineshape model.

Access parameters by name: `result[:center].value`, `result[:fwhm].err`
"""
struct PeakFitResult
    params::Vector{Symbol}
    values::Vector{Float64}
    errors::Vector{Float64}
    ci::Vector{Tuple{Float64,Float64}}
    r_squared::Float64
    rss::Float64
    mse::Float64
    npoints::Int
    region::Tuple{Float64,Float64}
    model::String
    sample_id::String
end

# Index by parameter name
function Base.getindex(r::PeakFitResult, name::Symbol)
    idx = findfirst(==(name), r.params)
    isnothing(idx) && error("Parameter :$name not found. Available: $(r.params)")
    return (value=r.values[idx], err=r.errors[idx], ci=r.ci[idx])
end

Base.haskey(r::PeakFitResult, name::Symbol) = name in r.params

function Base.show(io::IO, r::PeakFitResult)
    println(io, "PeakFitResult: $(r.sample_id)")
    println(io, "  Region: $(round(Int, r.region[1])) - $(round(Int, r.region[2])) cm⁻¹ ($(r.npoints) points)")
    println(io, "  Model:  $(r.model)")
    println(io)

    println(io, "  Parameter       Value       Std Err     95% CI")
    println(io, "  ─────────────────────────────────────────────────────────")

    function fmt_value(x)
        if abs(x) < 0.01 || abs(x) >= 10000
            return rpad(string(round(x, sigdigits=4)), 10)
        else
            return rpad(string(round(x, digits=2)), 10)
        end
    end
    fmt_ci(ci) = "($(round(ci[1], digits=2)), $(round(ci[2], digits=2)))"

    for (i, name) in enumerate(r.params)
        label = rpad(string(name), 14)
        println(io, "  $(label)  $(fmt_value(r.values[i]))  $(fmt_value(r.errors[i]))  $(fmt_ci(r.ci[i]))")
    end

    println(io)
    print(io, "  R² = $(round(r.r_squared, digits=5))    MSE = $(round(r.mse, sigdigits=4))")
end

# =============================================================================
# Multi-peak fitting result
# =============================================================================

"""
    MultiPeakFitResult

Result of multi-peak fitting (1 to N peaks with polynomial baseline).

Supports indexing by peak number: `result[1]` returns first peak's `PeakFitResult`.
"""
struct MultiPeakFitResult
    peaks::Vector{PeakFitResult}
    baseline_params::Vector{Float64}
    baseline_order::Int
    r_squared::Float64
    rss::Float64
    mse::Float64
    npoints::Int
    region::Tuple{Float64,Float64}
    sample_id::String
    _coef::Vector{Float64}
    _peak_fn::Function
    _n_peak_params::Int
    _x::Vector{Float64}
    _y::Vector{Float64}
end

Base.getindex(r::MultiPeakFitResult, i::Int) = r.peaks[i]
Base.length(r::MultiPeakFitResult) = length(r.peaks)
Base.iterate(r::MultiPeakFitResult) = iterate(r.peaks)
Base.iterate(r::MultiPeakFitResult, state) = iterate(r.peaks, state)
Base.firstindex(r::MultiPeakFitResult) = 1
Base.lastindex(r::MultiPeakFitResult) = length(r.peaks)

function Base.show(io::IO, r::MultiPeakFitResult)
    n = length(r.peaks)
    model = n > 0 ? r.peaks[1].model : "unknown"
    print(io, "MultiPeakFitResult: $n peak$(n == 1 ? "" : "s") ($model), R² = $(round(r.r_squared, digits=4))")
end

function Base.show(io::IO, ::MIME"text/plain", r::MultiPeakFitResult)
    n = length(r.peaks)
    model = n > 0 ? r.peaks[1].model : "unknown"

    println(io, "MultiPeakFitResult ($n peak$(n == 1 ? "" : "s"), $model)")
    println(io, "  Region: $(round(Int, r.region[1])) - $(round(Int, r.region[2])) ($(r.npoints) points)")
    if !isempty(r.sample_id)
        println(io, "  Sample: $(r.sample_id)")
    end
    println(io)

    if n > 0
        param_names = r.peaks[1].params
        header = "  Peak   " * join([rpad(string(p), 14) for p in param_names], "")
        println(io, header)
        println(io, "  " * "─"^(length(header) - 2))

        for (i, pk) in enumerate(r.peaks)
            row = "  $(lpad(i, 4))   "
            for (j, _) in enumerate(pk.params)
                v = round(pk.values[j], digits=2)
                row *= rpad(string(v), 14)
            end
            println(io, row)
        end
    end

    println(io)
    if r.baseline_order == 0
        println(io, "  Baseline: constant = $(round(r.baseline_params[1], sigdigits=4))")
    elseif r.baseline_order == 1
        println(io, "  Baseline: linear, c₀ = $(round(r.baseline_params[1], sigdigits=4)), c₁ = $(round(r.baseline_params[2], sigdigits=4))")
    else
        println(io, "  Baseline: order $(r.baseline_order), params = $(round.(r.baseline_params, sigdigits=4))")
    end

    println(io)
    print(io, "  R² = $(round(r.r_squared, digits=5))    MSE = $(round(r.mse, sigdigits=4))")
end

# =============================================================================
# TAMatrix: 2D Transient Absorption Data
# =============================================================================

"""
    TAMatrix <: AbstractSpectroscopyData

Two-dimensional transient absorption data (time x wavelength).

# Fields
- `time::Vector{Float64}`: Time axis (ps)
- `wavelength::Vector{Float64}`: Wavelength axis (nm) or wavenumber (cm⁻¹)
- `data::Matrix{Float64}`: ΔA signal matrix, size (n_time, n_wavelength)
- `metadata::Dict{Symbol,Any}`: Additional info

# Indexing
```julia
matrix[λ=800]     # Extract TATrace at λ ≈ 800 nm
matrix[t=1.0]     # Extract TASpectrum at t ≈ 1.0 ps
```
"""
struct TAMatrix <: AbstractSpectroscopyData
    time::Vector{Float64}
    wavelength::Vector{Float64}
    data::Matrix{Float64}
    metadata::Dict{Symbol,Any}
end

TAMatrix(time, wavelength, data; metadata=Dict{Symbol,Any}()) =
    TAMatrix(time, wavelength, data, metadata)

xdata(m::TAMatrix) = m.wavelength
ydata(m::TAMatrix) = m.time
zdata(m::TAMatrix) = m.data
is_matrix(::TAMatrix) = true
source_file(m::TAMatrix) = get(m.metadata, :source, "")
npoints(m::TAMatrix) = size(m.data)

function xlabel(m::TAMatrix)
    minval, maxval = extrema(m.wavelength)
    if minval > 1200 && maxval < 5000
        return "Wavenumber (cm⁻¹)"
    else
        return "Wavelength (nm)"
    end
end
ylabel(::TAMatrix) = "Time (ps)"
zlabel(::TAMatrix) = "ΔA"

function _detect_wavelength_unit(m::TAMatrix)
    minval, maxval = extrema(m.wavelength)
    (minval > 1200 && maxval < 5000) ? "cm⁻¹" : "nm"
end

function Base.show(io::IO, m::TAMatrix)
    n_time, n_wl = size(m.data)
    t_range = "$(round(minimum(m.time), digits=2)) to $(round(maximum(m.time), digits=2)) ps"
    wl_range = "$(round(minimum(m.wavelength), digits=1)) to $(round(maximum(m.wavelength), digits=1))"
    wl_unit = _detect_wavelength_unit(m)
    print(io, "TAMatrix: $n_time × $n_wl ($t_range, $wl_range $wl_unit)")
end

function Base.show(io::IO, ::MIME"text/plain", m::TAMatrix)
    n_time, n_wl = size(m.data)
    wl_unit = _detect_wavelength_unit(m)

    println(io, "TAMatrix")
    println(io, "  Time points:   $n_time ($(round(minimum(m.time), digits=2)) to $(round(maximum(m.time), digits=2)) ps)")
    println(io, "  Wavelengths:   $n_wl ($(round(minimum(m.wavelength), digits=1)) to $(round(maximum(m.wavelength), digits=1)) $wl_unit)")
    println(io, "  Data range:    $(round(minimum(m.data), sigdigits=3)) to $(round(maximum(m.data), sigdigits=3))")
    if haskey(m.metadata, :source)
        println(io, "  Source:        $(m.metadata[:source])")
    end
end

# =============================================================================
# TAMatrix indexing
# =============================================================================

function _find_nearest_idx(arr::AbstractVector, target::Real)
    _, idx = findmin(abs.(arr .- target))
    return idx
end

function Base.getindex(m::TAMatrix; λ=nothing, t=nothing)
    if !isnothing(λ) && isnothing(t)
        idx = _find_nearest_idx(m.wavelength, λ)
        actual_λ = m.wavelength[idx]
        signal = m.data[:, idx]

        metadata = Dict{Symbol,Any}(
            :extracted_from => get(m.metadata, :source, "TAMatrix"),
            :requested_wavelength => λ,
            :actual_wavelength => actual_λ,
            :wavelength_index => idx
        )

        return TATrace(m.time, signal, actual_λ, metadata)

    elseif !isnothing(t) && isnothing(λ)
        idx = _find_nearest_idx(m.time, t)
        actual_t = m.time[idx]
        signal = m.data[idx, :]

        metadata = Dict{Symbol,Any}(
            :extracted_from => get(m.metadata, :source, "TAMatrix"),
            :requested_time => t,
            :actual_time => actual_t,
            :time_index => idx
        )

        return TASpectrum(m.wavelength, signal, actual_t, metadata)

    elseif !isnothing(λ) && !isnothing(t)
        error("Cannot specify both λ and t. Use matrix[λ=...] or matrix[t=...]")
    else
        error("Must specify either λ or t for indexing. Use matrix[λ=...] or matrix[t=...]")
    end
end

# =============================================================================
# Global fitting result
# =============================================================================

"""
    GlobalFitResult

Result of global fitting multiple traces with shared parameters.

# Fields
- `tau`, `sigma`, `t0`: Shared parameters
- `amplitudes`, `offsets`, `labels`: Per-trace parameters
- `rsquared`, `rsquared_individual`, `residuals`: Fit quality
"""
struct GlobalFitResult
    tau::Float64
    sigma::Float64
    t0::Float64
    amplitudes::Vector{Float64}
    offsets::Vector{Float64}
    labels::Vector{String}
    rsquared::Float64
    rsquared_individual::Vector{Float64}
    residuals::Vector{Vector{Float64}}
end

function Base.show(io::IO, r::GlobalFitResult)
    print(io, "GlobalFitResult: τ = $(round(r.tau, digits=2)) ps, R² = $(round(r.rsquared, digits=4)) ($(length(r.amplitudes)) traces)")
end

function Base.show(io::IO, ::MIME"text/plain", r::GlobalFitResult)
    n = length(r.amplitudes)

    println(io, "GlobalFitResult ($(n) traces)")
    println(io)

    println(io, "  Shared Parameters")
    println(io, "  ─────────────────────────────")
    println(io, "  τ        $(lpad(round(r.tau, digits=3), 8)) ps")
    println(io, "  σ_IRF    $(lpad(round(r.sigma, digits=3), 8)) ps")
    println(io, "  t₀       $(lpad(round(r.t0, digits=3), 8)) ps")
    println(io)

    println(io, "  Trace         Amplitude      Offset       R²")
    println(io, "  ─────────────────────────────────────────────────")
    for i in 1:n
        label = rpad(r.labels[i], 12)
        amp = lpad(round(r.amplitudes[i], sigdigits=4), 12)
        off = lpad(round(r.offsets[i], sigdigits=4), 12)
        r2 = round(r.rsquared_individual[i], digits=4)
        println(io, "  $label $amp $off    $r2")
    end

    println(io)
    print(io, "  Global R² = $(round(r.rsquared, digits=5))")
end

# =============================================================================
# report() and format_results()
# =============================================================================

"""
    report(result)

Print a formatted summary of a fit result to the terminal.
"""
report(result) = (show(stdout, MIME("text/plain"), result); println(); nothing)

"""
    format_results(result) -> String

Return a Markdown-formatted string of fit results.
"""
function format_results end

function format_results(r::MultiPeakFitResult)
    io = IOBuffer()
    n = length(r.peaks)
    model = n > 0 ? r.peaks[1].model : "unknown"

    println(io, "## Peak Fit Results")
    println(io)

    for (i, pk) in enumerate(r.peaks)
        if n > 1
            println(io, "### Peak $i")
            println(io)
        end
        println(io, "| Parameter | Value | Uncertainty |")
        println(io, "|-----------|-------|-------------|")
        for (j, name) in enumerate(pk.params)
            val = round(pk.values[j], digits=4)
            err = round(pk.errors[j], digits=4)
            println(io, "| $(name) | $(val) | ± $(err) |")
        end
        println(io)
    end

    if r.baseline_order == 0
        println(io, "**Baseline:** constant = $(round(r.baseline_params[1], sigdigits=4))")
    elseif r.baseline_order == 1
        println(io, "**Baseline:** linear")
    else
        println(io, "**Baseline:** polynomial order $(r.baseline_order)")
    end

    println(io, "**Model:** $(model) | **R²:** $(round(r.r_squared, digits=5)) | **Region:** $(round(Int, r.region[1]))–$(round(Int, r.region[2]))")

    return String(take!(io))
end

function format_results(r::ExpDecayIRFFit)
    io = IOBuffer()
    signal_str = r.signal_type == :esa ? "ESA" : "GSB"
    has_irf = !isnan(r.sigma)

    println(io, "## Exponential Decay Fit ($signal_str)")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "|-----------|-------|")
    println(io, "| τ | $(round(r.tau, digits=3)) ps |")
    if has_irf
        println(io, "| σ_IRF | $(round(r.sigma, digits=3)) ps |")
    end
    println(io, "| t₀ | $(round(r.t0, digits=3)) ps |")
    println(io, "| Amplitude | $(round(r.amplitude, sigdigits=4)) |")
    println(io, "| Offset | $(round(r.offset, sigdigits=4)) |")
    println(io)
    println(io, "**R²:** $(round(r.rsquared, digits=5))")

    return String(take!(io))
end

function format_results(r::ExpDecayFit)
    io = IOBuffer()
    signal_str = r.signal_type == :esa ? "ESA" : "GSB"

    println(io, "## Exponential Decay Fit ($signal_str)")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "|-----------|-------|")
    println(io, "| τ | $(round(r.tau, digits=3)) |")
    println(io, "| Amplitude | $(round(r.amplitude, sigdigits=4)) |")
    println(io, "| Offset | $(round(r.offset, sigdigits=4)) |")
    println(io)
    println(io, "**R²:** $(round(r.rsquared, digits=5))")

    return String(take!(io))
end

function format_results(r::BiexpDecayFit)
    io = IOBuffer()
    signal_str = r.signal_type == :esa ? "ESA" : "GSB"
    has_irf = !isnan(r.sigma)

    total_amp = abs(r.amplitude1) + abs(r.amplitude2)
    w1 = total_amp > 0 ? round(100 * abs(r.amplitude1) / total_amp, digits=1) : 0.0
    w2 = total_amp > 0 ? round(100 * abs(r.amplitude2) / total_amp, digits=1) : 0.0

    println(io, "## Biexponential Decay Fit ($signal_str)")
    println(io)
    println(io, "| Component | τ (ps) | Amplitude | Weight |")
    println(io, "|-----------|--------|-----------|--------|")
    println(io, "| Fast | $(round(r.tau1, digits=3)) | $(round(r.amplitude1, sigdigits=4)) | $(w1)% |")
    println(io, "| Slow | $(round(r.tau2, digits=3)) | $(round(r.amplitude2, sigdigits=4)) | $(w2)% |")
    println(io)
    if has_irf
        println(io, "**σ_IRF:** $(round(r.sigma, digits=3)) ps | **t₀:** $(round(r.t0, digits=3)) ps | **Offset:** $(round(r.offset, sigdigits=4))")
    else
        println(io, "**t₀:** $(round(r.t0, digits=3)) ps | **Offset:** $(round(r.offset, sigdigits=4))")
    end
    println(io, "**R²:** $(round(r.rsquared, digits=5))")

    return String(take!(io))
end

function format_results(r::MultiexpDecayFit)
    io = IOBuffer()
    n = length(r.taus)
    signal_str = r.signal_type == :esa ? "ESA" : "GSB"
    has_irf = !isnan(r.sigma)
    w = weights(r)

    println(io, "## Multi-exponential Decay Fit ($signal_str, $n components)")
    println(io)
    println(io, "| Component | τ (ps) | Amplitude | Weight |")
    println(io, "|-----------|--------|-----------|--------|")
    for i in 1:n
        println(io, "| $i | $(round(r.taus[i], digits=3)) | $(round(r.amplitudes[i], sigdigits=4)) | $(round(100*w[i], digits=1))% |")
    end
    println(io)
    if has_irf
        println(io, "**σ_IRF:** $(round(r.sigma, digits=3)) ps | **t₀:** $(round(r.t0, digits=3)) ps | **Offset:** $(round(r.offset, sigdigits=4))")
    else
        println(io, "**t₀:** $(round(r.t0, digits=3)) ps | **Offset:** $(round(r.offset, sigdigits=4))")
    end
    println(io, "**R²:** $(round(r.rsquared, digits=5))")

    return String(take!(io))
end

function format_results(r::GlobalFitResult)
    io = IOBuffer()
    n = length(r.amplitudes)

    println(io, "## Global Fit Results ($n traces)")
    println(io)
    println(io, "### Shared Parameters")
    println(io)
    println(io, "| Parameter | Value |")
    println(io, "|-----------|-------|")
    println(io, "| τ | $(round(r.tau, digits=3)) ps |")
    println(io, "| σ_IRF | $(round(r.sigma, digits=3)) ps |")
    println(io, "| t₀ | $(round(r.t0, digits=3)) ps |")
    println(io)
    println(io, "### Per-Trace Parameters")
    println(io)
    println(io, "| Trace | Amplitude | Offset | R² |")
    println(io, "|-------|-----------|--------|----|")
    for i in 1:n
        label = r.labels[i]
        amp = round(r.amplitudes[i], sigdigits=4)
        off = round(r.offsets[i], sigdigits=4)
        r2 = round(r.rsquared_individual[i], digits=4)
        println(io, "| $(label) | $(amp) | $(off) | $(r2) |")
    end
    println(io)
    println(io, "**Global R²:** $(round(r.rsquared, digits=5))")

    return String(take!(io))
end
