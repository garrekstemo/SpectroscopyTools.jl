# PL/Raman spatial mapping types and analysis

# =============================================================================
# PLMap type
# =============================================================================

"""
    PLMap <: AbstractSpectroscopyData

Photoluminescence intensity map from a CCD raster scan.

A 2D spatial grid where each point has a full CCD spectrum. The `intensity`
field holds the integrated PL signal at each position; the full spectra are
stored in `spectra` for extraction at individual positions.

# Fields
- `intensity::Matrix{Float64}` — Integrated PL intensity `(nx, ny)`
- `spectra::Array{Float64,3}` — Raw CCD counts `(nx, ny, n_pixel)`
- `x::Vector{Float64}` — Spatial x positions (μm)
- `y::Vector{Float64}` — Spatial y positions (μm)
- `pixel::Vector{Float64}` — Pixel indices (or wavelength if calibrated)
- `metadata::Dict{String,Any}` — Source file, grid dims, step size, etc.
"""
struct PLMap <: AbstractSpectroscopyData
    intensity::Matrix{Float64}
    spectra::Array{Float64,3}
    x::Vector{Float64}
    y::Vector{Float64}
    pixel::Vector{Float64}
    metadata::Dict{String,Any}
end

# =============================================================================
# AbstractSpectroscopyData interface
# =============================================================================

xdata(m::PLMap) = m.x
ydata(m::PLMap) = m.y
zdata(m::PLMap) = m.intensity
xlabel(::PLMap) = "X (μm)"
ylabel(::PLMap) = "Y (μm)"
zlabel(::PLMap) = "PL Intensity"
is_matrix(::PLMap) = true
npoints(m::PLMap) = (length(m.x), length(m.y))
source_file(m::PLMap) = get(m.metadata, "source_file", "unknown")
title(m::PLMap) = source_file(m)

# Semantic accessor
"""
    intensity(m::PLMap) -> Matrix{Float64}

Return the PL intensity matrix.
"""
intensity(m::PLMap) = m.intensity

function Base.show(io::IO, m::PLMap)
    nx, ny = length(m.x), length(m.y)
    np = length(m.pixel)
    print(io, "PLMap($(nx)×$(ny) grid, $(np) pixels)")
end

function Base.show(io::IO, ::MIME"text/plain", m::PLMap)
    nx, ny = length(m.x), length(m.y)
    np = length(m.pixel)
    println(io, "PLMap")
    println(io, "  Grid:     $(nx) × $(ny) spatial points")
    println(io, "  Pixels:   $(np) per spectrum")
    println(io, "  X range:  $(round(m.x[1], digits=1)) to $(round(m.x[end], digits=1)) μm")
    println(io, "  Y range:  $(round(m.y[1], digits=1)) to $(round(m.y[end], digits=1)) μm")
    print(io, "  Source:   $(source_file(m))")
end

# =============================================================================
# Spectrum extraction
# =============================================================================

"""
    extract_spectrum(m::PLMap, ix::Int, iy::Int) -> NamedTuple

Extract the CCD spectrum at grid index `(ix, iy)`.

Returns `(pixel=..., signal=..., x=..., y=...)`.
"""
function extract_spectrum(m::PLMap, ix::Int, iy::Int)
    1 <= ix <= length(m.x) || error("ix=$ix out of range 1:$(length(m.x))")
    1 <= iy <= length(m.y) || error("iy=$iy out of range 1:$(length(m.y))")
    return (pixel=m.pixel, signal=vec(m.spectra[ix, iy, :]),
            x=m.x[ix], y=m.y[iy])
end

"""
    extract_spectrum(m::PLMap; x::Real, y::Real) -> NamedTuple

Extract the CCD spectrum at the spatial position nearest to `(x, y)`.

Returns `(pixel=..., signal=..., x=..., y=..., ix=..., iy=...)`.
"""
function extract_spectrum(m::PLMap; x::Real, y::Real)
    ix = argmin(abs.(m.x .- x))
    iy = argmin(abs.(m.y .- y))
    spec = extract_spectrum(m, ix, iy)
    return (pixel=spec.pixel, signal=spec.signal,
            x=spec.x, y=spec.y, ix=ix, iy=iy)
end

# =============================================================================
# Background subtraction
# =============================================================================

"""
    subtract_background(m::PLMap; positions=nothing, margin=5) -> PLMap

Subtract a background spectrum from every grid point.

The background is the average CCD spectrum over the reference positions.
After subtraction, the intensity map is recomputed from the corrected spectra
using the same `pixel_range` as the original load (if any).

# Arguments
- `positions`: Vector of `(x, y)` spatial coordinate tuples (μm) for background
  reference points. These should be off-flake positions with no PL signal.
- `margin`: Number of grid points from each edge used for auto-detection when
  `positions` is not given. Auto mode averages the corners of the bottom half
  of the map (avoids top-row artifacts). Default: 5.

# Example
```julia
m = load_pl_map("data.lvm"; nx=51, ny=51, step_size=2.16, pixel_range=(950, 1100))

# Explicit background positions
m_bg = subtract_background(m; positions=[(-40, -40), (40, -40), (-40, -20)])

# Auto mode (bottom corners)
m_bg = subtract_background(m)
```
"""
function subtract_background(m::PLMap; positions=nothing, margin::Int=5)
    nx, ny = length(m.x), length(m.y)

    if !isnothing(positions)
        # Explicit: average spectra at user-specified (x, y) positions
        bg_spectra = zeros(length(m.pixel))
        bg_positions = Tuple{Float64,Float64}[]
        for pos in positions
            spec = extract_spectrum(m; x=pos[1], y=pos[2])
            bg_spectra .+= spec.signal
            push!(bg_positions, (spec.x, spec.y))
        end
        bg_spectra ./= length(positions)
    else
        # Auto: average corners from the bottom half of the map
        bg_spectra = zeros(length(m.pixel))
        bg_positions = Tuple{Float64,Float64}[]
        for ix in [1:margin; (nx-margin+1):nx]
            for iy in 1:margin
                bg_spectra .+= vec(m.spectra[ix, iy, :])
                push!(bg_positions, (m.x[ix], m.y[iy]))
            end
        end
        bg_spectra ./= length(bg_positions)
        @info "Auto background: averaged $(length(bg_positions)) spectra from bottom corners " *
              "(ix=1:$margin and $(nx-margin+1):$nx, iy=1:$margin)"
    end

    # Subtract background spectrum from every grid point
    corrected = m.spectra .- reshape(bg_spectra, 1, 1, :)

    # Recompute intensity with the same pixel_range as the original
    pixel_range = get(m.metadata, "pixel_range", nothing)
    if !isnothing(pixel_range)
        p1, p2 = pixel_range
        intensity = dropdims(sum(corrected[:, :, p1:p2]; dims=3); dims=3)
    else
        intensity = dropdims(sum(corrected; dims=3); dims=3)
    end

    new_metadata = copy(m.metadata)
    new_metadata["background_positions"] = bg_positions
    return PLMap(intensity, corrected, m.x, m.y, m.pixel, new_metadata)
end

# =============================================================================
# Normalization
# =============================================================================

"""
    normalize(m::PLMap) -> PLMap

Return a new PLMap with intensity normalized to [0, 1].
"""
function normalize(m::PLMap)
    imin, imax = extrema(m.intensity)
    if imax == imin
        norm_intensity = zeros(size(m.intensity))
    else
        norm_intensity = (m.intensity .- imin) ./ (imax - imin)
    end
    return PLMap(norm_intensity, m.spectra, m.x, m.y, m.pixel, m.metadata)
end

# =============================================================================
# Peak center (centroid) map
# =============================================================================

"""
    peak_centers(m::PLMap; pixel_range=nothing, threshold=0.05) -> Matrix{Float64}

Compute the centroid (intensity-weighted average pixel) at each grid point.

Returns a `(nx, ny)` matrix of peak center positions in pixel units. Grid points
where the PL intensity (`m.intensity`) is below `threshold` × the map maximum
are set to `NaN` (renders as transparent in heatmaps with `nan_color=:transparent`).

Masking uses the PLMap's intensity field — the integrated PL signal already stored
in the map. This produces clean masks that match the intensity heatmap: off-flake
regions (low PL signal) are transparent, on-flake regions show centroid positions.

# Arguments
- `pixel_range`: `(start, stop)` pixel range to compute the centroid over.
  Falls back to the `pixel_range` stored in metadata, or all pixels if unset.
- `threshold`: Fraction of the maximum PL intensity below which a grid point is
  masked as `NaN`. Default `0.05` (5%). Set to `0` to disable masking.

# Example
```julia
m = load_pl_map("scan.lvm"; nx=51, ny=51, pixel_range=(950, 1100))
m = subtract_background(m)
centers = peak_centers(m)
heatmap(m.x, m.y, centers; colormap=:viridis, nan_color=:transparent)
```
"""
function peak_centers(m::PLMap; pixel_range::Union{Tuple{Int,Int},Nothing}=nothing,
                      threshold::Real=0.05)
    pr = !isnothing(pixel_range) ? pixel_range : get(m.metadata, "pixel_range", nothing)

    if !isnothing(pr)
        p1 = max(1, pr[1])
        p2 = min(length(m.pixel), pr[2])
        pixels = m.pixel[p1:p2]
        spectra_slice = @view m.spectra[:, :, p1:p2]
    else
        pixels = m.pixel
        spectra_slice = m.spectra
    end

    nx, ny = length(m.x), length(m.y)

    # Mask against the PL intensity map (integrated signal already in m.intensity)
    max_intensity = maximum(m.intensity)
    cutoff = max_intensity > 0 ? max_intensity * threshold : zero(max_intensity)

    centers = Matrix{Float64}(undef, nx, ny)
    for iy in 1:ny
        for ix in 1:nx
            if m.intensity[ix, iy] <= cutoff
                centers[ix, iy] = NaN
            else
                sig = @view spectra_slice[ix, iy, :]
                total = sum(sig)
                centers[ix, iy] = sum(pixels[k] * sig[k] for k in eachindex(sig)) / total
            end
        end
    end
    return centers
end
