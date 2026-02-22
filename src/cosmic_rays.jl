# Cosmic ray detection and removal
#
# Cosmic rays are single-pixel spectral artifacts from high-energy particles
# hitting the CCD detector. They appear as sharp, narrow spikes that don't
# correlate with spatial neighbors.
#
# Algorithm: Whitaker-Hayes modified z-score on first differences.
# For PLMap data, spatial validation removes false positives where
# neighboring pixels share the same spectral feature.

# =============================================================================
# Result types
# =============================================================================

"""
    CosmicRayResult

Result of cosmic ray detection on a 1D signal.

# Fields
- `indices::Vector{Int}` — flagged channel indices
- `count::Int` — number of flagged channels
"""
struct CosmicRayResult
    indices::Vector{Int}
    count::Int
end

"""
    CosmicRayMapResult

Result of cosmic ray detection on a PLMap.

# Fields
- `mask::BitArray{3}` — `(nx, ny, n_pixel)`, `true` = cosmic ray
- `count::Int` — total flagged voxels
- `affected_spectra::Int` — number of spectra with at least one cosmic ray
- `channel_counts::Vector{Int}` — cosmic ray count per spectral channel
"""
struct CosmicRayMapResult
    mask::BitArray{3}
    count::Int
    affected_spectra::Int
    channel_counts::Vector{Int}
end

# =============================================================================
# 1D detection and removal
# =============================================================================

"""
    detect_cosmic_rays(signal::AbstractVector; threshold=5.0) -> CosmicRayResult

Detect cosmic ray spikes in a 1D spectrum using modified z-scores on first
differences (Whitaker-Hayes method).

Computes first differences `Δ[k] = signal[k+1] - signal[k]`, then flags
channels where the modified z-score exceeds `threshold`. The modified z-score
uses the Median Absolute Deviation (MAD) for robust scale estimation.

# Arguments
- `signal`: 1D spectral signal
- `threshold`: z-score cutoff (default 5.0). Lower values detect more spikes
  but may flag real features.

# Returns
A [`CosmicRayResult`](@ref) with the flagged channel indices.

# Example
```julia
result = detect_cosmic_rays(spectrum; threshold=5.0)
println("Found \$(result.count) cosmic rays at indices \$(result.indices)")
```
"""
function detect_cosmic_rays(signal::AbstractVector; threshold::Real=5.0)
    n = length(signal)
    if n < 3
        return CosmicRayResult(Int[], 0)
    end

    # First differences
    diffs = diff(signal)

    med = median(diffs)
    mad_val = median(abs.(diffs .- med))

    if mad_val < eps(Float64)
        return CosmicRayResult(Int[], 0)
    end

    # Modified z-scores (0.6745 ≈ Φ⁻¹(0.75), makes MAD consistent with σ)
    z_scores = @. 0.6745 * (diffs - med) / mad_val

    # A cosmic ray spike creates a positive jump followed by a negative jump.
    # Flag the channel where the spike sits: if z[k] > threshold (upward jump),
    # the spike is at k+1; if z[k] < -threshold (downward jump from previous
    # spike), the spike was at k. We collect unique indices.
    flagged = Set{Int}()
    for k in eachindex(z_scores)
        if z_scores[k] > threshold
            push!(flagged, k + 1)  # spike at the "to" channel
        elseif z_scores[k] < -threshold
            push!(flagged, k)      # spike at the "from" channel
        end
    end

    indices = sort!(collect(flagged))
    # Clamp to valid range
    filter!(i -> 1 <= i <= n, indices)

    return CosmicRayResult(indices, length(indices))
end

"""
    remove_cosmic_rays(signal::AbstractVector, result::CosmicRayResult) -> Vector

Replace flagged cosmic ray channels with linearly interpolated values from
the nearest non-flagged neighbors.

Returns a new vector; does not mutate the input.

# Example
```julia
result = detect_cosmic_rays(spectrum)
cleaned = remove_cosmic_rays(spectrum, result)
```
"""
function remove_cosmic_rays(signal::AbstractVector, result::CosmicRayResult)
    cleaned = collect(Float64, signal)
    if result.count == 0
        return cleaned
    end

    n = length(cleaned)
    flagged_set = Set(result.indices)

    for idx in result.indices
        # Find nearest non-flagged neighbor to the left
        left_idx = idx - 1
        while left_idx >= 1 && left_idx in flagged_set
            left_idx -= 1
        end

        # Find nearest non-flagged neighbor to the right
        right_idx = idx + 1
        while right_idx <= n && right_idx in flagged_set
            right_idx += 1
        end

        if left_idx >= 1 && right_idx <= n
            # Linear interpolation
            t = (idx - left_idx) / (right_idx - left_idx)
            cleaned[idx] = cleaned[left_idx] * (1 - t) + cleaned[right_idx] * t
        elseif left_idx >= 1
            cleaned[idx] = cleaned[left_idx]
        elseif right_idx <= n
            cleaned[idx] = cleaned[right_idx]
        end
        # If both out of range (entire signal flagged), leave as-is
    end

    return cleaned
end

# =============================================================================
# PLMap detection and removal
# =============================================================================

"""
    detect_cosmic_rays(m::PLMap; threshold=5.0, pixel_range=nothing) -> CosmicRayMapResult

Detect cosmic ray spikes across all spectra in a PLMap.

Runs 1D detection on each spatial pixel's spectrum, then applies spatial
validation: if a flagged channel `k` at position `(ix, iy)` is also flagged
at the same channel in 2 or more 4-connected spatial neighbors, it is unmarked
as a real spectral feature rather than a cosmic ray.

# Arguments
- `m`: PLMap with 3D spectra array `(nx, ny, n_pixel)`
- `threshold`: z-score cutoff (default 5.0)
- `pixel_range`: `(start, stop)` channel indices to analyze. When set, only
  this subrange of each spectrum is checked for cosmic rays. Channels outside
  the range are never flagged. Defaults to the `pixel_range` in `m.metadata`,
  or the full spectrum if unset.

# Returns
A [`CosmicRayMapResult`](@ref) with a 3D boolean mask and summary statistics.

# Example
```julia
cr = detect_cosmic_rays(plmap; threshold=5.0)
println("Found \$(cr.count) cosmic rays in \$(cr.affected_spectra) spectra")
```
"""
function detect_cosmic_rays(m::PLMap; threshold::Real=5.0,
                            pixel_range::Union{Tuple{Int,Int},Nothing}=nothing)
    nx, ny, np = size(m.spectra)
    mask = falses(nx, ny, np)

    # Determine channel range for detection
    pr = !isnothing(pixel_range) ? pixel_range : get(m.metadata, "pixel_range", nothing)
    p1 = !isnothing(pr) ? max(1, Int(pr[1])) : 1
    p2 = !isnothing(pr) ? min(np, Int(pr[2])) : np

    # Step 1: Run 1D detection on each spatial pixel (within pixel_range only)
    for iy in 1:ny
        for ix in 1:nx
            signal = @view m.spectra[ix, iy, p1:p2]
            cr = detect_cosmic_rays(signal; threshold=Float64(threshold))
            for k in cr.indices
                mask[ix, iy, k + p1 - 1] = true  # offset back to full-spectrum index
            end
        end
    end

    # Step 2: Spatial validation — unmark if ≥2 neighbors also flagged at same channel
    # 4-connected neighbors: (ix±1, iy), (ix, iy±1)
    to_unmark = Tuple{Int,Int,Int}[]
    for iy in 1:ny
        for ix in 1:nx
            for k in 1:np
                mask[ix, iy, k] || continue

                neighbor_count = 0
                if ix > 1 && mask[ix-1, iy, k]
                    neighbor_count += 1
                end
                if ix < nx && mask[ix+1, iy, k]
                    neighbor_count += 1
                end
                if iy > 1 && mask[ix, iy-1, k]
                    neighbor_count += 1
                end
                if iy < ny && mask[ix, iy+1, k]
                    neighbor_count += 1
                end

                if neighbor_count >= 2
                    push!(to_unmark, (ix, iy, k))
                end
            end
        end
    end

    for (ix, iy, k) in to_unmark
        mask[ix, iy, k] = false
        # Also unmark the neighbors at the same channel (they're part of the same feature)
        if ix > 1 && mask[ix-1, iy, k]
            mask[ix-1, iy, k] = false
        end
        if ix < nx && mask[ix+1, iy, k]
            mask[ix+1, iy, k] = false
        end
        if iy > 1 && mask[ix, iy-1, k]
            mask[ix, iy-1, k] = false
        end
        if iy < ny && mask[ix, iy+1, k]
            mask[ix, iy+1, k] = false
        end
    end

    total = count(mask)
    affected = 0
    for iy in 1:ny
        for ix in 1:nx
            if any(@view mask[ix, iy, :])
                affected += 1
            end
        end
    end

    channel_counts = zeros(Int, np)
    for k in 1:np
        channel_counts[k] = count(@view mask[:, :, k])
    end

    return CosmicRayMapResult(mask, total, affected, channel_counts)
end

"""
    remove_cosmic_rays(m::PLMap, result::CosmicRayMapResult) -> PLMap

Remove cosmic rays from a PLMap by replacing flagged voxels.

For each flagged `(ix, iy, k)`: replaces with the median of `spectra[neighbors, k]`
from non-flagged 4-connected spatial neighbors. Falls back to spectral interpolation
(linear from nearest non-flagged channels) for edge pixels with no valid neighbors.

Returns a new PLMap with recomputed intensity.

# Example
```julia
cr = detect_cosmic_rays(plmap)
cleaned = remove_cosmic_rays(plmap, cr)
```
"""
function remove_cosmic_rays(m::PLMap, result::CosmicRayMapResult)
    if result.count == 0
        return m
    end

    nx, ny, np = size(m.spectra)
    cleaned = copy(m.spectra)

    for iy in 1:ny
        for ix in 1:nx
            for k in 1:np
                result.mask[ix, iy, k] || continue

                # Collect valid neighbor values at channel k
                neighbor_vals = Float64[]
                if ix > 1 && !result.mask[ix-1, iy, k]
                    push!(neighbor_vals, m.spectra[ix-1, iy, k])
                end
                if ix < nx && !result.mask[ix+1, iy, k]
                    push!(neighbor_vals, m.spectra[ix+1, iy, k])
                end
                if iy > 1 && !result.mask[ix, iy-1, k]
                    push!(neighbor_vals, m.spectra[ix, iy-1, k])
                end
                if iy < ny && !result.mask[ix, iy+1, k]
                    push!(neighbor_vals, m.spectra[ix, iy+1, k])
                end

                if !isempty(neighbor_vals)
                    cleaned[ix, iy, k] = median(neighbor_vals)
                else
                    # Fallback: spectral interpolation from the mask
                    flagged_channels = findall(@view result.mask[ix, iy, :])
                    cr_1d = CosmicRayResult(flagged_channels, length(flagged_channels))
                    signal = @view m.spectra[ix, iy, :]
                    full_cleaned = remove_cosmic_rays(signal, cr_1d)
                    cleaned[ix, iy, k] = full_cleaned[k]
                end
            end
        end
    end

    # Recompute intensity
    pixel_range = get(m.metadata, "pixel_range", nothing)
    if !isnothing(pixel_range)
        p1, p2 = pixel_range
        new_intensity = dropdims(sum(cleaned[:, :, p1:p2]; dims=3); dims=3)
    else
        new_intensity = dropdims(sum(cleaned; dims=3); dims=3)
    end

    return PLMap(new_intensity, cleaned, m.x, m.y, m.pixel, m.metadata)
end
