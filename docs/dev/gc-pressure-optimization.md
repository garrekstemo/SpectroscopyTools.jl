# GC Pressure in Per-Pixel Cosmic Ray Loops

Status: future optimization opportunity (not in current scope)

## Problem

`detect_cosmic_rays(m::PLMap)` and `remove_cosmic_rays(m::PLMap, result)` in `src/cosmic_rays.jl` iterate over every pixel in the spatial grid (e.g., 51x51 = 2,601 pixels). Each iteration allocates multiple temporary arrays. Under single-threaded execution, Julia's GC handles this without issue. Under `Threads.@threads` parallelism, all threads trigger GC simultaneously, causing stop-the-world contention that limits parallel speedup beyond 2-3x even with 6-8 threads available.

## Allocation Sites

### `detect_cosmic_rays(m::PLMap)` (two identical passes: noise estimation + flagging)

1. **`_most_similar_neighbor` centering arrays** (lines 254, 259):
   - `s_centered = signal .- s_mean` -- allocates a Vector{Float64}(n_ch) once per pixel
   - `nb_centered = nb .- nb_mean` -- allocates a Vector{Float64}(n_ch) per neighbor (up to 4x per pixel)
   - `s_centered .* nb_centered` (line 264) -- another temporary for the dot product broadcast

2. **Residual vector** (lines 367, 391):
   - `residual = Vector{Float64}(undef, n_ch)` -- allocated fresh each pixel, twice (noise pass + flagging pass)

3. **MAD intermediates** (lines 371, 398):
   - `residual .- median(residual)` -- temporary Vector{Float64}(n_ch)
   - `abs.(...)` on the above -- second temporary Vector{Float64}(n_ch)
   - Two arrays per MAD call, called twice per pixel (noise pass + flagging pass) = 8 temporary arrays per pixel

4. **`_neighbor_spectra` container** (lines 224, 361, 379):
   - Small vector of views, allocated per pixel. Minor compared to spectral-length arrays.

### `remove_cosmic_rays(m::PLMap, result)` (lines 463-529)

5. **`_most_similar_neighbor` again** (line 496): same centering allocations as above for each affected pixel.

### Total per-pixel allocation estimate (detect, both passes)

For n_ch = 512 channels, ~14 temporary Float64 vectors per pixel = ~57 KB per pixel. Over 2,601 pixels, that is ~144 MB of short-lived allocations triggering frequent GC.

## Potential Approaches

### 1. Thread-local workspace buffers

Pre-allocate a struct of scratch arrays per thread, indexed by `Threads.threadid()`. Reuse across pixels:

```julia
struct CRWorkspace
    residual::Vector{Float64}
    s_centered::Vector{Float64}
    nb_centered::Vector{Float64}
    abs_dev::Vector{Float64}  # for MAD calculation
end
```

Create `nthreads()` workspaces before the loop. Each thread grabs its own. Zero allocations in the hot loop.

### 2. In-place MAD calculation

Replace `median(abs.(residual .- med_r))` with a single-pass algorithm that reuses the residual buffer:

```julia
# Compute abs deviations in-place, then use partialsort! for median
for k in eachindex(buf)
    buf[k] = abs(residual[k] - med_r)
end
mad_val = partialsort!(buf, n_ch >> 1 + 1)  # destructive median via partial sort
```

Eliminates two intermediate arrays per MAD call (4 arrays saved per pixel in detect).

### 3. In-place Pearson correlation

Rewrite `_most_similar_neighbor` to accept workspace buffers and compute correlation without broadcast allocations:

```julia
function _most_similar_neighbor!(s_buf, nb_buf, signal, neighbors)
    s_buf .= signal .- mean(signal)
    s_ss = sqrt(sum(abs2, s_buf))
    for (i, nb) in enumerate(neighbors)
        nb_buf .= nb .- mean(nb)
        corr = dot(s_buf, nb_buf) / (s_ss * sqrt(sum(abs2, nb_buf)))
        ...
    end
end
```

Uses `dot()` from LinearAlgebra instead of `sum(a .* b)` to avoid the elementwise temporary.

### 4. Merge the two passes in detect

The noise-estimation pass (lines 360-375) and flagging pass (lines 378-421) both compute the MSN and residual for every pixel. Merging into a single pass with a deferred flagging step would halve the number of `_most_similar_neighbor` calls and residual allocations.

## Expected Benefit

- Single-threaded: modest improvement (less GC pressure, ~10-20% faster).
- Multi-threaded: significant improvement. GC contention is currently the bottleneck limiting scaling beyond ~2-3x. With zero-allocation inner loops, expect near-linear scaling up to available cores (target: 4-6x on 8-thread machines).
- Combined with `Threads.@threads` parallelization, estimated total speedup: 5-10x over current single-threaded baseline for a 51x51 map.
