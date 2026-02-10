# SpectroscopyTools.jl — TODO

## 1. Improve chirp detection R² on real data

**Status**: Working but inaccurate on noisy CCD data (R² = 0.492 on PB-film)

**Context**: `detect_chirp` works well on synthetic data (R² > 0.98) but struggles with real broadband TA from a CCD camera (2048 wavelengths, 151 time points, ~467–733 nm). The CCD data is noisy and the chirp onset is subtle relative to later dynamics.

**Current implementation** (`src/chirp.jl`):
- Two methods: `:xcorr` (cross-correlation of absolute onset gradients) and `:threshold` (half-maximum crossing)
- Wavelength binning (`bin_width`) averages adjacent columns before detection
- Savitzky-Golay smoothing before gradient/threshold calculation
- MAD-based outlier rejection after polynomial fitting
- `_auto_chirp_range()` restricts detection to the onset region (10% threshold → peak gradient)

**What to try**:
- Better wavelength binning: adaptive bin widths based on signal strength, or allow user to pass pre-binned data
- Weighted polynomial fitting: weight by signal strength or inverse noise estimate per bin
- Two-stage detection: coarse detection with wide bins first, then refine with narrow bins
- Alternative onset detection: fit a step function or error function to each bin's time trace rather than using raw gradient/threshold
- Pre-detection denoising: apply SVD filter (see task 2) before chirp detection

**Test data**: Real CCD data lives in `QPSTools.jl/data/broadbandTA_experiment/`. The example script is `QPSTools.jl/examples/broadband_ta_example.jl`. To test against real data, you'll need both repos checked out. Synthetic tests in `test/runtests.jl` cover correctness; real-data R² is the quality metric.

**Key gotcha**: The chirp polynomial gives `t_shift(λ)` = time the pump arrives at each wavelength relative to reference. To correct: evaluate `original(t + t_shift)`, NOT `original(t - t_shift)`.

---

## 2. Implement SVD filter for TA matrix denoising

**Status**: Not started

**Context**: SVD (singular value decomposition) filtering is a standard preprocessing step for broadband TA data. The time × wavelength matrix is decomposed into singular components, and only the first N components (containing real signal) are kept, discarding noise-dominated components.

**Where it goes**: `src/chirp.jl` (alongside `subtract_background`) or a new `src/preprocessing.jl` — either works, since it's a TA matrix preprocessing step.

**API design**:
```julia
# Filter matrix keeping only first n_components singular values
filtered = svd_filter(matrix::TAMatrix; n_components::Int=5) -> TAMatrix

# Diagnostic: inspect singular values to choose n_components
singular_values(matrix::TAMatrix) -> Vector{Float64}
```

**Implementation sketch**:
1. `U, S, V = svd(matrix.data)` (LinearAlgebra is already a dependency)
2. Zero out `S[n+1:end]`
3. Reconstruct: `U * Diagonal(S) * V'`
4. Return new `TAMatrix` with filtered data and metadata noting the filter

**Considerations**:
- Should go before chirp detection in the preprocessing pipeline (denoise → background subtract → detect chirp → correct chirp)
- The `n_components` choice matters: too few removes signal, too many keeps noise. A scree plot / elbow detection helper would be useful
- Plotting the singular value spectrum is a QPSTools concern (Makie-dependent), not SpectroscopyTools — but returning the values for the user to inspect is in scope

**Tests**: Synthetic matrix with known rank + noise. After filtering with correct n_components, the reconstructed matrix should be close to the noise-free original.

---

## 3. Multi-exponential global analysis for TAMatrix

**Status**: Not started (extends existing `fit_global`)

**Context**: Global analysis fits all wavelengths of a broadband TA matrix simultaneously with shared time constants. The per-wavelength amplitudes for each decay component form the **decay-associated spectra (DAS)**, which reveal which spectral features share dynamics. This is the standard analysis for broadband TA data.

**What already exists** (`src/fitting.jl`):
- `fit_global(traces::Vector{TATrace})` — single shared τ, per-trace amplitude + offset
- `GlobalFitResult` struct — single τ, single σ, per-trace amplitudes/offsets
- `predict`, `report`, `format_results` for `GlobalFitResult`
- `_exp_decay_irf_conv` and `_multiexp_irf_conv` helper functions
- `fit_exp_decay` already supports `n_exp` for single traces

**What needs to be built**:

### 3a. Extend `fit_global` to support multi-exponential

The current `fit_global` only fits a single shared τ. Extend to N shared time constants:

```julia
# Current (keep working):
result = fit_global(traces; irf_width=0.15)

# New — multi-exponential:
result = fit_global(traces; n_exp=2, irf_width=0.15)
```

Shared parameters: `τ₁, τ₂, ..., τₙ, σ, t₀`. Per-trace parameters: `A₁, A₂, ..., Aₙ, offset`.

The `GlobalFitResult` struct needs to change:
- `tau::Float64` → `taus::Vector{Float64}` (or keep `tau` for n_exp=1 backward compat)
- `amplitudes::Vector{Float64}` → `amplitudes::Matrix{Float64}` (n_traces × n_exp)

Target struct (replaces existing `GlobalFitResult`):

```julia
struct GlobalFitResult
    taus::Vector{Float64}              # Shared time constants [τ₁, τ₂, ...]
    taus_err::Vector{Float64}          # Uncertainties on τ values
    sigma::Float64                     # IRF width
    t0::Float64                        # Time zero
    amplitudes::Matrix{Float64}        # n_traces × n_exp (per-trace, per-component)
    offsets::Vector{Float64}           # Per-trace offsets
    labels::Vector{String}             # Trace labels
    das::Union{Nothing, Matrix{Float64}}  # n_exp × n_wavelengths (nothing for Vector{TATrace} input)
    wavelength::Union{Nothing, Vector{Float64}}  # Wavelength axis for DAS
    rsquared::Float64
    rsquared_individual::Vector{Float64}
    residuals::Vector{Vector{Float64}}
end
```

**Important**: The existing single-τ `fit_global` must still work — the `n_exp=1` case should produce equivalent results. Existing tests must pass.

### 3b. Add `TAMatrix` dispatch

```julia
# Fit all wavelengths (extracts traces internally)
result = fit_global(matrix::TAMatrix; n_exp=2, λ=nothing)

# Fit selected wavelengths only (faster)
result = fit_global(matrix::TAMatrix; n_exp=2, λ=[450, 500, 550, 600])
```

When `λ` is not specified, extract traces at every wavelength (or a reasonable subset — every Nth for very dense CCD data). When `λ` is a vector, extract traces at those wavelengths using `matrix[λ=val]`.

### 3c. DAS extraction

The decay-associated spectra are the per-wavelength amplitudes for each time constant:

```julia
# DAS is amplitudes as a function of wavelength for each τ
result.das        # Matrix: n_exp × n_wavelengths
result.wavelength # Vector of wavelengths used in the fit
result.taus       # Vector of shared time constants
```

DAS(λ, i) = amplitude of the i-th exponential component at wavelength λ. This falls out naturally from the fit — the per-wavelength amplitudes ARE the DAS.

### 3d. Update `predict`, `report`, `format_results`

- `predict(result, matrix)` → reconstructed TAMatrix
- `predict(result, traces)` → vector of fitted curves (existing, update for multi-exp)
- `report(result)` should show all τ values and a summary of DAS features
- `format_results(result)` markdown table with τ values

**Plotting note**: `plot_das()` is a QPSTools concern (Makie-dependent). SpectroscopyTools just provides the data via `result.das` and `result.wavelength`.

**Error estimation**: `taus_err` comes from `stderror(sol)` via CurveFit.jl (covariance-based). Bootstrap estimation (`errors=:bootstrap`) is a future enhancement, not needed for v1.

**Future extension — target analysis**: After global analysis works, the next level is `fit_target(matrix; model=:sequential)` which fits a kinetic model (A → B → C) to extract species-associated spectra (SAS) instead of DAS. This is out of scope for now but the `GlobalFitResult` design should not preclude it.

**Tests**:
- Synthetic TAMatrix with known 2-component dynamics → recovered τ values match input
- Single-component case (`n_exp=1`) produces same results as existing `fit_global`
- DAS shape matches expected spectral features
- Edge cases: single wavelength, very noisy data

---

## General notes for the bot

- All tests use synthetic data (no local file dependencies). See existing tests in `test/runtests.jl` for patterns.
- `_sg_filter` is the module-level alias for `savitzky_golay` from SavitzkyGolay.jl — use it directly, no `using` needed.
- `Interpolations` and `JSON` are available at module level.
- `LinearAlgebra` is available at module level (for SVD).
- `Statistics: mean, median, std` are available at module level.
- Run tests with: `cd SpectroscopyTools.jl && julia --project=. -e 'using Pkg; Pkg.test()'`
- This is a public package — add `[compat]` entries for any new dependencies.
- Follow existing docstring style (see `detect_chirp`, `correct_chirp` for examples).
