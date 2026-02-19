# SpectroscopyTools.jl

## Project Mission

**SpectroscopyTools.jl provides general-purpose spectroscopy analysis tools for the Julia ecosystem.** This is the public, registered package extracted from QPS.jl — containing peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, and unit conversions.

The package serves anyone doing spectroscopy in Julia: ultrafast, FTIR, Raman, UV-vis, fluorescence. Lab-specific features (instrument I/O, sample registries, eLabFTW integration) stay in QPS.jl.

---

## Package Scope

### In Scope

- Peak fitting (Gaussian, Lorentzian, Pseudo-Voigt) via CurveFit.jl + CurveFitModels.jl
- Peak detection via Peaks.jl
- Baseline correction (ALS, ARPLS, SNIP, polynomial, Whittaker)
- Exponential decay fitting (single/multi-exponential) with optional IRF convolution
- Global fitting with shared parameters across traces
- Unit conversions (wavenumber/wavelength/energy, linewidth/decay-time)
- Normalization, smoothing, spectral math
- Chirp correction for broadband TA (detection, correction, serialization)
- Background subtraction for TA matrices
- Typed spectroscopy data (`AbstractSpectroscopyData` hierarchy)
- Plotting via Makie extension (no themes — users set their own)

### Out of Scope (stays in QPS.jl)

- Instrument-specific I/O (LabVIEW .lvm, JASCO CSV)
- Sample registries and metadata lookup
- eLabFTW integration
- Lab themes (`qps_theme`, `publication_theme`)
- Wavenumber calibration constants

---

## Package Dependencies

**Target: < 15 direct dependencies.** Keep the dependency tree minimal.

**Compat entries required**: Since SpectroscopyTools.jl will be a public registered package, all non-stdlib dependencies MUST have `[compat]` entries in Project.toml. When adding a new dependency with `Pkg.add()`, verify that a compat entry was added automatically. This differs from QPSTools.jl, which uses loose constraints as a private lab package.

### Core Dependencies

| Package | Purpose |
|---------|---------|
| CurveFit.jl | Fitting backend |
| CurveFitModels.jl | Peak/decay model functions |
| Unitful.jl | Physical units |
| PhysicalConstants.jl | CODATA constants |
| SavitzkyGolay.jl | Smoothing |
| Peaks.jl | Peak detection |
| Interpolations.jl | Cubic spline interpolation (chirp correction) |
| JSON.jl | Chirp calibration serialization |
| SparseArrays (stdlib) | Baseline correction |
| LinearAlgebra (stdlib) | Fitting internals |
| Statistics (stdlib) | Basic stats |

### CurveFit.jl API

CurveFit.jl is the nonlinear fitting backend. Problem-solver pattern:

```julia
using CurveFit

# Model: fn(params, x) — parameters first
fn(p, x) = @. p[1] + p[2] * exp(-x / p[3])

prob = NonlinearCurveFitProblem(fn, p0, x, y)
sol = solve(prob)

# Solution stats (StatsAPI interface)
coef(sol)          # Fitted parameters
residuals(sol)     # y_data - y_fit (observed - predicted)
fitted(sol)        # Fitted values at original x points
predict(sol, x)    # Evaluate at new x points
rss(sol)           # Residual sum of squares
mse(sol)           # Mean squared error (rss / dof_residual)
stderror(sol)      # Standard errors of coefficients
confint(sol)       # 95% confidence intervals (vector of tuples)
nobs(sol)          # Number of observations
dof(sol)           # Model degrees of freedom
dof_residual(sol)  # Residual DOF (nobs - dof)
vcov(sol)          # Variance-covariance matrix
isconverged(sol)   # Convergence status
```

**Important:** SpectroscopyTools extends `CurveFit.residuals`, `CurveFit.predict`, and `CurveFit.fitted` with methods for its own fit result types. Use `import CurveFit: residuals, predict, fitted` (not just `using`) to enable method extension. CurveFit does NOT provide R² — compute as `1 - rss(sol) / ss_tot`.

### Extensions (Weak Dependencies)

| Extension | Trigger | Purpose |
|-----------|---------|---------|
| `SpectroscopyToolsMakieExt` | CairoMakie or GLMakie | Plotting functions |

Plotting is an extension, not a hard dependency. Users who only need fitting don't pay the Makie load time.

### NOT Dependencies (stay in QPS.jl)

- HTTP (eLabFTW only)
- JASCOFiles (instrument-specific)
- FileIO, Dates, DelimitedFiles (LVM parsing only)

---

## Type Hierarchy

```
AbstractSpectroscopyData (root interface)
├── TATrace          (kinetics: signal vs time at fixed wavelength)
├── TASpectrum       (spectrum: signal vs wavenumber at fixed time)
└── TAMatrix         (2D: time × wavelength heatmap)
```

Note: `AnnotatedSpectrum` (with metadata, FTIR/Raman subtypes) is defined in QPS.jl, not here.

### Uniform Interface

All types implement:

```julia
xdata(data)       # Primary x-axis
ydata(data)       # Signal (1D) or secondary axis (2D)
zdata(data)       # Matrix data (2D only)
xlabel(data)      # X-axis label
ylabel(data)      # Y-axis label
is_matrix(data)   # true for 2D data
source_file(data) # Source filename
npoints(data)     # Number of points
title(data)       # Display title
```

---

## API Design Principles

### Dual Interface: Typed Data or Raw Vectors

Functions accept either typed spectroscopy data (preferred — more information, better defaults) or raw vectors (always available for users with their own data structures).

```julia
# Typed path (preferred)
result = fit_exp_decay(trace)                    # TATrace input
result = fit_peaks(spectrum, (2000, 2100))       # AnnotatedSpectrum input

# Raw vector path
result = fit_exp_decay(time, signal; irf=true)   # Vector input
baseline = als_baseline(y; lambda=1e5, p=0.01)   # Vector input
```

### Model Functions from CurveFitModels.jl

Never define fitting functions inline. Use the standard models from CurveFitModels.jl:

- `lorentzian` — `p = [A, x₀, Γ, offset]`
- `gaussian` — `p = [A, x₀, Γ, offset]`
- `pseudo_voigt` — `p = [A, x₀, σ, η]`
- `single_exponential` — `p = [A, τ, offset]`

All models use signature `fn(p, x)` (parameters first) and are ForwardDiff-compatible.

### Fit Results are Structs, Not Tuples

Every fitting function returns a proper result type with accessors:

```julia
result = fit_exp_decay(trace)
result.tau        # Time constant
result.sigma      # IRF width
result.amplitude  # Amplitude
result.rsquared   # Fit quality

predict(result)   # Fitted curve
residuals(result) # Residuals
report(result)    # Formatted output
```

### No Themes in Plotting

The plotting extension provides functions (`plot_spectrum`, `plot_kinetics`, `plot_peaks`), not aesthetics. Users set their own Makie themes. Only use inline styling for semantic distinction (e.g., `color=:red` for fit vs data).

---

## Development Standards

### Julia Requirements

- **Minimum version**: Julia 1.10 (LTS)
- **Development workflow**: `julia --project=.` and `using Revise`

### Code Style

```julia
# Iteration
for i in eachindex(arr)          # ✓ Not 1:length(arr)

# Nil coalescing
something(a, default)            # ✓ Not ternary

# Array allocation (ForwardDiff compatible)
similar(p, N)                    # ✓ Not zeros(N)

# Number formatting
string(round(x, digits=2))       # ✓ Not @sprintf
```

### Type Flexibility

Use generic signatures, not `Float64`-locked:

```julia
# Good
function fit_peaks(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

# Bad (locks out Float32, Measurement, dual numbers)
function fit_peaks(x::Vector{Float64}, y::Vector{Float64})
```

### Docstrings

Every exported function needs:

```julia
"""
    fit_peaks(spectrum, region; kwargs...)

Fit peaks in the specified region of a spectrum.

# Arguments
- `spectrum`: Input spectrum (AbstractSpectroscopyData or x, y vectors)
- `region`: Tuple `(x_min, x_max)` defining the fitting region

# Keywords
- `model=lorentzian`: Peak shape function
- `n_peaks=:auto`: Number of peaks (`:auto` uses peak detection)

# Returns
`MultiPeakFitResult` with fitted parameters and uncertainties.

# Examples
```julia
result = fit_peaks(spec, (2000, 2100))
report(result)
```
"""
```

### Tests

- All tests use synthetic data — no local file dependencies
- Every exported function has at least one test
- Target > 80% code coverage
- Test edge cases: empty data, single point, NaN handling

### Issue Workflow

When working through bugs, cleanup, or design issues, address each one individually:

1. **Explain the issue** — what's wrong, with a concrete code example showing the failure
2. **Propose a fix** — what changes, why this approach
3. **Connect to the larger goal** — how does this fix align with the package's design principles (type flexibility, ForwardDiff compat, dual interface, etc.)
4. **Implement and test** — make the change, run the test suite

Do not batch unrelated fixes. One issue at a time keeps changes reviewable and reversible.

---

## Package Structure

```
SpectroscopyTools.jl/
├── Project.toml
├── LICENSE (MIT)
├── README.md
├── src/
│   ├── SpectroscopyTools.jl    # Main module
│   ├── types.jl                # AbstractSpectroscopyData + fit result types
│   ├── fitting.jl              # Exponential decay (single/multi/global + IRF)
│   ├── peakfitting.jl          # Multi-peak fitting + TA spectrum fitting
│   ├── peakdetection.jl        # Peak finding
│   ├── baseline.jl             # ALS, ARPLS, SNIP, polynomial, Whittaker
│   ├── spectroscopy.jl         # Normalize, conversions, smoothing
│   ├── units.jl                # Unitful conversions
│   └── chirp.jl                # Chirp detection, correction, serialization
├── ext/
│   └── SpectroscopyToolsMakieExt.jl  # (not yet implemented)
├── test/
│   └── runtests.jl             # All tests in one file (345 tests)
└── docs/
    ├── Project.toml
    ├── make.jl
    └── src/
        ├── index.md
        ├── tutorials/
        └── api.md
```

---

## Exports

```julia
# Types & Interface
export AbstractSpectroscopyData
export xdata, ydata, zdata, xlabel, ylabel, zlabel, is_matrix
export source_file, npoints, title
export TATrace, TASpectrum, TAMatrix

# Fit Results
export ExpDecayFit, MultiexpDecayFit, GlobalFitResult
export MultiPeakFitResult, PeakFitResult
export TAPeak, TASpectrumFit, anharmonicity, das

# Fitting
export fit_exp_decay, fit_global
export fit_peaks, predict_peak, predict_baseline
export fit_ta_spectrum
export report, polynomial

# Chirp Correction
export ChirpCalibration, detect_chirp, correct_chirp, subtract_background
export save_chirp, load_chirp
export svd_filter, singular_values

# Peak Detection
export PeakInfo, find_peaks, peak_table

# Baseline
export als_baseline, arpls_baseline, snip_baseline
export correct_baseline

# Spectroscopy Utilities
export normalize, subtract_spectrum
export smooth_data, calc_fwhm
export transmittance_to_absorbance, absorbance_to_transmittance

# Units
export wavenumber_to_wavelength, wavelength_to_wavenumber
export wavenumber_to_energy, wavelength_to_energy, energy_to_wavelength
export linewidth_to_decay_time, decay_time_to_linewidth

# Re-exports from CurveFit.jl
export solve, NonlinearCurveFitProblem
export coef, residuals, predict, fitted, stderror, confint, rss, mse, nobs, isconverged

# Re-exports from CurveFitModels.jl
export gaussian, lorentzian, pseudo_voigt, single_exponential
```

---

## Three-Package Boundary Design

The stack has three layers with clear separation of concerns:

| Package | Role | Owns |
|---------|------|------|
| **CurveFitModels.jl** | Model functions only | `gaussian`, `lorentzian`, `pseudo_voigt`, `single_exponential`, `n_exponentials`, `poly`, `combine`, etc. All `fn(p, x)` signature. Zero dependencies. |
| **SpectroscopyTools.jl** | Analysis engine | Types, fitting (`fit_exp_decay`, `fit_peaks`, `fit_global`), baseline correction, peak detection, chirp correction (`detect_chirp`, `correct_chirp`, `subtract_background`), unit conversions, smoothing. Uses CurveFitModels for model functions. |
| **QPS.jl** | Lab-specific layer | Instrument I/O (`load_lvm`, `load_ftir`, `load_raman`), sample registry, eLabFTW integration, Makie plotting (themes, layers, layouts, `plot_chirp`). Re-exports SpectroscopyTools. |

### Boundary Rules

- **CurveFitModels** defines model functions. It has no analysis logic, no types, no I/O.
- **SpectroscopyTools** calls CurveFitModels functions internally (e.g., `n_exponentials(n)` inside `fit_exp_decay`). It does NOT re-export model constructors that are only used internally.
- **QPS.jl** `import`s SpectroscopyTools functions to extend them with `AnnotatedSpectrum` dispatches. It re-exports user-facing SpectroscopyTools names.
- Lab members use `using QPS` — SpectroscopyTools.jl is invisible to them.

### Functions QPS Extends via `import`

These must remain exported from SpectroscopyTools:

```julia
find_peaks, fit_peaks
transmittance_to_absorbance, absorbance_to_transmittance
subtract_spectrum, correct_baseline
xdata, ydata, xlabel, ylabel, source_file
normalize
```

### Functions QPS Re-exports via `import` (no extension)

These are imported so QPSTools can re-export them:

```julia
n_exp, weights, anharmonicity, format_results
ChirpCalibration, polynomial, detect_chirp, correct_chirp
subtract_background, save_chirp, load_chirp
```

---

## Pre-Registration Checklist

- [x] Package loads: `using SpectroscopyTools`
- [x] All tests pass: `Pkg.test("SpectroscopyTools")` (391 tests)
- [x] Docs structure (make.jl, index.md, api.md, tutorials)
- [x] Docs build and deploy (CI)
- [x] CI workflow exists (.github/workflows/CI.yml)
- [x] CI green on Linux, macOS, Windows (Julia 1.10, 1.12, nightly)
- [x] README has badges, installation, quickstart
- [x] LICENSE file present
- [x] No local path dependencies in Project.toml
- [x] All direct dependencies registered (CurveFitModels.jl registered)
- [x] < 15 direct dependencies (11 total, 8 non-stdlib)

---

## References

- `TODO.md` — Active development tasks (chirp R² improvement)
