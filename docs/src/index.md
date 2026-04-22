# SpectroscopyTools.jl

A general-purpose spectroscopy analysis toolkit for Julia.

SpectroscopyTools provides peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, chirp correction, and unit conversions for spectroscopic data. It is designed to work with any spectroscopy discipline --- ultrafast, FTIR, Raman, UV-vis, fluorescence --- while keeping a minimal dependency footprint.

## Installation

```julia
using Pkg
Pkg.add("SpectroscopyTools")
```

## Quick Start

```julia
using SpectroscopyTools
using CurveFitModels  # for model functions like `lorentzian`
using Unitful         # for unit-aware conversions

# Simulate a noisy Lorentzian peak
x = collect(range(1900, 2200, length=600))
y = lorentzian([0.45, 2062.0, 24.0, 0.01], x) .+ 0.002 .* randn(length(x))

# Fit peaks in a region
result = fit_peaks(x, y, (2000, 2150))
report(result)

# Baseline correction
bl = correct_baseline(x, y; method=:arpls, λ=1e5)
y_corrected = bl.y

# Unit conversions
wavelength_to_wavenumber(1500u"nm")  # -> wavenumber in cm⁻¹
```

## Package Overview

| Module | What it does |
|--------|-------------|
| **Peak fitting** | Gaussian, Lorentzian, Pseudo-Voigt via CurveFit.jl |
| **Peak detection** | Automatic peak finding via Peaks.jl |
| **Baseline correction** | ALS, ARPLS, SNIP, polynomial, Whittaker |
| **Exponential decay** | Single/multi-exponential with optional IRF convolution |
| **Global fitting** | Shared parameters across multiple traces |
| **Chirp correction** | GVD detection and correction for broadband TA |
| **Unit conversions** | Wavenumber, wavelength, energy, linewidth interconversion |
| **Smoothing** | Savitzky-Golay filtering |

## Documentation Layout

This documentation follows the [Diátaxis](https://diataxis.fr/) framework:

- **Tutorials** — step-by-step walkthroughs for complete workflows (FTIR/Raman peak fitting, PL/Raman map analysis, chirp correction)
- **How-To Guides** — focused recipes for specific tasks (tuning detection sensitivity, fitting overlapping peaks, choosing a peak model, etc.)
- **Reference** — complete API documentation grouped by topic
- **Explanation** — background theory and design rationale for fit statistics and baseline algorithms

## Related Packages

- [CurveFitModels.jl](https://garrekstemo.github.io/CurveFitModels.jl/stable/) — lineshape and temporal model functions used for fitting
- [CurveFit.jl](https://github.com/garrekstemo/CurveFit.jl) — the nonlinear least-squares solver
- [QPSTools.jl](https://github.com/garrekstemo/QPSTools.jl) — lab-specific loaders, plotting themes, and eLabFTW integration built on top of SpectroscopyTools
