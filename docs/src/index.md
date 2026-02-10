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

# Fit peaks in a spectrum
result = fit_peaks(wavenumber, absorbance, (2000, 2100))
report(result)

# Baseline correction
baseline = als_baseline(y; lambda=1e5, p=0.01)
y_corrected = y .- baseline

# Unit conversions
wavelength_to_wavenumber(1500u"nm")  # -> wavenumber in cm^-1
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

## Tutorials

```@contents
Pages = ["tutorials/chirp_correction.md"]
Depth = 2
```

## API Reference

See the full [API Reference](@ref api-reference) for documentation of all exported functions and types.
