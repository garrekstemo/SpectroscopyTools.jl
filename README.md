# SpectroscopyTools.jl

[![CI](https://github.com/garrekstemo/SpectroscopyTools.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/garrekstemo/SpectroscopyTools.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://garrekstemo.github.io/SpectroscopyTools.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://garrekstemo.github.io/SpectroscopyTools.jl/dev/)

A general-purpose spectroscopy analysis toolkit for Julia.

SpectroscopyTools provides peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, chirp correction, and unit conversions for spectroscopic data. It works with any spectroscopy discipline: ultrafast, FTIR, Raman, UV-vis, fluorescence.

> **Note:** This package is currently being tested internally in our lab.

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

# Exponential decay fitting
trace = TATrace(time, signal; wavelength=800.0)
fit = fit_exp_decay(trace; n_exp=2, irf=true)
report(fit)

# Baseline correction
baseline = als_baseline(y; lambda=1e5, p=0.01)
y_corrected = y .- baseline

# Unit conversions (with Unitful)
wavelength_to_wavenumber(1500u"nm")
decay_time_to_linewidth(1.0u"ps")
```

## Features

| Module | Description |
|--------|-------------|
| **Peak fitting** | Gaussian, Lorentzian, Pseudo-Voigt via [CurveFit.jl](https://github.com/garrekstemo/CurveFit.jl) |
| **TA spectrum fitting** | N-peak model with ESA/GSB/SE labels and sign convention |
| **Peak detection** | Automatic peak finding with prominence filtering |
| **Baseline correction** | ALS, ARPLS, SNIP |
| **Exponential decay** | Single/multi-exponential with optional IRF convolution |
| **Global fitting** | Shared parameters across multiple traces |
| **Chirp correction** | GVD detection and correction for broadband TA |
| **SVD filtering** | Matrix denoising for broadband TA data |
| **Unit conversions** | Wavenumber, wavelength, energy, linewidth interconversion |
| **Smoothing** | Savitzky-Golay filtering |

## Data Types

SpectroscopyTools provides typed containers for spectroscopy data:

```julia
# Kinetics at a single wavelength
trace = TATrace(time, signal; wavelength=800.0)

# Spectrum at a fixed time delay
spectrum = TASpectrum(wavenumber, signal; time_delay=1.0)

# 2D time x wavelength matrix with semantic indexing
matrix = TAMatrix(time, wavelength, data)
matrix[λ=800]   # -> TATrace at nearest wavelength
matrix[t=1.0]   # -> TASpectrum at nearest time delay
```

## Dependencies

SpectroscopyTools uses [CurveFit.jl](https://github.com/garrekstemo/CurveFit.jl) as the nonlinear fitting backend and [CurveFitModels.jl](https://github.com/garrekstemo/CurveFitModels.jl) for model functions. Unit conversions use [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

## License

MIT
