# SpectroscopyTools.jl

[![CI](https://github.com/garrekstemo/SpectroscopyTools.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/garrekstemo/SpectroscopyTools.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://garrekstemo.github.io/SpectroscopyTools.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://garrekstemo.github.io/SpectroscopyTools.jl/dev/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/garrekstemo/SpectroscopyTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/garrekstemo/SpectroscopyTools.jl)

**Spectroscopy analysis for Julia — steady-state and ultrafast.**

For steady-state work, SpectroscopyTools provides peak fitting (Gaussian, Lorentzian, Voigt, Fano), baseline correction, peak detection, spectral transforms (Kramers-Kronig, Tauc, Kubelka-Munk), and unit conversions for FTIR, Raman, and UV-vis data. For ultrafast work, it provides typed data structures (`TATrace`, `TASpectrum`, `TAMatrix`), global fitting with decay-associated spectra, IRF-convolved exponential fitting, and chirp correction for broadband pump-probe experiments.

> **Note:** This package is currently being tested internally in our lab.

## Installation

```julia
using Pkg
Pkg.add("SpectroscopyTools")
```

## Quick Start

```julia
using SpectroscopyTools

# --- Ultrafast: kinetics with IRF-convolved biexponential ---
trace = TATrace(time, signal; wavelength=800.0)
fit = fit_exp_decay(trace; n_exp=2, irf=true)
report(fit)

# --- Ultrafast: global fit across a TA matrix ---
matrix = TAMatrix(time, wavelength, data)
gfit = fit_global(matrix; n_exp=2)   # shared τ, decay-associated spectra
spectra = das(gfit)                  # n_exp × n_wavelengths matrix

# --- Steady-state: peak fitting (FTIR, Raman, UV-vis) ---
result = fit_peaks(x, y, (2000, 2100))
report(result)

# --- Steady-state: baseline correction ---
bl = correct_baseline(x, y; method=:arpls)
y_corrected = bl.y

# --- Unit conversions (with Unitful) ---
wavelength_to_wavenumber(1500u"nm")
decay_time_to_linewidth(1.0u"ps")
```

## Features

| Module | Description |
|--------|-------------|
| **TA data types** | `TATrace`, `TASpectrum`, `TAMatrix` with semantic axis indexing (`m[λ=800]`, `m[t=1.0]`) |
| **Exponential decay** | Single/multi-exponential with optional IRF convolution |
| **Global fitting** | Shared parameters across traces, decay-associated spectra |
| **TA spectrum fitting** | N-peak model with ESA/GSB/SE labels and sign convention |
| **Chirp correction** | GVD detection (cross-correlation, threshold) and correction for broadband TA |
| **SVD filtering** | Matrix denoising for broadband TA data |
| **Peak fitting** | Gaussian, Lorentzian, Voigt, Pseudo-Voigt, Fano via [CurveFit.jl](https://github.com/garrekstemo/CurveFit.jl) |
| **Peak detection** | Automatic peak finding with prominence filtering |
| **Baseline correction** | arPLS, SNIP, rubber band, iModPoly, rolling ball |
| **Spectral math** | Smoothing, derivatives, band area, normalization, spectral arithmetic |
| **Transforms** | Kramers-Kronig, Kubelka-Munk, Tauc plot, SNV, Urbach tail |
| **Unit conversions** | Wavenumber, wavelength, energy, linewidth interconversion |
| **PL/Raman mapping** | Spatial maps, peak fitting, cosmic ray detection and removal |

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
