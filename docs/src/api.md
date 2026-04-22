# [Additional API Reference](@id api-reference)

This page collects exported names that aren't covered by the grouped reference pages (peak detection, peak fitting, baseline correction, preprocessing, PL/Raman mapping).

## Module

```@docs
SpectroscopyTools.SpectroscopyTools
```

## Transient Absorption Types

```@docs
AbstractSpectroscopyData
TATrace
TASpectrum
TAMatrix
TASpectrumFit
TAPeak
fit_ta_spectrum
anharmonicity
```

## Chirp Correction

```@docs
ChirpCalibration
detect_chirp
correct_chirp
polynomial
save_chirp
load_chirp
```

## SVD Filtering

```@docs
svd_filter
singular_values
```

## Exponential Decay Fitting

```@docs
fit_exp_decay
fit_global
ExpDecayFit
MultiexpDecayFit
GlobalFitResult
das
```

## Spectroscopy Utilities

```@docs
normalize
calc_fwhm
transmittance_to_absorbance
absorbance_to_transmittance
npoints
source_file
title
```

## Unit Conversions

```@docs
wavenumber_to_wavelength
wavelength_to_wavenumber
wavenumber_to_energy
wavelength_to_energy
energy_to_wavelength
linewidth_to_decay_time
decay_time_to_linewidth
```
