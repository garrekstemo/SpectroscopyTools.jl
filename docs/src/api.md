# [API Reference](@id api-reference)

## Types

```@docs
SpectroscopyTools.SpectroscopyTools
AbstractSpectroscopyData
TATrace
TASpectrum
TAMatrix
TASpectrumFit
TAPeak
```

## [Chirp Correction](@id chirp-api)

```@docs
ChirpCalibration
detect_chirp
correct_chirp
subtract_background
polynomial
save_chirp
load_chirp
```

## SVD Filtering

```@docs
svd_filter
singular_values
```

## Cosmic Ray Detection

```@docs
CosmicRayResult
CosmicRayMapResult
detect_cosmic_rays
remove_cosmic_rays
```

## PL / Raman Mapping

```@docs
PLMap
extract_spectrum
integrated_intensity
intensity_mask
peak_centers
normalize_intensity
```

## Peak Fitting

```@docs
fit_peaks
fit_ta_spectrum
MultiPeakFitResult
PeakFitResult
anharmonicity
```

## Peak Detection

```@docs
PeakInfo
find_peaks
peak_table
```

## Exponential Decay Fitting

```@docs
fit_exp_decay
fit_global
ExpDecayFit
MultiexpDecayFit
GlobalFitResult
das
report
```

## Baseline Correction

```@docs
als_baseline
arpls_baseline
snip_baseline
correct_baseline
```

## Spectroscopy Utilities

```@docs
normalize
smooth_data
calc_fwhm
subtract_spectrum
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
