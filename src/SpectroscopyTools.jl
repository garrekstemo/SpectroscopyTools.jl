"""
SpectroscopyTools.jl — General-purpose spectroscopy analysis toolkit.

Provides data types, fitting routines, baseline correction, peak detection,
unit conversions, and utility functions for spectroscopic data analysis.

Extracted from QPS.jl to serve as a standalone, reusable foundation.
"""
module SpectroscopyTools

# Dependencies
using CurveFit
import CurveFit: residuals, predict, fitted
using CurveFitModels
using Statistics
using LinearAlgebra
using SparseArrays
using Unitful
using Peaks: findmaxima, peakproms!, peakwidths!
using SavitzkyGolay: savitzky_golay as _sg_filter
using Interpolations
using JSON
import PhysicalConstants.CODATA2022: h, c_0, ħ

# Source files (order matters: types before functions that use them)
include("types.jl")
include("units.jl")
include("spectroscopy.jl")
include("baseline.jl")
include("peakdetection.jl")
include("peakfitting.jl")
include("fitting.jl")
include("chirp.jl")

# ==========================================================================
# Exports — Types & Interface
# ==========================================================================
export AbstractSpectroscopyData
export xdata, ydata, zdata, xlabel, ylabel, zlabel, is_matrix
export source_file, npoints, title
export TATrace, TASpectrum, TAMatrix

# ==========================================================================
# Exports — Fit Results
# ==========================================================================
export ExpDecayFit, MultiexpDecayFit, GlobalFitResult
export MultiPeakFitResult, PeakFitResult
export TAPeak, TASpectrumFit, anharmonicity, das

# ==========================================================================
# Exports — Fitting
# ==========================================================================
export fit_exp_decay, fit_global
export fit_peaks, predict_peak, predict_baseline
export fit_ta_spectrum
export report, polynomial

# ==========================================================================
# Exports — Chirp correction
# ==========================================================================
export ChirpCalibration, detect_chirp, correct_chirp, subtract_background
export save_chirp, load_chirp
export svd_filter, singular_values

# ==========================================================================
# Exports — Peak detection
# ==========================================================================
export PeakInfo, find_peaks, peak_table

# ==========================================================================
# Exports — Baseline correction
# ==========================================================================
export als_baseline, arpls_baseline, snip_baseline
export correct_baseline

# ==========================================================================
# Exports — Spectroscopy utilities
# ==========================================================================
export normalize, subtract_spectrum
export smooth_data, calc_fwhm
export transmittance_to_absorbance, absorbance_to_transmittance

# ==========================================================================
# Exports — Unit conversions
# ==========================================================================
export wavenumber_to_wavelength, wavelength_to_wavenumber
export wavenumber_to_energy, wavelength_to_energy, energy_to_wavelength
export linewidth_to_decay_time, decay_time_to_linewidth

# ==========================================================================
# Re-exports from CurveFit.jl
# ==========================================================================
export solve, NonlinearCurveFitProblem
export coef, residuals, predict, fitted, stderror, confint, rss, mse, nobs, isconverged

# ==========================================================================
# Re-exports from CurveFitModels.jl
# ==========================================================================
export gaussian, lorentzian, pseudo_voigt, single_exponential

end # module SpectroscopyTools
