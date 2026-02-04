"""
SpectroscopyTools.jl — General-purpose spectroscopy analysis toolkit.

Provides data types, fitting routines, baseline correction, peak detection,
unit conversions, and utility functions for spectroscopic data analysis.

Extracted from QPS.jl to serve as a standalone, reusable foundation.
"""
module SpectroscopyTools

# Dependencies
using CurveFit
using CurveFitModels
using Statistics
using LinearAlgebra
using SparseArrays
using Unitful
import PhysicalConstants.CODATA2022: h, c_0, ħ

# Source files (order matters: types before functions that use them)
include("types.jl")
include("units.jl")
include("spectroscopy.jl")
include("baseline.jl")
include("peakdetection.jl")
include("peakfitting.jl")
include("fitting.jl")

# ==========================================================================
# Exports — Types
# ==========================================================================
export AbstractSpectroscopyData
export xdata, ydata, zdata, xlabel, ylabel, zlabel, is_matrix
export source_file, npoints, title
export AxisType, time_axis, wavelength_axis
export PumpProbeData, xaxis, xaxis_label
export TATrace, TASpectrum, TAMatrix
export TASpectrumFit, GlobalFitResult
export ExpDecayFit, ExpDecayIRFFit, BiexpDecayFit, MultiexpDecayFit
export n_exp, weights
export PeakFitResult, MultiPeakFitResult
export PumpProbeResult
export report, format_results

# ==========================================================================
# Exports — Fitting
# ==========================================================================
export fit_exp_decay, fit_biexp_decay, fit_global, fit_ta_spectrum
export fit_decay, fit_decay_irf
export irf_fwhm, pulse_fwhm

# ==========================================================================
# Exports — Peak fitting
# ==========================================================================
export fit_peaks, predict_peak, predict_baseline

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
export normalize, time_index
export calc_ΔA, fit_decay_trace, extract_tau, fit_global_decay
export subtract_spectrum, linear_baseline_correction
export smooth_data, savitzky_golay, calc_fwhm
export transmittance_to_absorbance, absorbance_to_transmittance

# ==========================================================================
# Exports — Unit conversions
# ==========================================================================
export parse_concentration, parse_time
export wavenumber_to_wavelength, wavelength_to_wavenumber
export wavenumber_to_energy, wavelength_to_energy, energy_to_wavelength
export linewidth_to_decay_time, decay_time_to_linewidth

# ==========================================================================
# Re-exports from CurveFit
# ==========================================================================
export solve, NonlinearCurveFitProblem
export coef, stderror, confint
export residuals, rss, mse, nobs, predict, isconverged

# ==========================================================================
# Re-exports from CurveFitModels
# ==========================================================================
export gaussian, lorentzian, single_exponential, pseudo_voigt
export n_exponentials

end # module SpectroscopyTools
