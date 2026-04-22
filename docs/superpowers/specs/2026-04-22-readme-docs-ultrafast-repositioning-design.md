# README and Docs Repositioning — Ultrafast Emphasis

**Date:** 2026-04-22
**Status:** Design
**Scope:** Cosmetic only (no new tutorials, no new reference pages)

## Motivation

SpectroscopyTools.jl currently presents itself as "a general-purpose spectroscopy analysis toolkit," which puts it in direct competition with the more mature, registered Spectra.jl on 1D preprocessing (baseline, smoothing, peak detection).
The package's real differentiation is in ultrafast / time-resolved spectroscopy: typed TA data structures (`TATrace`, `TASpectrum`, `TAMatrix`), IRF-convolved exponential fitting, global fitting with decay-associated spectra, and chirp correction for broadband pump-probe data — none of which exist in Spectra.jl.

New users landing on the README should immediately see the ultrafast capabilities. Steady-state spectroscopy support remains first-class but should be framed as a complementary capability, not the headline.

## Non-goals

- No rename of the package.
- No new tutorials in this pass (real TA data not yet available for a global-analysis walkthrough).
- No new reference pages in this pass.
- No changes to source code, public API, or exports.
- No removal of existing steady-state capabilities.

## Design

### 1. Headline and opening paragraph

Replace `README.md:9–11` and the parallel lines in `docs/src/index.md:3–5`.

**Before:**
> A general-purpose spectroscopy analysis toolkit for Julia.
>
> SpectroscopyTools provides peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, chirp correction, and unit conversions for spectroscopic data. While it works with any spectroscopy discipline (FTIR, Raman, UV-vis, fluorescence), there is an emphasis on ultrafast spectroscopy, including transient absorption data types, global analysis with decay-associated spectra, and chirp correction.

**After:**
> **Spectroscopy analysis for Julia — steady-state and ultrafast.**
>
> For steady-state work, SpectroscopyTools provides peak fitting (Gaussian, Lorentzian, Voigt, Fano), baseline correction, peak detection, spectral transforms (Kramers-Kronig, Tauc, Kubelka-Munk), and unit conversions for FTIR, Raman, and UV-vis data. For ultrafast work, it provides typed data structures (`TATrace`, `TASpectrum`, `TAMatrix`), global fitting with decay-associated spectra, IRF-convolved exponential fitting, and chirp correction for broadband pump-probe experiments.

Both sides of the paragraph get parallel length and specificity — equal framing footing.

### 2. Quick Start reorder

Replace the Quick Start blocks in `README.md:22–43` and `docs/src/index.md:14–35`.

Lead with the ultrafast workflow (the more distinctive capability); follow with 1D peak fitting and baseline correction. The examples remain schematic (no synthetic-data preamble) — same convention as the current README. The `docs/src/index.md` version keeps its existing runnable-synthetic-data preamble for the 1D example, since that was already a choice there.

**README (schematic):**
```julia
using SpectroscopyTools

# --- Ultrafast: kinetics with IRF-convolved biexponential ---
trace = TATrace(time, signal; wavelength=800.0)
fit = fit_exp_decay(trace; n_exp=2, irf=true)
report(fit)

# --- Ultrafast: global fit across a TA matrix ---
matrix = TAMatrix(time, wavelength, data)
gfit = fit_global(matrix; n_exp=2)   # shared τ, decay-associated spectra
das = gfit.das

# --- Steady-state: peak fitting (FTIR, Raman, UV-vis) ---
result = fit_peaks(x, y, (2000, 2100))
report(result)

# --- Steady-state: baseline correction ---
bl = correct_baseline(x, y; method=:arpls)
```

**docs/src/index.md Quick Start:** mirror the same ultrafast-first ordering. Include a short synthetic-data preamble for **both** the ultrafast and 1D examples so both blocks are pasteable-runnable — this matches the existing convention of `docs/src/index.md`, which currently has a runnable 1D Lorentzian preamble. For the ultrafast example, generate a synthetic TA trace and TAMatrix inline (e.g., a biexponential decay plus noise, and a small 2D matrix with two Gaussian decay components) so `fit_exp_decay` and `fit_global` calls will actually execute for a reader pasting the block.

### 3. Features table reorder

Replace the table in `README.md:47–60` and `docs/src/index.md:39–48`. Reorder rows into three groups — ultrafast first, steady-state second, specialty last. Current rows are preserved; only order changes (plus the addition of a new "TA data types" row at the top that surfaces the typed hierarchy).

**New order:**

| Module | Description |
|--------|-------------|
| **TA data types** | `TATrace`, `TASpectrum`, `TAMatrix` with semantic axis indexing (`m[λ=800]`, `m[t=1.0]`) |
| **Exponential decay** | Single/multi-exponential with optional IRF convolution |
| **Global fitting** | Shared parameters across traces, decay-associated spectra |
| **TA spectrum fitting** | N-peak model with ESA/GSB/SE labels and sign convention |
| **Chirp correction** | GVD detection (cross-correlation, threshold) and correction for broadband TA |
| **SVD filtering** | Matrix denoising for broadband TA data |
| **Peak fitting** | Gaussian, Lorentzian, Voigt, Pseudo-Voigt, Fano via CurveFit.jl |
| **Peak detection** | Automatic peak finding with prominence filtering |
| **Baseline correction** | arPLS, SNIP, rubber band, iModPoly, rolling ball |
| **Spectral math** | Smoothing, derivatives, band area, normalization, spectral arithmetic |
| **Transforms** | Kramers-Kronig, Kubelka-Munk, Tauc plot, SNV, Urbach tail |
| **Unit conversions** | Wavenumber, wavelength, energy, linewidth interconversion |
| **PL/Raman mapping** | Spatial maps, peak fitting, cosmic ray detection and removal |

The `docs/src/index.md` table gets the same reorder; its current row set is a subset of the README's, so its new content should match the README's new rows exactly (so the two surfaces stay in sync).

### 4. Tutorial nav reorder

Edit `docs/make.jl:14–19`. Move `Chirp Correction` to the top of the Tutorials group — it is currently the only TA-adjacent tutorial and should appear first in the nav. Rename the entry to `Chirp Correction (TA)` to signal its ultrafast context.

**Before:**
```julia
"Tutorials" => [
    "Peak Fitting (FTIR)" => "tutorials/ftir_peak_fitting.md",
    "Peak Fitting (Raman)" => "tutorials/raman_peak_fitting.md",
    "PL / Raman Map Analysis" => "tutorials/plmap_analysis.md",
    "Chirp Correction" => "tutorials/chirp_correction.md",
],
```

**After:**
```julia
"Tutorials" => [
    "Chirp Correction (TA)" => "tutorials/chirp_correction.md",
    "Peak Fitting (FTIR)" => "tutorials/ftir_peak_fitting.md",
    "Peak Fitting (Raman)" => "tutorials/raman_peak_fitting.md",
    "PL / Raman Map Analysis" => "tutorials/plmap_analysis.md",
],
```

### 5. Docs index layout section

`docs/src/index.md:50–57` describes the Diátaxis layout. Update the tutorials bullet to reflect the reordered content and the TA emphasis, e.g., lead the list of example tutorials with chirp correction. No structural change to Diátaxis sections.

### 6. README "testing internally" note

Leave `README.md:13` ("This package is currently being tested internally in our lab") unchanged. Removing it is a registration-readiness decision separate from this cosmetic pass.

## Out of scope / followups

These are flagged explicitly so they are not lost, but they are not part of this pass:

- **TA reference pages missing.** There are no reference pages for `fit_exp_decay`, `fit_global`, chirp correction API, TA types, unit conversions, or spectral transforms. Before registration, a reference-page buildout pass is needed.
- **Hero TA tutorial missing.** "Global analysis of a transient absorption matrix" would be the ideal `aha` tutorial. Blocked on real TA data the user plans to collect.
- **Additional TA how-to recipes.** Handling coherent artifact, SVD diagnostics, cleaning cosmic rays from TA data.

## Validation

- `julia --project=docs docs/make.jl` builds cleanly (no new warnings introduced by the nav reorder or content changes).
- `README.md` renders correctly on GitHub (table syntax, code fences).
- Package scope lines in `CLAUDE.md` remain consistent with the new framing (spot-check; no changes expected since `CLAUDE.md` already acknowledges the ultrafast focus).
- No changes to `src/`, no tests affected.

## Files changed

- `README.md` — headline/paragraph, Quick Start, features table.
- `docs/src/index.md` — headline/paragraph, Quick Start, features table, layout section bullet.
- `docs/make.jl` — tutorial nav ordering.
