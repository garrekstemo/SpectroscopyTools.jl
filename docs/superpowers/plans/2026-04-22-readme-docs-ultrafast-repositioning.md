# README and Docs Ultrafast Repositioning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reposition `README.md` and `docs/src/index.md` to put ultrafast / time-resolved spectroscopy forward as a primary capability alongside steady-state, and reorder the Documenter tutorial nav to surface the TA-adjacent tutorial first.

**Architecture:** Pure documentation edits. Three files change: `README.md`, `docs/src/index.md`, `docs/make.jl`. No source code, no exports, no tests affected. Verification is a clean local Documenter build and visual inspection of rendered Markdown.

**Tech Stack:** Documenter.jl for the docs site, Markdown for README and `docs/src/*.md`, Julia 1.10+ for the doc build.

**Spec:** `docs/superpowers/specs/2026-04-22-readme-docs-ultrafast-repositioning-design.md`

---

## Task 1: Update `README.md`

**Files:**
- Modify: `README.md` (lines 9–11 headline/intro, lines 22–43 Quick Start, lines 47–60 features table)

- [ ] **Step 1: Read current `README.md`**

Run: `cat README.md`
Expected: confirm current line numbers for the three sections we are editing (headline/intro paragraph around lines 9–11, Quick Start code block around lines 22–43, features table around lines 47–60). If the file has drifted, re-locate by section heading before editing.

- [ ] **Step 2: Replace headline + opening paragraph**

Replace the current tagline and intro paragraph:

```markdown
A general-purpose spectroscopy analysis toolkit for Julia.

SpectroscopyTools provides peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, chirp correction, and unit conversions for spectroscopic data. While it works with any spectroscopy discipline (FTIR, Raman, UV-vis, fluorescence), there is an emphasis on ultrafast spectroscopy, including transient absorption data types, global analysis with decay-associated spectra, and chirp correction.
```

with:

```markdown
**Spectroscopy analysis for Julia — steady-state and ultrafast.**

For steady-state work, SpectroscopyTools provides peak fitting (Gaussian, Lorentzian, Voigt, Fano), baseline correction, peak detection, spectral transforms (Kramers-Kronig, Tauc, Kubelka-Munk), and unit conversions for FTIR, Raman, and UV-vis data. For ultrafast work, it provides typed data structures (`TATrace`, `TASpectrum`, `TAMatrix`), global fitting with decay-associated spectra, IRF-convolved exponential fitting, and chirp correction for broadband pump-probe experiments.
```

- [ ] **Step 3: Replace the Quick Start code block**

Replace the entire fenced code block beginning at `using SpectroscopyTools` in the Quick Start section (currently lines 24–43) with this schematic block — ultrafast first, then steady-state:

````markdown
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
````

Note: variables like `time`, `signal`, `data`, `wavelength`, `x`, `y` are assumed pre-existing in the user's session — this matches the schematic convention of the current README. The docs/src/index.md version in Task 2 will be runnable.

- [ ] **Step 4: Replace the features table**

Replace the entire Markdown table currently at `README.md` lines 48–60 with this table (ultrafast rows first, steady-state second, specialty last; includes a new "TA data types" row at the top):

```markdown
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
```

- [ ] **Step 5: Verify render**

Run: `git diff README.md | head -200`
Expected: the diff shows exactly the three replacements (headline/intro, Quick Start, features table), no unrelated changes. The changed sections use valid GitHub-flavored Markdown: code fence uses triple backticks with `julia` language tag, table uses `|` separators with a header row, bold uses `**`.

- [ ] **Step 6: Commit**

```bash
git add README.md
git commit -m "docs: reposition README around ultrafast and steady-state spectroscopy"
```

---

## Task 2: Update `docs/src/index.md`

**Files:**
- Modify: `docs/src/index.md` (lines 3–5 intro, lines 14–35 Quick Start, lines 39–48 features table, lines 50–57 layout section)

- [ ] **Step 1: Read current `docs/src/index.md`**

Run: `cat docs/src/index.md`
Expected: confirm the four sections to edit: headline/intro (lines 3–5), Quick Start (lines 14–35), features table (lines 39–48), "Documentation Layout" bullets (lines 50–57).

- [ ] **Step 2: Replace headline + opening paragraph**

Replace the current line-3 headline and line-5 paragraph:

```markdown
# SpectroscopyTools.jl

A general-purpose spectroscopy analysis toolkit for Julia.

SpectroscopyTools provides peak fitting, baseline correction, exponential decay fitting with IRF deconvolution, chirp correction, and unit conversions for spectroscopic data. It is designed to work with any spectroscopy discipline --- ultrafast, FTIR, Raman, UV-vis, fluorescence --- while keeping a minimal dependency footprint.
```

with (keep the `# SpectroscopyTools.jl` heading, replace the tagline and intro paragraph):

```markdown
# SpectroscopyTools.jl

**Spectroscopy analysis for Julia — steady-state and ultrafast.**

For steady-state work, SpectroscopyTools provides peak fitting (Gaussian, Lorentzian, Voigt, Fano), baseline correction, peak detection, spectral transforms (Kramers-Kronig, Tauc, Kubelka-Munk), and unit conversions for FTIR, Raman, and UV-vis data. For ultrafast work, it provides typed data structures (`TATrace`, `TASpectrum`, `TAMatrix`), global fitting with decay-associated spectra, IRF-convolved exponential fitting, and chirp correction for broadband pump-probe experiments.
```

- [ ] **Step 3: Replace the Quick Start code block with runnable synthetic-data version**

Replace the existing Quick Start code block (currently beginning `using SpectroscopyTools` around line 16 and ending with `wavelength_to_wavenumber(1500u"nm")` around line 34) with this runnable, pasteable version:

````markdown
```julia
using SpectroscopyTools
using CurveFitModels  # for lineshape functions
using Unitful

# ============================================================
# Ultrafast example: synthetic biexponential decay + TA matrix
# ============================================================
time = collect(range(0, 20.0, length=400))
signal = 0.6 .* exp.(-time ./ 0.5) .+ 0.4 .* exp.(-time ./ 5.0) .+ 0.01 .* randn(length(time))

trace = TATrace(time, signal; wavelength=800.0)
fit = fit_exp_decay(trace; n_exp=2)
report(fit)

# Synthetic 2D TA matrix: rows = time, cols = wavelength
matrix_time = collect(range(0, 15.0, length=200))
matrix_wl   = collect(range(500.0, 700.0, length=40))
matrix_data = zeros(length(matrix_time), length(matrix_wl))
for (j, λ) in pairs(matrix_wl)
    a_fast = exp(-((λ - 550) / 20)^2)
    a_slow = -0.8 * exp(-((λ - 640) / 30)^2)
    matrix_data[:, j] = a_fast .* exp.(-matrix_time ./ 0.8) .+ a_slow .* exp.(-matrix_time ./ 4.0)
end
matrix_data .+= 0.005 .* randn(size(matrix_data))

matrix = TAMatrix(matrix_time, matrix_wl, matrix_data)
gfit = fit_global(matrix; n_exp=2)
spectra = das(gfit)   # n_exp × n_wavelengths decay-associated spectra

# ============================================================
# Steady-state example: noisy Lorentzian peak + baseline
# ============================================================
x = collect(range(1900, 2200, length=600))
y = lorentzian([0.45, 2062.0, 24.0, 0.01], x) .+ 0.002 .* randn(length(x))

result = fit_peaks(x, y, (2000, 2150))
report(result)

bl = correct_baseline(x, y; method=:arpls, λ=1e5)
y_corrected = bl.y

# Unit conversions
wavelength_to_wavenumber(1500u"nm")  # -> wavenumber in cm⁻¹
decay_time_to_linewidth(1.0u"ps")    # -> linewidth in cm⁻¹
```
````

Variable names in the TAMatrix synthetic data (`matrix_time`, `matrix_wl`, `matrix_data`) are prefixed so they do not shadow the earlier `time` / `wavelength` scalars used for the `TATrace`. This keeps the whole block pasteable-runnable.

- [ ] **Step 4: Replace the features table**

Replace the existing table (currently lines 39–48) with the same 13-row table from Task 1 Step 4 (identical content so README and docs home stay in sync):

```markdown
| Module | What it does |
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
```

Note: the column header reads "What it does" in `index.md` vs "Description" in README; keep the `docs/src/index.md` version's "What it does" column header, since it matches the prose tone of the docs home.

- [ ] **Step 5: Update the Documentation Layout bullet for Tutorials**

Replace the Tutorials bullet in the "Documentation Layout" section (currently around line 54) to lead with chirp correction (TA) in the example list:

Before:
```markdown
- **Tutorials** --- step-by-step walkthroughs for complete workflows (FTIR/Raman peak fitting, PL/Raman map analysis, chirp correction)
```

After:
```markdown
- **Tutorials** --- step-by-step walkthroughs for complete workflows (chirp correction for broadband TA, FTIR/Raman peak fitting, PL/Raman map analysis)
```

- [ ] **Step 6: Verify render**

Run: `git diff docs/src/index.md | head -300`
Expected: the diff shows exactly four changes — intro paragraph, Quick Start block, features table, tutorials bullet. No unrelated changes.

- [ ] **Step 7: Commit**

```bash
git add docs/src/index.md
git commit -m "docs: reposition docs home around ultrafast and steady-state"
```

---

## Task 3: Reorder tutorial nav in `docs/make.jl`

**Files:**
- Modify: `docs/make.jl:14–19`

- [ ] **Step 1: Read current `docs/make.jl`**

Run: `cat docs/make.jl`
Expected: confirm the `"Tutorials" => [...]` block currently lists `Peak Fitting (FTIR)`, `Peak Fitting (Raman)`, `PL / Raman Map Analysis`, `Chirp Correction` in that order.

- [ ] **Step 2: Reorder tutorials**

Replace:
```julia
"Tutorials" => [
    "Peak Fitting (FTIR)" => "tutorials/ftir_peak_fitting.md",
    "Peak Fitting (Raman)" => "tutorials/raman_peak_fitting.md",
    "PL / Raman Map Analysis" => "tutorials/plmap_analysis.md",
    "Chirp Correction" => "tutorials/chirp_correction.md",
],
```

with:
```julia
"Tutorials" => [
    "Chirp Correction (TA)" => "tutorials/chirp_correction.md",
    "Peak Fitting (FTIR)" => "tutorials/ftir_peak_fitting.md",
    "Peak Fitting (Raman)" => "tutorials/raman_peak_fitting.md",
    "PL / Raman Map Analysis" => "tutorials/plmap_analysis.md",
],
```

- [ ] **Step 3: Verify diff**

Run: `git diff docs/make.jl`
Expected: only the four tutorial lines changed — first entry moved from last to first, and renamed from `"Chirp Correction"` to `"Chirp Correction (TA)"`.

- [ ] **Step 4: Commit**

```bash
git add docs/make.jl
git commit -m "docs: move chirp correction tutorial to top of nav and label as TA"
```

---

## Task 4: Build docs locally and verify

**Files:**
- No edits; only a doc build.

- [ ] **Step 1: Build the docs**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate(); include("docs/make.jl")'`

Expected: Documenter runs to completion. Any `missing_docs` or `cross_references` warnings that existed before this change will remain (they are set to `warnonly` in `docs/make.jl`), but **no new warnings** should appear. Pay particular attention to warnings about the Home page Markdown (new intro/table/code block) — those would indicate a rendering issue we need to fix.

If the initial `Pkg.develop` + `Pkg.instantiate` step takes too long because Julia is resolving a fresh docs environment, an equivalent shorter invocation is:

```bash
julia --project=docs docs/make.jl
```

(use this if `docs/Manifest.toml` already resolves cleanly).

- [ ] **Step 2: Inspect build output**

Expected build output includes a line like `[ Info: Deployment complete` or, for a local build without CI, the build finishes printing generated file counts.
Expected: `docs/build/index.html` exists and is freshly modified. No error lines.

Run: `ls -la docs/build/index.html`

- [ ] **Step 3: Visually verify the built HTML**

Run: `grep -c "steady-state and ultrafast" docs/build/index.html`
Expected: at least `1` — confirms the new tagline rendered into the HTML.

Run: `grep -c "Chirp Correction (TA)" docs/build/index.html`
Expected: at least `1` — confirms the renamed tutorial entry appears in the nav.

- [ ] **Step 4: Verify README renders on GitHub-flavored Markdown**

Spot-check by visual inspection. Run: `head -70 README.md`
Expected: headline bold, intro paragraph parallel, Quick Start fenced code block with `julia` tag, features table has 13 rows plus a header row.

- [ ] **Step 5: If everything renders cleanly, no commit needed**

The doc build produces artifacts in `docs/build/` which is already gitignored. If the build succeeded and the grep checks pass, this task is complete with no new commit.

If the build introduced a new warning or the HTML is missing expected content, stop and surface the error — do not proceed with a commit that breaks docs.

---

## Self-Review Checklist (for the plan author, not the executing agent)

After writing the plan, the writing-plans skill requires a self-review pass. Results of that pass, done inline:

**1. Spec coverage.** Every section of the spec maps to a task:
- Spec §1 (headline/intro) → Task 1 Step 2 + Task 2 Step 2 ✓
- Spec §2 (Quick Start reorder) → Task 1 Step 3 + Task 2 Step 3 ✓
- Spec §3 (features table reorder) → Task 1 Step 4 + Task 2 Step 4 ✓
- Spec §4 (tutorial nav reorder) → Task 3 ✓
- Spec §5 (docs index layout bullet) → Task 2 Step 5 ✓
- Spec §6 (leave "testing internally" note alone) → no-op, satisfied by not including any edit ✓
- Spec validation (doc build cleanly) → Task 4 ✓

**2. Placeholder scan.** No "TBD", "TODO", or vague instructions. All code blocks show the exact content to paste.

**3. Type consistency.** `TATrace(time, signal; wavelength=800.0)` matches `src/types.jl:99`. `fit_exp_decay(trace; n_exp=2, irf=true)` matches `src/fitting.jl:148`. `fit_global(matrix; n_exp=2)` matches `src/fitting.jl:496`. `das(gfit)` matches `src/types.jl:852` and is exported in `src/SpectroscopyTools.jl:61`. `TAMatrix(time, wavelength, data)` matches `src/types.jl:687`. Variable names inside the docs Quick Start preamble (`matrix_time`, `matrix_wl`, `matrix_data`) are fresh and do not shadow earlier `time` / `wavelength` scalars.
