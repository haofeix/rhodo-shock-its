# Opposing Responses of Rhodobacteraceae to Two Marine Environmental Shocks

**Authors:** Yirong Xu, Haofei Xu (Washington University in St. Louis)

## Overview

Interrupted time-series (ITS) + event-study analysis of two exogenous marine
shocks on the abundance of *Rhodobacteraceae* in a 93-day daily 16S rRNA time
series from Nahant, MA (Martin-Platero et al. 2018, *Nat. Commun.*
9:266).

- **Shocks.** Coastal upwelling (days 215–221, 7 d); Hurricane Earl
  (days 247–249, 3 d).
- **Main outcome.** Log relative abundance of Rhodobacteraceae.
- **Headline result.** Rhodobacteraceae responds in opposite directions to the
  two shocks: –2.51 σ peak during upwelling; +3.25 σ delayed bloom peaking at
  event time k = +6 after hurricane. Pre-trend clean on both sides
  (0/6 |t| > 1.96 violations).

## Directory layout

```
final_project/
├── README.md                        <- this file
├── code/                            <- analysis scripts (run 01 → 06)
│   ├── 01_clean.R
│   ├── 02_baseline_acf.R
│   ├── 03_its_main.R
│   ├── 04_event_study.R
│   ├── 05_headline_fig.R
│   └── 06_robustness.R
├── data/
│   └── raw/                         <- original inputs (do not modify)
│       ├── daily_alpha_diversity-1.0.xlsx
│       ├── daily_mean_relative_abundance.xlsx
│       └── 环境数据.xlsx             <- environmental metadata
└── poster/                          <- final poster (LaTeX source + PDF)
    ├── poster.pdf                   <- 48 × 36 in landscape deliverable
    ├── poster.tex
    ├── poster.bib
    ├── beamerthemegemini.sty        <- Gemini theme (A. Athalye)
    ├── beamercolorthemewashu.sty    <- WashU brand color palette
    ├── figures/                     <- fig08, fig11 embedded in poster
    └── logos/                       <- WashU logo
```

After running the pipeline, `data/clean/` and `output/` will also appear
(both regenerable from `code/` + `data/raw/`).

## Requirements

- **R** ≥ 4.0 with packages:
  `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggtext`, `splines`, `sandwich`,
  `lmtest`, `patchwork`, `readxl`
- **LaTeX** with the `lualatex` engine, plus:
  `beamerposter`, `type1cm`, `raleway`, `lato`
  (install with `tlmgr install beamerposter type1cm raleway lato`)

## Reproduce the analysis

From the project root, run the six R scripts in order:

```bash
Rscript code/01_clean.R
Rscript code/02_baseline_acf.R
Rscript code/03_its_main.R
Rscript code/04_event_study.R
Rscript code/05_headline_fig.R
Rscript code/06_robustness.R
```

Full pipeline < 1 minute on a modern laptop. Each script prints a console
summary and writes figures/tables to `output/` (auto-created on first run).

| Script | Purpose | Key outputs |
|---|---|---|
| `01_clean.R` | Merge 3 raw Excel files; aggregate OTUs to 5 indicator taxa | `data/clean/daily_panel.csv`, `output/figs/fig01–02.pdf` |
| `02_baseline_acf.R` | Log-transform taxa; fit baseline trend; residual ACF for HAC lag selection | `data/clean/panel_transformed.csv`, `output/figs/fig03–04.pdf`, `output/tables/tab_pre_baseline.md` |
| `03_its_main.R` | Main ITS regression, HAC SE, primary/secondary outcome tables | `output/figs/fig05_coef_forest.pdf`, `output/tables/tab_main_primary.md`, `tab_main_secondary.md` |
| `04_event_study.R` | Dynamic event studies for both shocks × 8 outcomes | `output/figs/fig06–07.pdf`, `output/tables/tab_pretrend.md` |
| `05_headline_fig.R` | Poster headline figure (Rhodo × 2 shocks) | `output/figs/fig08_rhodo_headline.pdf` |
| `06_robustness.R` | Specification curve: Rhodo × Upwelling across 9 specs | `output/figs/fig11_spec_curve.pdf` |

## Compile the poster

The two poster figures (`fig08_rhodo_headline.pdf`, `fig11_spec_curve.pdf`)
are already in `poster/figures/`, so the poster can be rebuilt without
rerunning R:

```bash
cd poster
lualatex poster.tex
bibtex poster
lualatex poster.tex
lualatex poster.tex
```

Output: `poster/poster.pdf` (48 × 36 in landscape).

If you re-run the R pipeline first, copy the refreshed figures before
compiling:

```bash
cp output/figs/fig08_rhodo_headline.pdf output/figs/fig11_spec_curve.pdf poster/figures/
```

## Key design choices

- **Main spec is ITS, not TWFE.** Data are aggregated to daily means
  (N = 1 time series), so unit fixed effects are not identified. A
  natural cubic spline with df = 6 replaces two-way fixed effects to
  model the smooth trend.
- **No environmental controls in the main spec.** Water temperature,
  salinity, and nutrients are *post-treatment mediators* of the shocks;
  conditioning on them would trigger bad-control bias (violates the
  backdoor criterion: do not condition on colliders / post-treatment
  variables).
- **HAC standard errors with outcome-specific lag.** Newey–West lag
  chosen from the residual ACF of the trend-only fit
  (lag = 14 for Rhodobacteraceae to cover one spring–neap tidal cycle).
- **Pre-trend diagnostic uses individual |t| > 1.96 count.** HAC vcov on
  one-day event-time indicator dummies is structurally rank-deficient,
  so a joint F/Wald test is not well-defined. We report the count of
  pre-period coefficients with individual |t| > 1.96 (out of 6) and the
  max |t|.

## Data provenance

- `daily_alpha_diversity-1.0.xlsx`: daily Shannon, Richness, Simpson.
- `daily_mean_relative_abundance.xlsx`: OTU-level relative abundances with
  Greengenes taxonomy strings (used to aggregate Vibrionaceae,
  Flavobacteriaceae, Cyanobacteria, SAR11/Pelagibacteraceae, Rhodobacteraceae).
- `环境数据.xlsx`: daily environmental metadata (air/water temperature,
  salinity, atmospheric pressure, wave height, wind speed, nutrients).

Original source: Martin-Platero A.M., Cleary B., Kauffman K., et al. (2018).
"High resolution time series reveals cohesive but short-lived communities in
coastal plankton." *Nature Communications* 9:266.
