# 3rd Infodengue-Mosqlimate Dengue Challenge (IMDC) 2026 — State-level SARIMAX

Submission repository for the **2026 3rd Infodengue-Mosqlimate Dengue Challenge**: state-level (UF) weekly probable-case forecasts for **dengue** and **chikungunya** in Brazil, produced with SARIMAX models fit independently per state (all 26 UFs, both diseases) and per focal city (15 dengue cities, 10 chikungunya cities — the Optional City-Level Challenges).

**GitHub repository:** [EzequielEBS/3rd_imdc_emap_epidematicos_sarimax_fixed](https://github.com/EzequielEBS/3rd_imdc_emap_epidematicos_sarimax_fixed)

**Related repositories:**
- [EzequielEBS/3rd_imdc_emap_epidematicos_prophet_fixed](https://github.com/EzequielEBS/3rd_imdc_emap_epidematicos_prophet_fixed) — Prophet models

## 1. Team and Contributors

**Team:** Epidemáticos

**Affiliation:** School of Applied Mathematics, Fundação Getulio Vargas (FGV/EMAp), Rio de Janeiro, Brazil

**Contributors:**

- [Eduardo Adame, M.Sc. — FGV/EMAp](https://github.com/adamesalles)
- [Ezequiel Braga, M.Sc. — FGV/EMAp](https://github.com/ezequielebs)
- [Iara Cristina, M.Sc. — FGV/EMAp](https://github.com/iaracastro)
- [Isaque Pim, Ph.D. — FGV/EMAp](https://github.com/isaquepim)

## 2. Repository Structure

```
.
├── README.md                   # this file
├── LICENSE                     # GNU GPLv3
├── pyproject.toml              # Python dependencies (Poetry), used only by sub_pred.r via reticulate
├── .gitattributes              # marks *.csv / *.gz as Git LFS objects
├── .gitignore
│
├── data_prep/                  # data preparation pipeline (R), run in this order
│   ├── prep_climate.r              # 1. fills climate gaps for island municipalities
│   ├── merge_dengue.r               # 2a. joins raw dengue cases with climate/env/ocean/population covariates
│   ├── merge_chikungunya.r          # 2b. same as above, for chikungunya
│   ├── handle_na.r                  # 3. interpolates missing ENSO/IOD/PDO values
│   ├── agg_data_uf.r                # 4. aggregates city-level data to state (UF) level + lag/rolling features
│   ├── sel_cities.r                 # optional: extracts a few focal cities for diagnostics
│   ├── summarise_sel_cities.r       # placeholder, currently empty / unused
│   └── README.md                    # details on each script above
│
├── raw_data/                   # original input datasets (Git LFS)
│   ├── dengue.csv.gz, chikungunya.csv.gz        # case counts by geocode/epiweek + train/target split flags
│   ├── climate.csv.gz                            # weather variables by geocode/epiweek
│   ├── environ_vars.csv.gz                       # Köppen climate classification, biome (static, by geocode)
│   ├── ocean_climate_oscillations.csv.gz         # ENSO, IOD, PDO indices by epiweek (national scale)
│   ├── datasus_population_2001_2025.csv.gz       # DATASUS population by geocode/year
│   ├── map_regional_health.csv                   # geocode → regional/macroregional health mapping
│   ├── shape_muni.gpkg, shape_regional_health.gpkg, shape_macroregional_health.gpkg  # spatial boundaries (currently unused as model inputs)
│   ├── forecasting_climate.csv.gz                # forecasted climate variables (currently unused by the pipeline)
│   └── *_update_2026.csv.gz                      # optional incremental extracts (dengue/chikungunya/climate/ocean_climate_oscillations) that extend coverage to a later epiweek; merged in automatically by data_prep/ when present — see data_prep/README.md, "Incremental _update_2026 raw files"
│
├── processed_data/              # output of data_prep/ — the modeling-ready tables
│   ├── climate/climate.csv.gz
│   ├── dengue/
│   │   ├── dengue_merged.csv.gz             # city-level, all covariates merged
│   │   ├── dengue_<UF>_agg.csv.gz           # one file per state — read directly by the SARIMAX scripts
│   │   └── sel_cities/                      # focal-city extracts (diagnostics only)
│   └── chikungunya/                         # same structure as dengue/
│
└── sarimax/                     # modeling code and results
    ├── README.md                            # modeling workflow and how-to-run instructions
    ├── src/
    │   ├── utils.r                          # shared helpers: fitting, covariate selection, metrics, PCA
    │   ├── model_sel.r                      # cross-validated covariate + SARIMAX-order selection (training)
    │   ├── fit.r                            # refits the best model per state and generates forecasts
    │   └── sub_pred.r                       # uploads forecasts to the Mosqlimate platform
    └── results/
        ├── concluded_states_dengue.csv, concluded_states_chikungunya.csv  # resumability checkpoints
        ├── metrics/                          # per-state/disease CV metrics, incl. best_wis_*_all_states.csv
        └── preds/                            # final forecast CSVs, one per state x target window
```

## 3. Libraries and Dependencies

**R** (used throughout `data_prep/` and `sarimax/src/`):

| Package | Used for |
|---|---|
| `tidyverse` (`dplyr`, `tidyr`, `purrr`, `readr`) | data wrangling |
| `lubridate`, `aweek` | date and epiweek handling |
| `zoo` | linear interpolation of missing ocean-climate index values |
| `runner` | rolling-window (3/6/9/12-month) covariate means |
| `sf`, `geobr` | municipality geometries/centroids, used to gap-fill island climate data |
| `forecast` | `Arima()` / `forecast()` — the SARIMAX engine |
| `scoringutils` | Weighted Interval Score (WIS) and interval-coverage scoring |
| `parallel`, `pbapply` | parallelized, progress-tracked grid search |
| `reticulate`, `data.table` | Python interop and fast I/O for the submission step |
| `dotenv` | loads the Mosqlimate API key from a local `.env` file (never hardcoded) |

**Python** (declared in `pyproject.toml`, only used via `reticulate` in `sarimax/src/sub_pred.r` to submit forecasts):

- `python = ">=3.10,<3.13"`
- `mosqlient==1.5.2` — Mosqlimate Predictions Registry client
- `epiweeks`, `python-dotenv`
- `jupyter`, `pmdarima` are declared but not currently invoked anywhere in this repository's code path — they are leftovers from the challenge template and can be removed if unused in your own workflow.

Install the R packages with, e.g., `install.packages(c("tidyverse","lubridate","aweek","zoo","runner","sf","geobr","forecast","scoringutils","pbapply","reticulate","data.table","dotenv"))`, and the Python side with `poetry install` (or `pip install mosqlient==1.5.2 epiweeks python-dotenv`).

## 4. Data and Variables

### Datasets

All disease and most covariate data are distributed by the Mosqlimate platform (Infodengue feed, probable cases — the `casprov`/`target_*` columns already reflect this); population comes from DATASUS.

- **Case data** — `raw_data/dengue.csv.gz`, `raw_data/chikungunya.csv.gz`: weekly probable case counts by municipality (`geocode`) and `epiweek`, with the official `train_1..4` / `target_1..4` split indicators used for both retrospective validation and the real submission (see Section 6).
- **Climate** — `raw_data/climate.csv.gz`: minimum/mean/maximum temperature, precipitation, atmospheric pressure, relative humidity, thermal range, and rainy-day count, by municipality and epiweek.
- **Environmental attributes** — `raw_data/environ_vars.csv.gz`: Köppen climate classification and biome, static per municipality.
- **Ocean-climate indices** — `raw_data/ocean_climate_oscillations.csv.gz`: ENSO, IOD, and PDO, indexed by epiweek (national scale, not municipality-specific).
- **Population** — `raw_data/datasus_population_2001_2025.csv.gz`: municipality population by year (2026 obtained by carrying the 2025 value forward — see `data_prep/merge_dengue.r`).
- **Regional-health mapping and shapefiles** — `raw_data/map_regional_health.csv` and the `shape_*.gpkg` files are read in by the merge scripts / available for spatial work but are **not** used as direct model covariates in the current pipeline; `sf`/`geobr` are used only in `data_prep/prep_climate.r` to fill climate gaps for island municipalities by borrowing their nearest mainland neighbor's series.
- `raw_data/forecasting_climate.csv.gz` is included but currently unused by any script.
- **Incremental updates** — when Infodengue/Mosqlimate publishes a later data drop (extending coverage to a more recent epiweek, or revising already-published weeks), it's dropped in as `raw_data/<name>_update_2026.csv.gz` (currently supported for `dengue`, `chikungunya`, `climate`, and `ocean_climate_oscillations`) rather than replacing the base file in place. `data_prep/prep_climate.r`, `merge_dengue.r`, and `merge_chikungunya.r` each merge the update into its base file automatically when present: the update's row wins on any `(geocode, epiweek)` (or `date`, for the ocean-index table) it shares with the base file, rows only in the update are appended, and untouched base rows are left alone. If no update file exists, these scripts behave exactly as if this mechanism weren't there. See `data_prep/README.md` for the full merge semantics.

### Pre-processing

Implemented in `data_prep/` (see `data_prep/README.md` for full detail), in order:

1. `prep_climate.r` — fills missing climate series for island municipalities.
2. `merge_dengue.r` / `merge_chikungunya.r` — left-joins case data with climate, environmental, ocean-climate, and population covariates.
3. `handle_na.r` — linearly interpolates the remaining gaps in ENSO/IOD/PDO.
4. `agg_data_uf.r` — sums cases and averages weather covariates from city to state (UF) level, and engineers 4/8/12/16-week lagged and 3/6/9/12-month rolling-mean versions of every weather covariate. Writes the final `processed_data/<disease>/<disease>_<UF>_agg.csv.gz` files used for modeling. Espírito Santo (ES) is excluded.

Inside model fitting itself, covariates are additionally standardized (z-scored) using the **training split's mean/SD only**, to avoid leaking information from the forecast window into the scaling — see `fit_sarimax()` in `sarimax/src/utils.r` (around the "Build & standardize regressor matrices" step).

### Variable selection

Variable selection is fully automated per state/disease inside `run_model_selection()` (`sarimax/src/model_sel.r`), using only training-split rows at every filtering step to avoid leakage:

1. **`get_candidates()`** (`sarimax/src/utils.r`) enumerates the contemporaneous, lagged, and rolling-mean versions of the 8 base covariates (`temp_med_mean`, `precip_med_mean`, `rel_humid_med_mean`, `thermal_range_mean`, `rainy_days_mean`, `enso`, `iod`, `pdo`).
2. **`filter_low_variance()`** drops near-constant covariates (SD below a threshold).
3. **`filter_by_correlation()`** drops covariates whose absolute correlation with `log1p(cases)` is below a minimum threshold.
4. **Dimensionality reduction** — by default (`pca = TRUE`): **`pca_all()`** reduces the surviving local-weather covariates to the smallest number of principal components that explain ≥ 90% of variance in the most restrictive of the four CV folds. ENSO/IOD/PDO are deliberately excluded from the PCA and tested separately with **`filter_redundant_indices()`**, which keeps an index only if it still correlates with the response after the local-weather PCs' linear effect is removed — so ocean-climate signal isn't silently absorbed (or discarded) by the PCA. Candidate formulas use the leading 1..k PCs, with and without the surviving indices.
   An alternative, non-PCA path (`pca = FALSE`) is also implemented: **`select_best_per_variable()`** keeps only the best lag/rolling version of each base variable, and **`build_covariate_combinations()`** enumerates de-correlated covariate subsets (pairwise |correlation| below a threshold) up to a maximum size, from which a random sample is grid-searched.
5. **`run_grid_search()`** (`sarimax/src/utils.r`) cross-validates every candidate formula against a grid of SARIMAX `(p,d,q)(P,D,Q)` orders across the four train/target splits, and **`run_model_selection()`** picks the formula/order combination with the lowest mean cross-validated score (WIS by default).

## 5. Model Training

Each state and disease is modeled independently with `forecast::Arima()` (R package `forecast`), called through `fit_sarimax()` / `run_grid_search()` in `sarimax/src/utils.r` and orchestrated by `run_model_selection()` in `sarimax/src/model_sel.r`. The response is `log1p(cases)` (back-transformed with `expm1()` at forecast time and clipped at zero); the regressors are the PCA components / de-correlated covariate subset selected as described in Section 4, standardized using **training-split mean/SD only**.

**Hyperparameter search.** For every state and disease, and for every candidate covariate formula produced by the variable-selection step, `run_grid_search()` cross-validates a grid of SARIMAX `(p,d,q)(P,D,Q)` orders bounded by `max_order` (default `p <= 2, d <= 1, q <= 2, P <= 1, D = 0, Q <= 1`, seasonal period 52 weeks), fit with `method = "CSS-ML"` and `optim.method = "BFGS"`. Each (formula, order) combination is fit on each of the four `train_i` windows and scored against the matching `target_i` window (see Section 6) using `compute_metrics()` -- Weighted Interval Score (WIS), MAE, RMSE, MAPE, and empirical interval coverage at the 50/80/90/95% levels. `run_model_selection()` ranks every combination by mean cross-validated score (WIS by default, configurable via its `metric` argument) and keeps, per state and disease, the single formula/order pair with the lowest mean score. The full ranking is written to `sarimax/results/metrics/metrics_all_formulas_<disease>_<state>.csv`; the winners (one row per state) are written to `sarimax/results/metrics/best_wis_<disease>_all_states.csv` and are what `sarimax/src/fit.r` reads to know which formula/order to refit.

**Where the code lives and how to run it.** From the repository root, after the data-prep pipeline (Section 4) has produced `processed_data/<disease>/<disease>_<UF>_agg.csv.gz` for every state:

```r
# 1. Training / model selection for every state (and a few focal cities), both
#    diseases. Long-running (grid search x 26 states x 2 diseases). Safe to
#    interrupt and re-run: finished states are skipped via
#    results/concluded_states_*.csv.
source("sarimax/src/model_sel.r")

# 2. Refit the winning formula/order per state on all 4 train/target windows
#    and write one forecast CSV per state x window to results/preds/.
source("sarimax/src/fit.r")

# 3. Submit the forecasts to the Mosqlimate Predictions Registry.
source("sarimax/src/sub_pred.r")
```

See `sarimax/README.md` for a file-by-file breakdown and for instructions on calling `run_model_selection()` directly on a single state instead of running the full loop.

## 6. Data Usage Restriction

The IMDC rules require that each validation split's forecast -- and, separately, the Forecast Phase submission -- be produced using only data available up to that window's own EW25 cutoff. This repository encodes that rule directly in the data rather than leaving it to be remembered at modeling time: every `processed_data/<disease>/<disease>_<UF>_agg.csv.gz` file carries four pairs of boolean columns, `train_1..train_4` and `target_1..target_4` (supplied upstream by the Mosqlimate/Infodengue data feed and carried through unchanged by `data_prep/agg_data_uf.r`), one pair per Validation Test (challenge rules, Section 5.1):

- **`train_1`-`train_4` / `target_1`-`target_4`** are the four Validation Test splits. All four are retrospective backtests that exist to cross-validate the modeling pipeline against past/current seasons (this is what Section 5's grid search scores models on) -- including split 4: Validation Test 4 is a fixed window (the 2025-2026 season -- train through EW25 2025, forecast EW41 2025..EW40 2026), the same kind of already-defined backtest as splits 1-3, not a placeholder for "whatever this pipeline is currently submitting."
- **The Forecast Phase** (challenge rules, Section 5.2) is a fifth, separate window -- the 2026-2027 season, train through EW25 2026, forecast EW41 2026..EW40 2027 -- with no `train_i`/`target_i` counterpart in the processed data at all. **It is not Validation Test 4 under another name.** It is computed directly from the rules (`SUBMISSION_FORECAST_YEAR`/`FORECAST_WINDOW`, configured at the top of `fit.r`) rather than read from any upstream flag column, and is fit/submitted under its own `target_forecast` identifier (`pred_<disease>_<id>_target_forecast.csv`), completely independent of `target_4` -- fitting or submitting one never touches or overwrites the other, and both can exist on disk (and in the Mosqlimate registry) at the same time.

The unified `fit_sarimax()` (`sarimax/src/fit.r` calls it via `window_for_split()`, which resolves each split -- `1`..`4`, or the separate string `"forecast"` -- to its own `train_start`/`train_end`/`forecast_start`/`forecast_end` bounds: splits 1-4 always via the fixed challenge-rule formula for that test's season, `"forecast"` via the `FORECAST_WINDOW` computed from `SUBMISSION_FORECAST_YEAR`) is the function used for every validation split and the Forecast Phase alike. Two mechanisms keep every split leak-free with respect to its own cutoff:

1. **No future case data.** It fits `Arima()` only on `data[epiweek >= train_start & epiweek <= train_end, ]` -- rows after each split's cutoff are never passed to the model, so no future case counts can influence the fit.
2. **No future covariate data.** Each split's forecast horizon extends past the last epiweek for which real weather/ocean-index values exist as of that split's own cutoff. Rather than requiring those unobserved values (or borrowing another year's realized trajectory as a proxy), `fit_sarimax()` calls `forecast_covariates()` (`sarimax/src/forecast_covariates.r`) to produce a genuine, out-of-sample forecast for every regressor the selected formula uses (PCA components and/or the surviving ocean-climate indices): for each genuinely distinct series, every future epiweek is forecast as the mean of that same calendar week-of-year across the most recent 5 years of the series' own history up to that split's `train_end` (a seasonal, "epiweek-of-year" climatology; the 5-year window is configurable via `n_years`). At state level, some candidate columns aren't genuinely distinct series at all -- `data_prep/agg_data_uf.r` derives lagged (`_lag4/8/12/16`) and rolling-mean (`_mean_3/6/9/12mo`) versions of every weather aggregate as fixed transforms (a shift, or a trailing average) of that one contemporaneous series -- so `forecast_covariates()` forecasts the underlying series exactly once and mechanically re-derives its lag/rolling columns from that single forecast with the same arithmetic `agg_data_uf.r` used to build them from real data, rather than independently forecasting up to 9 columns that are mechanically one series. This keeps the covariate side of every validation split's forecast ex-ante the same way the case-count model itself is -- no real future values and no other year's realized trajectory are read, only out-of-sample predictions computed from data through that split's own cutoff. `run_grid_search()` (`sarimax/src/utils.r`) uses the same `fit_sarimax()` during model selection, precomputing each fold's covariate forecasts once (rather than once per formula/order combination) so that CV scoring is evaluated under the same genuinely-forecasted-covariate regime the final submission uses. Every one of these root-series forecasts is also cached in a persistent, shared table (`processed_data/covariate_forecasts_<disease>_<state|city>.csv`, keyed on unit/series/fold window; see "Covariate forecasting" in `sarimax/README.md`), so `model_sel.r`, `fit.r`, and `evaluate_covariate_forecasts.r` reuse the same forecast for a given unit/fold instead of each independently recomputing it -- `sarimax/src/build_covariate_forecast_table.r` can warm this table for every unit up front.

## 7. Predictive Uncertainty

Prediction intervals come from `forecast::forecast()`'s bootstrap-simulation path (`bootstrap = TRUE`, `npaths = 1000`, set in `fit_sarimax()` in `sarimax/src/utils.r`) rather than the default Gaussian/normal-theory intervals: 1,000 forecast paths are simulated from the fitted SARIMAX model and empirical quantiles of those paths are taken as the interval bounds. This is more robust to the non-normality introduced by the `log1p`/`expm1` response transform than a Gaussian approximation would be. Intervals are produced at four nominal coverage levels -- 50%, 80%, 90%, and 95% (the `levels` argument, default `c(50, 80, 90, 95)`) -- and reported as `lower_<level>`/`upper_<level>` columns alongside the point forecast `pred` in every `sarimax/results/preds/pred_<disease>_<UF>_target_<n>.csv` file; all bounds are clipped at zero after back-transforming with `expm1()`.

Because `forecast()`'s point forecast (`mean`, an analytic recursive calculation) and its bootstrap intervals (empirical quantiles of the 1,000 simulated paths) come from two different mechanisms, they aren't automatically guaranteed to agree -- for a low-signal series (sparse counts, no surviving covariates), the point forecast can occasionally fall outside its own narrowest interval, or an interval can come out wider than a nominally broader one. `enforce_interval_nesting()` (`sarimax/src/utils.r`, called at the end of every `fit_sarimax()` invocation) corrects this before anything is written to disk: it widens -- never narrows -- whichever bounds are violated, from the narrowest level outward, so every row satisfies `lower_95 <= lower_90 <= lower_80 <= lower_50 <= pred <= upper_50 <= upper_80 <= upper_90 <= upper_95`. This matters beyond internal consistency: the Mosqlimate submission API rejects any prediction row where the bounds aren't properly ordered, so this guarantee is what keeps every submitted row valid.

Interval calibration is checked, not just assumed: `compute_metrics()` (`sarimax/src/utils.r`) computes empirical coverage at each of the four levels (the fraction of held-out `target_i` observations that actually fall inside the predicted interval) every time a candidate model is cross-validated. Critically, model selection (Section 5) ranks candidates by mean Weighted Interval Score rather than by point-forecast error alone -- WIS is a proper score that rewards sharp *and* well-calibrated intervals simultaneously (see Reference 1 below), so the chosen model per state/disease is selected for probabilistic quality, not just point accuracy.

## 8. References

This submission does not implement one specific previously published forecasting method; it is a per-state/disease SARIMAX pipeline built on the `forecast` R package, with cross-validated covariate selection and Weighted-Interval-Score-based model ranking. The following works underpin the modeling methodology and the data/evaluation/submission infrastructure this repository builds on:

1. Hyndman, R. J., & Khandakar, Y. (2008). Automatic time series forecasting: the forecast package for R. *Journal of Statistical Software*, 27(3), 1-22. https://doi.org/10.18637/jss.v027.i03 -- the `forecast` package providing `Arima()`/`forecast()`, the SARIMAX fitting and bootstrap-interval engine used throughout `sarimax/src/`.
2. Bracher, J., Ray, E. L., Gneiting, T., & Reich, N. G. (2021). Evaluating epidemic forecasts in an interval format. *PLOS Computational Biology*, 17(2), e1008618. https://doi.org/10.1371/journal.pcbi.1008618 -- defines the Weighted Interval Score that `compute_metrics()` and `run_model_selection()` use to rank candidate models and assess interval calibration (Sections 5 and 7).
3. Ganem, F., et al. (2024). Mosqlimate: a platform to providing automatable access to data and forecasting models for arbovirus disease. *arXiv:2410.18945*. https://arxiv.org/abs/2410.18945 -- describes the Mosqlimate data store and Predictions Registry that `raw_data/` originates from and that `sarimax/src/sub_pred.r` submits forecasts to (via the `mosqlient` package).
4. Araujo, E. C., Carvalho, L. M., Ganem, F., Coelho, F. C., et al. (2025). Leveraging probabilistic forecasts for dengue preparedness and control: the 2024 Dengue Forecasting Sprint in Brazil. *medRxiv* 2025.05.12.25327419. https://doi.org/10.1101/2025.05.12.25327419 -- reports on the first edition of the Infodengue-Mosqlimate Dengue Challenge that this repository's submission (the 3rd edition) continues.
