<img src="spaCraft_logo.png" align="right" height="250" alt="spaCraft logo" />

# spaCraft

> **Analysis-aware power calculation and sample-size planning for multi-sample spatial transcriptomics.**

![version](https://img.shields.io/badge/version-1.0.1-5E8CB6)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-5E8CB6)
![license](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%203)-5E8CB6)
![platforms](https://img.shields.io/badge/platforms-Visium%20%7C%20Visium%20HD%20%7C%20Stereo--seq-5E8CB6)

`spaCraft` answers the first question of every comparative spatial-transcriptomics
study — **how many samples per group do I need?** — by learning a cohort-level
generative model from a small **pilot** dataset and running a full
*generate → recover → test* Monte-Carlo pipeline. Power is reported per endpoint and
per group sample size (*K*), and summarized into a concrete sample-size recommendation.

Unlike single-cell power tools, `spaCraft` is built for the structure that defines
spatial data: it models **domain composition**, **spatial geometry**, and **spatially
correlated expression** jointly, re-discovers spatial domains in every synthetic
replicate, and evaluates the *same* tests an analyst would actually run — so the
reported power reflects the real analysis pipeline, not an idealized oracle.

---

## Why spaCraft

- **Pilot-based, single-pilot friendly.** Learns the generator from real pilot
  sections; designs with as little as one section per group are supported through
  automatic pseudo-replication (`generate_pseudo_reps`).
- **A three-layer generative model** that respects spatial biology:
  composition (baseline-anchored Binomial–logit), geometry
  (Fisher–Gaussian kernel mixture, FGKMM), and expression
  (Gaussian process with nearest-neighbor approximation, NNGP).
- **Honest *generate → recover → test* power.** Spatial domains are recovered in every
  replicate with a pilot-guided BANKSY extension (**pBANKSY**), so clustering error is
  *inside* the power estimate rather than assumed away.
- **Two complementary endpoints**, both with a minimum-effect **TREAT** margin:
  a spatially-adjusted differential-expression test (**SaLFC**) and a baseline-anchored
  compositional test (**LOR**).
- **Validity-guarded hyperparameter tuning.** `tuneSpaDesignLambdas()` maximizes a
  *regularized noncentrality* (not raw variance) under a null-calibration guard, so
  signal can never be bought by type-I inflation.
- **Calibrated, monotone reporting.** Per-group sample sizes are recommended from
  Monte-Carlo estimates smoothed with shape-constrained additive models (SCAM),
  guaranteeing power curves that are monotone in *K*.
- **Cross-platform.** Validated on **10x Visium**, **Visium HD**, and **Stereo-seq**.

---

## Installation

`spaCraft` links to compiled code (`Rcpp` / `RcppArmadillo`) and uses a few
Bioconductor packages, so install those first.

```r
# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "SingleCellExperiment", "SpatialExperiment", "spatialLIBD"))

# spaCraft (from GitHub)
# install.packages("remotes")
remotes::install_github("c16267/spaCraft")
```

Or from a local source tarball:

```r
install.packages("spaCraft_1.0.1.tar.gz", repos = NULL, type = "source")
```

**Requirements.** R ≥ 4.1.0 and a C++ toolchain (for `Rcpp`/`RcppArmadillo`).
Heavy steps are parallelized via `n_cores` (forked; falls back to serial on Windows).

---


## Input format

`createspaCraftObject()` expects a `list` of pilot samples. Each element is a `list` with:

| Field | Type | Description |
|-------|------|-------------|
| `counts` | genes × spots matrix | Raw counts (sparse or dense); **rownames required** |
| `coords` | `data.frame` | Columns `x`, `y`, and `domain` (annotated spatial domain) |
| `group` | scalar | Experimental group: control `= 0`, case `= 1` |
| `sample_id` | character | Sample identifier |

---

## Method at a glance

A `spaCraft` analysis has two phases: **fit** the generator once from pilot data,
then **sweep** sample sizes and effect sizes to map power.

**1. Cohort-level generative model (learned from the pilot).**

| Layer | What it captures | Model | Function |
|-------|------------------|-------|----------|
| Composition | Per-domain abundance, case vs. control | Baseline-anchored **Binomial–logit** GLM (intercept β₀, effect β₁) | `estimateCompositionParams()` |
| Geometry | Shape and placement of spatial domains | **Fisher–Gaussian kernel mixture (FGKMM)**, pooled across samples | `estimateGeometryParams()` |
| Expression | Domain means + spatial covariance + between-sample variance | **Gaussian process (NNGP)** with per-gene parameters | `estimateExpressionParams()` |

**2. Generate → recover → test (per Monte-Carlo replicate).**

```
simulatespaCraft()  ─►  pBANSKY()  ─►  [rearrangeSyntheticToPilot()]  ─►  SaLFC() / LOR()
  synthetic cohort       recover         morph onto pilot geometry          endpoint test
  at effect size δ, K    domains (d̂)     (SaLFC only)                       + TREAT margin τ
```

Repeating this across a grid of (`K`, effect size) yields a rejection-rate surface;
the **null row** (effect = 0) returns the empirical type-I error / FDP.

**3. Summarize → recommend.** `summary()` and `plotPowerCurve*()` fit a monotone SCAM
(monotone-increasing in *K*) to the surface and report the **minimum *K*** reaching a
target power (default 0.8) for each effect size.

---

## Quick start

```r
library(spaCraft)

## 0 — Example pilot data: full Visium human-brain (DLPFC) pilot list ----------
##     downloaded on-the-fly from the chunglab server (publicly reachable)
data_url <- "https://chunglab.bmi.osumc.edu/spaCraft/visium_human_brain_pilot_data_list.RData"

## cache the ~350 MB file so it downloads only once
cache_dir <- tools::R_user_dir("spaCraft", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dest <- file.path(cache_dir, basename(data_url))
if (!file.exists(dest) || file.info(dest)$size < 1e6) {
  options(timeout = max(3600, getOption("timeout")))       # ~350 MB: raise 60s default
  utils::download.file(data_url, destfile = dest, mode = "wb")   # mode = "wb" = binary
}

## the file stores a single object named `pilot_data_list`
load(dest)                       # -> creates `pilot_data_list` in the workspace
stopifnot(exists("pilot_data_list"))

## 1 — Pilot data  ->  spaCraft object -----------------------------------------
obj <- createspaCraftObject(pilot_data_list)   # auto pseudo-reps if K = 1

## 2 — Feature selection (domain-informative + stable/null genes) --------------
obj <- featureSelection(obj,       max_num_gene = 20)
obj <- featureSelectionStable(obj, max_num_gene = 50)

## 3 — Labeling / null / spike-in gene sets ------------------------------------
sets <- makeCustomGeneSets(
  obj,
  G_svg_base       = obj@params_expression$top_genes,
  G_stable         = obj@params_expression$stable_genes,
  target_domain    = "WM",
  reference_domain = "Layer6",
  n_de             = 5L
)

## 4 — Fit the three-layer generative model ------------------------------------
obj <- estimateCompositionParams(obj, target_domain = "WM", reference_domain = "Layer6")
obj <- estimateGeometryParams   (obj, n_cores = 4)
obj <- estimateExpressionParams (obj, target_domain = "WM", reference_domain = "Layer6",
                                 genes_to_use = sets$G_svg, n_cores = 4)

plot(obj, type = "pilot", n_per_group = 2)          # inspect the pilot geometry

## 5 — (optional) tune spatial hyperparameters as a stable power oracle --------
tuning <- tuneSpaDesignLambdas(
  obj, target_domain = "WM", reference_domain = "Layer6",
  G_DE        = sets$G_spike,
  scenario_H1 = list(DE_lfc = c(0.30, 0.35, 0.40),
                     target_prop_case = c(0.65, 0.70, 0.75)),
  B = 10, K_syn = 10, n_cores = 4
)
summary(tuning)
plot(tuning)

## 6a — Power for the spatial DE endpoint (SaLFC) ------------------------------
powSaLFC <- evaluatePowerSaLFC(
  obj, K_grid = 3:8, lfc_grid = c(0.3, 0.4, 0.5),
  G_svg = sets$G_svg, G_null = sets$G_null, G_spike = sets$G_spike,
  target_domain = "WM", reference_domain = "Layer6",
  n_sim = 30, tuning_obj = tuning, n_cores = 4
)
summary(powSaLFC, target_power = 0.8)   # min K per effect + FDR diagnostic
plotPowerCurveSaLFC(powSaLFC)

## 6b — Power for the composition endpoint (LOR) ------------------------------
powLOR <- evaluatePowerLOR(
  obj, K_grid = 3:8, delta_grid = c(0.60, 0.70, 0.80),
  genes = sets$G_svg, target_domain = "WM", reference_domain = "Layer6",
  effect_type = "target_prop_case",
  n_sim = 30, tuning_obj = tuning, n_cores = 4
)
summary(powLOR, target_power = 0.8)     # min K per effect + empirical size
plotPowerCurveLOR(powLOR)
```

> **Tip.** To prototype a single design point, run the chain directly:
>
> ```r
> sim <- simulatespaCraft(obj, n_sample_per_group = 5,
>                         scenario_settings = list(DE_lfc = 0.4, lambda_cond = 0.5),
>                         target_domain = "WM",
>                         genes_to_simulate = sets$G_svg, de_genes = sets$G_spike)
> syn <- pBANSKY(sim, obj@pilot_data, lambda = 0.5, do_hungarian = TRUE)
> syn <- rearrangeSyntheticToPilot(syn, obj@pilot_data, match_by = "hat_d")
> res <- SaLFC(syn, obj, genes = sets$G_spike,
>               target_domain = "WM", reference_domain = "Layer6",
>               domain_col = "hat_d", lfc_threshold = "pilot")
> ```

---

## The two endpoints

| | **SaLFC** — spatial differential expression | **LOR** — domain composition |
|---|---|---|
| **Question** | Is a gene's Target-vs-Reference log-fold-change different between case and control? | Does the target domain's abundance differ between case and control? |
| **Per-sample statistic** | z = x̄_T − x̄_R with spatial sampling variance τ², empirical-Bayes moderated (`limma::squeezeVar`) | Lₖ = log((yₖ + ε)/(bₖ + ε)) on target/reference spot counts |
| **Group comparison** | Welch *t* or Wilcoxon, Satterthwaite df | Welch *t* or Wilcoxon |
| **Null / margin** | TREAT against `lfc_threshold` (pilot plug-in τ_g) | TREAT centered on pilot baseline δ₀ = β₁, margin `tau` |
| **Multiplicity** | per-gene; BH (FDR) or Bonferroni | single test per (Target, Reference) pair |
| **Reported** | power & FDP over (K, `DE_lfc`) | power & empirical size over (K, effect) |

Both TREAT thresholds are estimated once from the pilot
(`estimate_lfc_threshold_from_pilot`, `estimate_lor_tau_from_pilot`; reachable via
`lfc_threshold = "auto"` / `tau = "auto"`) and shared across replicates, so type-I error
is controlled relative to a biologically meaningful minimum effect rather than a point null.

---

## Hyperparameter tuning

`tuneSpaDesignLambdas()` tunes two weights and returns a `spaCraftTuning` object:

- **`lambda_p`** — the pBANKSY recovery weight, tuned once by maximizing pairwise
  clustering reproducibility (ARI) across replicates.
- **`lambda_cond`** — the conditional-texture weight, tuned by maximizing a
  *regularized noncentrality* δ̂ = T̄ / (sd(T) + ε·s₀). This is scale-free and cannot be
  gamed by uniform shrinkage, unlike minimizing Var(T) alone. A null-calibration guard
  down-weights any `lambda_cond` whose empirical type-I error is inflated.

By default (`lambda_cond_selection = "per_effect"`) a pointwise-optimal `lambda_cond` is
returned for *every* effect size in `by_effect`, and `evaluatePowerSaLFC()` /
`evaluatePowerLOR()` match it per effect automatically when you pass `tuning_obj`.
Use `"shared"` for a single worst-case-optimal value.

---

## Function reference

| Stage | Functions |
|-------|-----------|
| **Object & accessors** | `createspaCraftObject`, `spaCraft-class`, `syntheticData()`, `syntheticData<-()`, `testingResult()`, `testingResult<-()`, `addTestingResult()` |
| **Feature / gene-set selection** | `featureSelection`, `featureSelectionStable`, `makeCustomGeneSets` |
| **Generative model (3 layers)** | `estimateCompositionParams`, `estimateGeometryParams`, `estimateExpressionParams` |
| **Simulate & recover domains** | `simulatespaCraft`, `pBANSKY`, `rearrangeSyntheticToPilot` |
| **Endpoints & TREAT thresholds** | `SaLFC`, `LOR`, `estimate_lfc_threshold_from_pilot`, `estimate_lor_tau_from_pilot` |
| **Power evaluation** | `evaluatePowerSaLFC`, `evaluatePowerLOR` |
| **Reporting** | `summary()`, `print()`, `plotPowerCurveSaLFC`, `plotPowerCurveLOR`, `as_powerSaLFC`, `as_powerLOR` |
| **Hyperparameter tuning** | `tuneSpaDesignLambdas` (+ `summary`/`plot` for `spaCraftTuning`) |
| **Visualization** | `plot(<spaCraft>, type = c("pilot", "synthetic_raw", "synthetic_mapped"))` |

`evaluatePowerSaLFC()` / `evaluatePowerLOR()` return a tagged data frame
(`powerSaLFC` / `powerLOR`), so `print()`, `summary()`, and `plotPowerCurve*()` work on
the returned object directly. If a `dplyr` pipeline strips the class, restore it with
`as_powerSaLFC()` / `as_powerLOR()`.

---

## Citation

If you use `spaCraft`, please cite the package (a methods manuscript is in
preparation):

```
Shin J, Chung D (2026). spaCraft: Analysis-Aware Power Calculation and
Sample-Size Planning for Multi-Sample Spatial Transcriptomics.
R package version 1.0.1.
```

```bibtex
@Manual{spaCraft,
  title  = {spaCraft: Analysis-Aware Power Calculation and Sample-Size Planning
            for Multi-Sample Spatial Transcriptomics},
  author = {Jungmin Shin and Dongjun Chung},
  year   = {2026},
  note   = {R package version 1.0.1}
}
```

---

## References

- Mukhopadhyay S., Li D., Dunson D.B. (2020). Fisher–Gaussian kernels.
  *JRSS-B* 82(5):1249–1271. [doi:10.1111/rssb.12390](https://doi.org/10.1111/rssb.12390)
- Saha A., Datta A. (2018). BRISC: nearest-neighbor Gaussian processes.
  *Stat* 7:e184. [doi:10.1002/sta4.184](https://doi.org/10.1002/sta4.184)
- Singhal V. *et al.* (2024). BANKSY. *Nature Genetics* 56:431–441.
  [doi:10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)
- Smyth G.K. (2004). Linear models and empirical Bayes methods. *SAGMB* 3(1):Article 3.
- McCarthy D.J., Smyth G.K. (2009). Testing relative to a fold-change threshold
  (TREAT). *Bioinformatics* 25(6):765–771. [doi:10.1093/bioinformatics/btp053](https://doi.org/10.1093/bioinformatics/btp053)
- Pya N., Wood S.N. (2015). Shape constrained additive models.
  *Statistics and Computing* 25(3):543–559.
- Maynard K.R. *et al.* (2021). DLPFC spatial transcriptomics.
  *Nature Neuroscience* 24:425–436. (data via `spatialLIBD`, Pardo et al. 2022)

---

## Authors & license

Jungmin Shin (aut, cre · `c16267@gmail.com`) · Dongjun Chung (aut)

Released under **GPL (≥ 3)**.
