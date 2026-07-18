<img src="spaCraft_logo.png" align="right" height="250" alt="spaCraft logo" />

# spaCraft

> **Analysis-aware power calculation and sample-size planning for multi-sample spatial transcriptomics.**

![version](https://img.shields.io/badge/version-1.0.1-5E8CB6)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-5E8CB6)
![license](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%203)-5E8CB6)
![platforms](https://img.shields.io/badge/platforms-Visium%20%7C%20Visium%20HD%20%7C%20Stereo--seq-5E8CB6)

`spaCraft` answers the first question of every comparative spatial study — **how many
samples per group do I need?** It learns a cohort-level generative model from a small
**pilot** dataset, then runs a full *generate → recover → test* Monte-Carlo pipeline to
report power per endpoint and per group size (*K*), summarized into a concrete
sample-size recommendation.

Unlike single-cell power tools, it models **domain composition**, **spatial geometry**,
and **spatially correlated expression** jointly, re-discovers spatial domains in every
synthetic replicate, and runs the *same* tests an analyst would actually run — so the
reported power reflects the real analysis pipeline, not an idealized oracle.

▶ **Explore it live (R Shiny):** [chunglab.bmi.osumc.edu/spaCraft](https://chunglab.bmi.osumc.edu/spaCraft/)

---

## Highlights

- **Pilot-based, single-pilot friendly** — learns the generator from real pilot
  sections; even one section per group works via automatic pseudo-replication
  (`generate_pseudo_reps`).
- **Three-layer spatial generative model** — composition (baseline-anchored
  Binomial–logit), geometry (Fisher–Gaussian kernel mixture, FGKMM), and expression
  (nearest-neighbor Gaussian process, NNGP).
- **Honest *generate → recover → test*** — domains are re-clustered in every replicate
  with a pilot-guided BANKSY extension (**pBANKSY**), so clustering error lives *inside*
  the power estimate rather than being assumed away.
- **Two endpoints, each with a TREAT minimum-effect margin** — a spatially-adjusted
  differential-expression test (**SaLFC**) and a baseline-anchored compositional test
  (**LOR**).
- **Validity-guarded tuning & calibrated reporting** — hyperparameters maximize a
  *regularized noncentrality* under a null-calibration guard (signal can't be bought by
  type-I inflation); power curves are SCAM-smoothed and monotone in *K*.
- **Cross-platform** — validated on **10x Visium**, **Visium HD**, and **Stereo-seq**.

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

**Requirements.** R ≥ 4.1.0 and a C++ toolchain (for `Rcpp`/`RcppArmadillo`). Heavy
steps are parallelized via `n_cores` (forked; falls back to serial on Windows).

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
obj <- featureSelection(object = obj, logfc_cutoff = 0.5, mean_in_cutoff = 0.8, max_num_gene = 3)
obj <- featureSelectionStable(object = obj, mean_global_cut = 0.1,min_domain_prop = 0.2,          
  max_abs_logfc_cut = 0.1, cv_cut = 0.2, max_num_gene = 10, use_low = FALSE)


## 3 — Labeling / null / spike-in gene sets ------------------------------------
sets <- makeCustomGeneSets(
  obj,
  G_svg_base       = obj@params_expression$top_genes,
  G_stable         = obj@params_expression$stable_genes,
  target_domain    = "WM",
  reference_domain = "Layer6",
  n_de             = 5L #number of spike genes
)

## 4 — Fit the three-layer generative model ------------------------------------
obj <- estimateCompositionParams(obj, target_domain = "WM", reference_domain = "Layer6")
obj <- estimateGeometryParams   (obj, n_cores = 4) #please check n_cores for your environment
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
  n_sim = 50, tuning_obj = tuning, n_cores = 4
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

## How it works

**Fit once, then sweep.** Learn the generator from pilot data, then map power over a
grid of group sample size (*K*) and effect size.

**1. Cohort-level generative model** (learned from the pilot):

| Layer | Captures | Model | Function |
|-------|----------|-------|----------|
| Composition | Per-domain abundance, case vs. control | Baseline-anchored **Binomial–logit** GLM (β₀, β₁) | `estimateCompositionParams()` |
| Geometry | Shape and placement of spatial domains | **Fisher–Gaussian kernel mixture (FGKMM)**, pooled across samples | `estimateGeometryParams()` |
| Expression | Domain means + spatial covariance + between-sample variance | **Gaussian process (NNGP)**, per-gene | `estimateExpressionParams()` |

**2. Generate → recover → test** (per Monte-Carlo replicate):

```
simulatespaCraft()  ─►  pBANSKY()  ─►  [rearrangeSyntheticToPilot()]  ─►  SaLFC() / LOR()
  synthetic cohort       recover         morph onto pilot geometry          endpoint test
  at effect size δ, K    domains (d̂)     (SaLFC only)                       + TREAT margin τ
```

Sweeping over the (`K`, effect) grid yields a rejection-rate surface; the **null row**
(effect = 0) returns empirical type-I error / FDP. `summary()` and `plotPowerCurve*()`
then fit a monotone SCAM (increasing in *K*) and report the **minimum *K*** reaching a
target power (default 0.8) for each effect size.

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

**Hyperparameter tuning.** `tuneSpaDesignLambdas()` returns a `spaCraftTuning` object
with two weights: **`lambda_p`** (pBANKSY recovery, tuned by pairwise clustering
reproducibility / ARI) and **`lambda_cond`** (conditional texture, tuned by a scale-free
*regularized noncentrality* δ̂ = T̄ / (sd(T) + ε·s₀) — it cannot be gamed by uniform
shrinkage, and a null-calibration guard down-weights any value with inflated type-I
error). With `lambda_cond_selection = "per_effect"` (default), `evaluatePowerSaLFC()` /
`evaluatePowerLOR()` match the per-effect optimum automatically when you pass
`tuning_obj`; use `"shared"` for a single worst-case-optimal value.

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

If you use `spaCraft`, please cite the package (a methods manuscript is in preparation):

```
Shin J, Xie J, and Chung D (2026). spaCraft: Analysis-Aware Power Calculation and
Sample-Size Planning for Multi-Sample Spatial Transcriptomics.
R package version 1.0.1.
```

```bibtex
@Manual{spaCraft,
  title  = {spaCraft: Analysis-Aware Power Calculation and Sample-Size Planning
            for Multi-Sample Spatial Transcriptomics},
  author = {Jungmin Shin, Juan Xie, and Dongjun Chung},
  year   = {2026},
  note   = {R package version 1.0.1}
}
```

---

## References

- Mukhopadhyay S., Li D., Dunson D.B. (2020). Fisher–Gaussian kernels. *JRSS-B*
  82(5):1249–1271. [doi:10.1111/rssb.12390](https://doi.org/10.1111/rssb.12390)
- Saha A., Datta A. (2018). BRISC: nearest-neighbor Gaussian processes. *Stat* 7:e184.
  [doi:10.1002/sta4.184](https://doi.org/10.1002/sta4.184)
- Singhal V. *et al.* (2024). BANKSY. *Nature Genetics* 56:431–441.
  [doi:10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)
- Smyth G.K. (2004). Linear models and empirical Bayes methods. *SAGMB* 3(1):Article 3.
- McCarthy D.J., Smyth G.K. (2009). Testing relative to a fold-change threshold (TREAT).
  *Bioinformatics* 25(6):765–771. [doi:10.1093/bioinformatics/btp053](https://doi.org/10.1093/bioinformatics/btp053)
- Pya N., Wood S.N. (2015). Shape constrained additive models. *Statistics and Computing*
  25(3):543–559.
- Maynard K.R. *et al.* (2021). DLPFC spatial transcriptomics. *Nature Neuroscience*
  24:425–436. (data via `spatialLIBD`, Pardo et al. 2022)

---

## Data downloads

All **21 pilot datasets** and their pre-computed results are hosted on the lab server and
downloadable over HTTPS under `https://chunglab.bmi.osumc.edu/spaCraft/data/…` (48 files,
**17.86 GB** total). Every file loads directly in R — for example:

```r
url  <- "https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/visium_human_brain_pilot_data_list.RData"
dest <- tempfile(fileext = ".RData"); options(timeout = 3600)
download.file(url, dest, mode = "wb"); load(dest)   # -> pilot_data_list
```

Power grids are plain CSV (`read.csv(url(...))`) and rendered curves are PDF. File-role
labels: `pilot data` = annotated pilot input · `pilot data (raw)` = raw (not-yet-analyzed)
spatial object · `fitted params` = estimated three-layer generator · `tuning oracle` =
`spaCraftTuning` object · `SaLFC/LOR grid (CSV)` = `evaluatePower*` tables · `SaLFC/LOR
curve (PDF)` = `plotPowerCurve*` figures.

### Visium

| Dataset | Species · Tissue · Condition | Files (click to download) · size |
|---|---|---|
| `human_brain` | human · brain (DLPFC) · case vs control · **example dataset** (PMID 33558695) | `pilot data` [visium_human_brain_pilot_data_list.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/visium_human_brain_pilot_data_list.RData) · **349.6 MB**<br>`fitted params` [visium_human_brain_parameter_est_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/visium_human_brain_parameter_est_obj.RData) · **349.6 MB**<br>`tuning oracle` [visium_human_brain_lambda_tuning_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/visium_human_brain_lambda_tuning_obj.RData) · **3.1 KB**<br>`SaLFC grid (CSV)` [Visium_human_brain_power_results_SaLFC.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/Visium_human_brain_power_results_SaLFC.csv) · **49.4 KB**<br>`LOR grid (CSV)` [Visium_human_brain_power_results_LOR.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/Visium_human_brain_power_results_LOR.csv) · **330.0 KB**<br>`SaLFC curve (PDF)` [Visium_human_barin_SaLFC_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/Visium_human_barin_SaLFC_power_curve.pdf) · **8.5 KB**<br>`LOR curve (PDF)` [Visium_human_barin_LOR_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain/Visium_human_barin_LOR_power_curve.pdf) · **8.2 KB** |
| `human_brain_MS` | human · brain (white matter) · MS lesion vs control (PMID 39501035) | `pilot data` [sskind_39501035_human_MS_pilot_annotated.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/sskind_39501035_human_MS_pilot_annotated.RData) · **579.6 MB**<br>`fitted params` [visium_human_brain_multiple_sclerosis_parameter_est_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/visium_human_brain_multiple_sclerosis_parameter_est_obj.RData) · **579.6 MB**<br>`tuning oracle` [visium_MS_human_lambda_tuning_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/visium_MS_human_lambda_tuning_obj.RData) · **2.6 KB**<br>`SaLFC grid (CSV)` [power_results_SaLFC.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/power_results_SaLFC.csv) · **49.0 KB**<br>`LOR grid (CSV)` [power_results_LOR.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/power_results_LOR.csv) · **148.8 KB**<br>`SaLFC curve (PDF)` [Visium_human_MS_SaLFC_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/Visium_human_MS_SaLFC_power_curve.pdf) · **8.6 KB**<br>`LOR curve (PDF)` [Visium_human_MS_LOR_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_MS/Visium_human_MS_LOR_power_curve.pdf) · **8.3 KB** |
| `human_brain_AD` | human · brain · Alzheimer's (PMID 36544231) | `pilot data (raw)` [sskind_36544231_human_AD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_AD/sskind_36544231_human_AD.RData) · **375.8 MB** |
| `human_brain_PD` | human · brain · Parkinson's (PMID 37667091) | `pilot data (raw)` [sskind_37667091_human_PD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_brain_PD/sskind_37667091_human_PD.RData) · **50.4 MB** |
| `human_breast` | human · breast · breast cancer (PMID 38114474) | `pilot data (raw)` [loki_3_human_breast_cancer.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_breast/loki_3_human_breast_cancer.RData) · **161.9 MB** |
| `human_heart_1` | human · heart · control vs case (PMID 35948637) | `pilot data (raw)` [loki_6_human_heart.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_heart_1/loki_6_human_heart.RData) · **772.8 MB** |
| `human_heart_2` | human · heart (Heart Cell Atlas) · benchmark split (PMID 37438528) | `pilot data (raw)` [loki_8_human_heart.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/human_heart_2/loki_8_human_heart.RData) · **794.9 MB** |
| `organoid_kidney` | human · kidney organoid · benchmark split (PMID 36776149) | `pilot data (raw)` [loki_2_organoid_kidney.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/organoid_kidney/loki_2_organoid_kidney.RData) · **9.9 MB** |
| `organoid_lung` | human · lung organoid · benchmark split (PMID 36776149) | `pilot data (raw)` [loki_2_organoid_lung.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/organoid_lung/loki_2_organoid_lung.RData) · **15.8 MB** |
| `5XFAD` | mouse · brain (cortex & hippocampus) · 5xFAD vs WT | `fitted params` [visium_5XFAD_parameter_est_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/visium_5XFAD_parameter_est_obj.RData) · **611.5 MB**<br>`tuning oracle` [visium_5XFAD_lambda_tuning_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/visium_5XFAD_lambda_tuning_obj.RData) · **2.6 KB**<br>`SaLFC grid (CSV)` [5XFAD_power_results_SaLFC.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/5XFAD_power_results_SaLFC.csv) · **83.7 KB**<br>`LOR grid (CSV)` [5XFAD_power_results_LOR.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/5XFAD_power_results_LOR.csv) · **334.7 KB**<br>`SaLFC curve (PDF)` [Visium_5XFAD_SaLFC_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/Visium_5XFAD_SaLFC_power_curve.pdf) · **8.1 KB**<br>`LOR curve (PDF)` [Visium_5XFAD_LOR_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/5XFAD/Visium_5XFAD_LOR_power_curve.pdf) · **8.3 KB** |
| `mouse_brain_AD` | mouse · brain · Alzheimer's (PMID 38036733) | `pilot data (raw)` [sskind_38036733_mouse_AD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain_AD/sskind_38036733_mouse_AD.RData) · **207.3 MB** |
| `mouse_brain_PD` | mouse · brain · Parkinson's (PMID 37667091) | `pilot data (raw)` [sskind_37667091_mouse_PD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain_PD/sskind_37667091_mouse_PD.RData) · **482.0 MB** |
| `mouse_brain_HD` | mouse · brain · Huntington's (PMID 40482637) | `pilot data (raw)` [sskind_40482637_mouse_HD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain_HD/sskind_40482637_mouse_HD.RData) · **486.1 MB** |
| `mouse_brain` | mouse · brain · benchmark split (PMID 36720873) | `pilot data (raw)` [loki_1_mouse_brain.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain/loki_1_mouse_brain.RData) · **578.6 MB** |
| `mouse_brain_2` | mouse · brain (hippocampus) · benchmark split (PMID 36776149) | `pilot data (raw)` [loki_2_mouse_brain.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain_2/loki_2_mouse_brain.RData) · **111.9 MB** |
| `mouse_brain_3` | mouse · brain · benchmark split (PMID 35027729) | `pilot data (raw)` [loki_5_mouse_brain.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_brain_3/loki_5_mouse_brain.RData) · **298.6 MB** |
| `mouse_bone` | mouse · bone marrow · benchmark split (PMID 36720873) | `pilot data (raw)` [loki_1_mouse_bone.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_bone/loki_1_mouse_bone.RData) · **71.5 MB** |
| `mouse_skeleton_muscle` | mouse · skeletal muscle · ALS (PMID 39932195) | `pilot data (raw)` [sskind_39932195_mouse_ALS.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium/mouse_skeleton_muscle/sskind_39932195_mouse_ALS.RData) · **58.9 MB** |

### Stereo-seq

| Dataset | Species · Tissue · Condition | Files (click to download) · size |
|---|---|---|
| `mouse_embroyo` | mouse · E12.5 whole embryo · developmental · **analyzed** | `fitted params` [Stereo_seq_mouse_parameter_est_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_mouse_parameter_est_obj.RData) · **3.51 GB**<br>`tuning oracle` [Stereo_seq_mouse_lambda_tuning_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_mouse_lambda_tuning_obj.RData) · **2.0 KB**<br>`SaLFC grid (CSV)` [Stereo_seq_embryo_power_results_SaLFC.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_embryo_power_results_SaLFC.csv) · **48.8 KB**<br>`LOR grid (CSV)` [Stereo_seq_embryo_power_results_LOR.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_embryo_power_results_LOR.csv) · **89.8 KB**<br>`SaLFC curve (PDF)` [Stereo_seq_SaLFC_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_SaLFC_power_curve.pdf) · **8.3 KB**<br>`LOR curve (PDF)` [Stereo_seq_LOR_power_curve.pdf](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_embroyo/Stereo_seq_LOR_power_curve.pdf) · **8.2 KB** |
| `mouse_brain_AD` | mouse · brain · Alzheimer's (PMID 38819990) | `pilot data (raw)` [sskind_38819990_mouse_AD.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Stereo_seq/mouse_brain_AD/sskind_38819990_mouse_AD.RData) · **1.13 GB** |

### Visium HD

| Dataset | Species · Tissue · Condition | Files (click to download) · size |
|---|---|---|
| `Visium_HD` | human · colorectal (CRC) · tumor vs normal · **analyzed** | `fitted params` [visium_HD_parameter_est_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/visium_HD_parameter_est_obj.RData) · **2.14 GB**<br>`tuning oracle` [Visium_HD_lambda_tuning_obj.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/Visium_HD_lambda_tuning_obj.RData) · **3.0 KB**<br>`spaCraft object` [spaDesign2_obj_visium_hd.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/spaDesign2_obj_visium_hd.RData) · **2.14 GB**<br>`fitted object` [spaDesign2_obj_after_fitting_visium_hd.RData](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/spaDesign2_obj_after_fitting_visium_hd.RData) · **2.14 GB**<br>`SaLFC grid (CSV)` [power_results.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/power_results.csv) · **33.0 KB**<br>`LOR grid (CSV)` [lor_results.csv](https://chunglab.bmi.osumc.edu/spaCraft/data/Visium_HD/lor_results.csv) · **67.0 KB** |

**Notes.**
- The `human_brain` example file is also reachable at the short URL used in **Quick start**
  (`…/spaCraft/visium_human_brain_pilot_data_list.RData`).
- Some legacy filenames keep a `barin` / `spaDesign2` spelling — download URLs must match
  the exact filename shown above.
- The Stereo-seq `mouse_embroyo` folder also holds a stray copy of
  `sskind_38819990_mouse_AD.RData`, identical to the `mouse_brain_AD` file (omitted above).
- Sizes are exact on-disk byte counts (1 KB = 1024 B); large `parameter_est` / `spaDesign2`
  objects are multi-GB, so the first download may take a while.

<details>
<summary><i>Server maintainer — how these files are published</i></summary>

Shiny Server serves the app's <code>www/</code> folder as static files, so the catalog is
mirrored there (structure preserved to avoid filename clashes):

```bash
sudo mkdir -p /srv/shiny-server/spaCraft/www/data
sudo cp -a /srv/shiny-server/spaCraft/catalog/. /srv/shiny-server/spaCraft/www/data/
sudo find /srv/shiny-server/spaCraft/www/data -name .DS_Store -delete
sudo chown -R shiny:shiny /srv/shiny-server/spaCraft/www
sudo chmod -R a+rX /srv/shiny-server/spaCraft/www
sudo systemctl restart shiny-server.service
# verify one file (expect HTTP 200):
curl -sSI http://127.0.0.1:3838/spaCraft/data/Visium/human_brain/visium_human_brain_pilot_data_list.RData | head -3
```
</details>

---

## Authors & license

Jungmin Shin (aut, cre · `c16267@gmail.com`) · Dongjun Chung (aut) — released under **GPL (≥ 3)**.
