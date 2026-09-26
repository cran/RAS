# RAS: Regional Association Score for GWAS

The **RAS** package implements the Regional Association Score method for genome-wide association studies. For every SNP it measures the strength of association in the surrounding region, arranges these regional scores along the chromosome into a −log₁₀(*p*) profile, and locates association regions on that profile. The method supports continuous and binary traits.

If you use this package in your research, please cite:

> Y. Jiang & H. Zhang, Empowering genome-wide association studies via a visualizable test based on the regional association score, *Proc. Natl. Acad. Sci. U.S.A.* 122(9) e2419721122 (2025). https://doi.org/10.1073/pnas.2419721122

---

## What is new in 1.1.0

* **One name, one engine.** `ras()`, `ras_scan()` and `ras_detect()` are now the compiled, disk-backed implementations; the pure-R implementations of 1.0.x are available as `ras_original()`, `ras_scan_original()` and `ras_detect_original()`. `ras()` still accepts an in-memory genotype matrix.
* **Two region detectors.** `detector = "box"` (default) is a box-scan detector that reports intervals `[tau_L, tau_R]`, delimits broad plateau-shaped regions, and calibrates its own threshold on permuted phenotypes. `detector = "changepoint"` is the published two-pass changepoint detector of 1.0.x. See *Choosing a detector* below.
* **A disk-backed engine.** `ras()` streams genotypes from a chunked on-disk file (`.rasbin`), so peak memory is set by a chunk size rather than by the chromosome. Profiles agree with `ras_original()` to machine precision.
* **External summary statistics as weights.** `ras_harmonize_sumstats()` aligns an independent GWAS to your genotype map; `ras_scan_external()` then scans all your samples with those fixed weights, without a train/holdout split.

---

## Installation

```r
install.packages("RAS")                       # CRAN release

# development version
install.packages("remotes")
remotes::install_github("hepingzhangyale/RAS")
```

Installing from source needs a C compiler ([Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows).

---

## Quick start

```r
library(RAS)

result <- ras(
  geno, phenotype, covariates,
  covariate_cols = c("age", "sex", paste0("pc", 1:10)),
  is_continuous  = TRUE,
  chrom          = 1,
  save_dir       = "results/"
)

print(result)              # detected regions
plot(result)               # full-chromosome scan profile
plot(result, zoom = TRUE)  # zoomed view around each region
result$detection$regions   # one row per interval: pos_L, pos_R, T_box, anchor_pos, ...
```

`geno` may be an in-memory matrix (converted for the run) or, for anything large, a `.rasbin` file. By default the box-scan detector calibrates its threshold by repeating the scan on `box_null = 100` permuted phenotypes (exact order statistic); `box_method = "gumbel"` gets a stable threshold from about 20 permutations at the price of being slightly liberal. For a large cohort run the null scans on a random subset of individuals (`box_null_n = 20000`; the null maximum depends on the linkage structure, not on the sample size), or calibrate once with `ras_box_calibrate()` and pass `box_calibration` (see below).

### Large cohorts: convert once, then run from the file

```r
bed_to_rasbin("cohort_chr1.bed", "cohort_chr1.rasbin")   # from PLINK 1 .bed/.bim/.fam
# or: geno_to_rasbin(geno_matrix, "cohort_chr1.rasbin") # from an in-memory matrix / .rds

result <- ras(
  "cohort_chr1.rasbin", phenotype, covariates,
  covariate_cols = c("age", "sex", paste0("pc", 1:10)),
  is_continuous  = TRUE,
  chrom          = 1,
  chunk_snps     = 5000,     # SNP columns held in memory at a time
  cores          = 4,        # repetitions run in parallel
  save_dir       = "results/"
)
```

Sample order in `phenotype` and `covariates` must match the `.rasbin` file (the `.fam` order for `bed_to_rasbin()`).

### Binary traits

For binary traits `ras()` uses a Rao score test that fits the covariate-only null model once and evaluates each window in closed form. The per-window logistic-regression Wald scan that was the default in 1.0.x is still available as `scan_test = "glm"`; it runs the in-memory `ras_original()` route (so it needs a genotype matrix) and gives the same detected regions.

---

## Choosing a detector

Both detectors work on the same scan profile; the choice only affects Stages 2 and 3.

| | `detector = "box"` (default) | `detector = "changepoint"` |
|---|---|---|
| What it reports | intervals `[tau_L, tau_R]` with an anchor | single changepoint positions |
| Method | window mean above both flanks and above the profile background, in raw −log₁₀(*p*) units (`ras_box_detect()`) | segmented regression + Davies test in a sliding window, then local re-validation (`ras_detect()`, `ras_validate()`) |
| Broad plateau-shaped regions | delimited at their edges | may be missed or placed inside the plateau |
| Threshold | calibrated on permuted phenotypes (automatic, or `ras_box_calibrate()`) | fixed p-value defaults (published) |
| Backwards compatibility | new in 1.1.0 | same results as RAS 1.0.x |

```r
## 1. null profiles: the same scan on permuted phenotypes
nulls <- lapply(1:200, function(b) {
  set.seed(b)
  ras_scan("cohort_chr1.rasbin", sample(phenotype), covariates,
                covariate_cols = covs, is_continuous = TRUE,
                save_dir = tempfile())$y
})

## 2. threshold at family-wise error 0.05
cal <- ras_box_calibrate(nulls, alpha = 0.05)

## 3. run with that calibration (skips the automatic one)
result <- ras("cohort_chr1.rasbin", phenotype, covariates,
              covariate_cols = covs, is_continuous = TRUE,
              box_calibration = cal, chrom = 1, save_dir = "results/")

## the changepoint detector of 1.0.x
result_cp <- ras("cohort_chr1.rasbin", phenotype, covariates,
                 covariate_cols = covs, is_continuous = TRUE,
                 detector = "changepoint", chrom = 1, save_dir = "results/")
```

The scan profile is never affected by the detector: `result$scan` is the same either way, so the two detectors can be compared on one scan.

---

## External summary statistics as weights

When an independent GWAS of the same trait is available, its effect sizes can replace the within-sample training split. The scan then uses every individual of the target cohort once.

```r
h <- ras_harmonize_sumstats(
  sumstats = "external_gwas.tsv",           # rsID / chr / pos / A1 / A2 / beta / p
  map      = "cohort_chr1.bim"              # target variant map
)
h$qc                                        # one row per harmonisation stage

scan <- ras_scan_external(
  "cohort_chr1.rasbin", phenotype, covariates,
  covariate_cols = covs, weights = h$weights, is_continuous = TRUE
)
det <- ras_box_detect(scan$x, scan$y, calibration = cal)   # or ras_detect() + ras_validate()
```

`ras_sumstats_report()` gives a non-fatal preview of how a summary statistics file lines up with a map (ID convention, genome build, coverage along the chromosome) before anything is committed.

---

## Step-by-step (advanced)

```r
# Step 1: averaged -log10(p) profile
scan <- ras_scan(geno, phenotype, covariates,
                 covariate_cols = c("age", "sex", paste0("pc", 1:10)),
                 is_continuous = TRUE, chrom = 1, save_dir = "results/")

# Step 2 + 3, changepoint detector
detected <- ras_detect(scan$x, scan$y, window_size = 3000,
                       slope.p.values.threshold.left  = 1e-10,
                       slope.p.values.threshold.right = 1e-20)
final <- ras_validate(detected, x = scan$x, y = scan$y,
                      this.skip = 10, p.value.threshold = 1e-10)

# Step 2 + 3, box-scan detector
final_box <- ras_box_detect(scan$x, scan$y, calibration = cal)

# Step 4: plot either one
result <- structure(list(scan = scan, detection = final, chrom = 1, save_dir = "results/"),
                    class = "ras")
plot(result)
```

## Documentation

```r
?RAS                     # package overview and pipeline description
?ras                     # one-call entry point (matrix or .rasbin file)
?ras_original            # the pure-R in-memory implementation of 1.0.x
?ras_box_detect          # box-scan region detector
?ras_box_calibrate       # its threshold
?ras_detect              # changepoint detector, first pass
?ras_validate            # changepoint detector, second pass
?ras_harmonize_sumstats  # external summary statistics
?plot.ras                # plotting
?ras_memory              # memory and CPU diagnostics for the in-memory route
```

## Authors

- Jiahe Jin &lt;jiahe.jin@yale.edu&gt;
- Yiran Jiang &lt;yiran.jiang@uky.edu&gt;
- Heping Zhang &lt;heping.zhang@yale.edu&gt; (maintainer)
