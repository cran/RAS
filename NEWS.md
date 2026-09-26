# RAS 1.1.2

This release adds a second region detector, a disk-backed engine in compiled
code for large cohorts, and support for external GWAS summary statistics as
regional weights. The compiled engine is now the implementation behind the
familiar names; the 1.0.3 implementations remain available under `_original`
names.

## Function names

* The compiled, disk-backed implementations are now the default and carry the
  plain names: `ras()`, `ras_scan()`, `ras_detect()`,
  `compute_gwas_weights()`, `screen_forward_max_region()`,
  `ras_scan_external()`. The pure-R, in-memory implementations that 1.0.x
  shipped under those names are available, unchanged, as `ras_original()`,
  `ras_scan_original()`, `ras_detect_original()`,
  `compute_gwas_weights_original()`, `screen_forward_max_region_original()`
  and `ras_scan_external_original()`.
* `ras()` and `ras_scan()` keep the 1.0.x calling convention: the first
  argument is still `geno` and still accepts an in-memory genotype matrix
  (converted to a temporary `.rasbin` file for the run); it now also accepts
  the path to a `.rasbin` file. For continuous traits the profile is identical
  to 1.0.x to machine precision. For binary traits the default scan is now the
  Rao score test (`scan_test = "score"`); the 1.0.x per-window logistic
  regression is still available as `scan_test = "glm"`, which runs
  `ras_original()` and needs a genotype matrix.
* The `_fast` names of the development versions (`ras_fast()`,
  `ras_scan_fast()`, ...) are kept as aliases.

## Two region detectors

* `ras()` gains `detector = c("box", "changepoint")`, default `"box"`, the
  new box-scan region detector; `"changepoint"` is the original two-pass
  detector, `ras_detect()` followed by `ras_validate()`, which reproduces the
  published method and the results of 1.0.3. `ras_original()` keeps
  `"changepoint"` as its default.
* `ras()` calibrates the box-scan threshold on the data: unless
  `box_threshold` or `box_calibration` is given, the scan is repeated on
  `box_null = 100` permuted phenotypes (Freedman-Lane residual permutation
  for continuous traits) and `ras_box_calibrate()` turns their maxima into a
  threshold at family-wise level `box_alpha = 0.05` (exact order statistic,
  or a Gumbel fit with `box_method = "gumbel"` when only 20 to 50
  permutations are affordable). `box_null_n` runs the null scans on a random
  subset of individuals, which makes the calibration affordable at biobank
  scale. The calibration is returned as `$box_calibration`.
* `ras_scan()` gains `rows`, the subset of samples within which the
  train/hold-out splits are drawn.
* New `ras_box_detect()`: scores every window of a set of widths by how far
  its mean stands above both adjacent flanks and above the profile
  background, in the raw -log10(p) units of the profile, and reports the
  intervals [tau_L, tau_R] that clear a threshold. Because the window must
  stand above each flank, a window inside a plateau or on the shoulder of a
  peak scores about zero, so broad plateau-shaped regions are delimited at
  their edges and neighbouring peaks stay separate. No noise scale is
  estimated. Regions at the ends of a profile are reported (`edge = "open"`).
* New `ras_box_calibrate()`: turns null scan profiles (the same scan on
  permuted phenotypes) into the threshold `ras_box_detect()` needs, as an
  order statistic of the null maxima, with an optional held-out error rate and
  a Gumbel extrapolation to longer profiles. New `ras_box_stat()` exposes the
  statistic itself.
* `print()` and `plot()` for `"ras"` objects understand both detectors: for a
  box-scan result the shaded areas are the detected intervals and the overlay
  plot shows the box-scan score in place of the Davies test significance.

## Disk-backed C engine

* New `.rasbin` genotype format: a flat column-major binary file from which an
  arbitrary SNP range is one seek and one read. `geno_to_rasbin()` converts an
  in-memory matrix or an `.rds`/`.csv`/`.tsv` file; `bed_to_rasbin()` converts
  a PLINK 1 `.bed/.bim/.fam` fileset in fixed-size SNP chunks without ever
  holding the dense matrix (int8 storage, 1 byte per genotype);
  `rasbin_header()` and `rasbin_read_chunk()` read it back.
* `ras_scan()` and `ras()`: the Stage-1 scan streams SNP chunks from a
  `.rasbin` file through C (`compute_gwas_weights()`,
  `screen_forward_max_region()`), so peak memory is bounded by `chunk_snps`
  rather than by the chromosome. The statistics are identical to
  `ras_scan_original()` (exact FWL for continuous traits, Rao score test for
  binary traits), and the profiles agree to machine precision.
* `ras_detect()`: a C port of the first-pass changepoint detector (Muggeo's
  segmented fit with bootstrap restarts and the Davies test). Resampling draws
  come from R's RNG stream, so results are statistically, not bit-for-bit,
  equal to `ras_detect_original()`.
* `ras_scan()`/`ras()` take `cores` (repetitions dispatched to a
  PSOCK cluster) and `keep_reps` (return the per-repetition profiles as well
  as their mean).
* Binary-trait GWAS weights on the fast route handle SNPs with missing
  dosages by downdating shared covariate cross-products in blocks instead of
  refitting each such SNP on its own, a 20-fold saving on post-QC sequencing
  data where almost every SNP carries at least one missing value.

## External summary statistics as weights

* New `ras_harmonize_sumstats()` (and the building blocks
  `ras_weights_from_sumstats()`, `ras_sumstats_report()`,
  `ras_read_variant_dictionary()`, `ras_map_ids_to_rsid()`,
  `ras_liftover_sumstats()`): match an external GWAS to the target genotype
  map by rsID or by position and allele pair, check the genome build, align
  effect alleles, drop strand-ambiguous variants, and return aligned weights
  together with a one-row-per-stage QC table.
* New `ras_scan_external()` (and the in-memory `ras_scan_external_original()`):
  the Stage-1 scan with fixed external weights over every individual, with no
  train/holdout split and no repetitions.

## Fixes

* The package now builds on Linux systems whose C library is not glibc
  (for example Alpine Linux / musl): `malloc_trim()` is a glibc extension and
  is now compiled only where glibc is present; `release_memory()` returns `NA`
  elsewhere, as it already did on Windows and macOS; the unit tests check
  0/1 or `NA` according to what the build compiled in. (1.1.0 and 1.1.1 were submitted and
  withdrawn for this reason; neither was released.)
* `ras_validate()` no longer stops on a non-finite Davies p-value from a
  degenerate window, and no longer keeps a candidate whose profile value is
  `NA`.

# RAS 1.0.3

Resubmission addressing the CRAN reviewer's comments on the 1.0.0 submission,
with a bug fix folded in.

* Removed the `\dontrun{}` examples that called the unexported internal
  functions `plot_ras_scan()` and `plot_ras_zoom_regions()` through `:::`.
  These functions are internal (`@keywords internal`) and are already exercised
  by the runnable example of the exported `plot()` method for `"ras"` objects.
* Made the `ras_memory()` "abort" example runnable by wrapping the intentional
  `stop()` in `try()`, instead of hiding it in `\dontrun{}`.
* The plotting functions now capture the caller's graphics parameters and
  restore them via an immediate `on.exit(par(oldpar))`, so plotting a `"ras"`
  object no longer leaves the user's `par()` settings modified.
* Fixed a crash in `screen_forward_max_region()` (and therefore `ras()`) when
  covariates contained missing values. Incomplete cases are now dropped before
  the model matrix is built, matching the documented "removes incomplete cases"
  behaviour, on all three scan paths (continuous exact-FWL, binary score, and
  the legacy glm path).
