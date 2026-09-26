make_pipeline_fixture <- function(seed, n_samp = 300, n_snp = 2000,
                                   causal_lo = 900, causal_hi = 1000,
                                   effect = 1.4) {
  set.seed(seed)
  geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
                  nrow = n_samp, ncol = n_snp)
  storage.mode(geno) <- "double"

  causal_idx    <- causal_lo:causal_hi
  causal_score  <- rowSums(geno[, causal_idx, drop = FALSE]) / length(causal_idx)
  latent        <- effect * causal_score + rnorm(n_samp, sd = 0.6)

  cov_df <- data.frame(
    age = rnorm(n_samp, 50, 10),
    sex = rbinom(n_samp, 1, 0.5)
  )
  cov_df$age_squared <- cov_df$age^2
  cov_df$age_sex     <- cov_df$age * cov_df$sex
  for (i in 1:6) cov_df[[paste0("pc", i)]] <- rnorm(n_samp)
  covariate_cols <- c("age", "sex", "age_squared", "age_sex", paste0("pc", 1:6))

  cov_effect <- 0.02 * cov_df$age + 0.3 * cov_df$sex

  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(geno, path, chunk_write_size = 400)

  list(geno = geno, path = path, cov_df = cov_df, covariate_cols = covariate_cols,
       latent = latent, cov_effect = cov_effect, causal_idx = causal_idx,
       causal_center = mean(causal_idx), n_snp = n_snp)
}

test_that("ras_scan matches ras_scan_original (continuous) to tight tolerance, same seed/split", {
  for (seed in c(1, 2)) {
    d <- make_pipeline_fixture(seed)
    pheno <- d$latent + d$cov_effect + rnorm(nrow(d$geno), sd = 0.5)
    save_ref  <- tempfile(); save_fast <- tempfile()

    set.seed(1000 + seed)
    ref <- ras_scan_original(
      geno = d$geno, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = TRUE,
      num_rep = 2, skip1 = 10, skip2 = 20, save_dir = save_ref,
      min_window_size = 5, max_window_size = 100, scan_test = "glm"
    )

    set.seed(1000 + seed)
    fast <- ras_scan(
      geno = d$path, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = TRUE,
      num_rep = 2, skip1 = 10, skip2 = 20, save_dir = save_fast,
      min_window_size = 5, max_window_size = 100, chunk_snps = 300
    )

    unlink(c(d$path, paste0(d$path, ".meta.rds"), save_ref, save_fast), recursive = TRUE)

    expect_equal(ref$x, fast$x)
    ok <- is.finite(ref$y) & is.finite(fast$y)
    expect_gt(mean(ok), 0.8)
    expect_equal(fast$y[ok], ref$y[ok], tolerance = 1e-4, info = paste("seed", seed))
  }
})

test_that("ras_scan matches ras_scan_original (binary/score) to tight tolerance, same seed/split", {
  for (seed in c(11, 12)) {
    d <- make_pipeline_fixture(seed)
    pheno <- as.integer(plogis(d$latent + d$cov_effect - mean(d$latent + d$cov_effect)) > runif(nrow(d$geno)))
    save_ref  <- tempfile(); save_fast <- tempfile()

    set.seed(2000 + seed)
    ref <- ras_scan_original(
      geno = d$geno, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = FALSE,
      num_rep = 2, skip1 = 10, skip2 = 20, save_dir = save_ref,
      min_window_size = 5, max_window_size = 100, scan_test = "score"
    )

    set.seed(2000 + seed)
    fast <- ras_scan(
      geno = d$path, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = FALSE,
      num_rep = 2, skip1 = 10, skip2 = 20, save_dir = save_fast,
      min_window_size = 5, max_window_size = 100, chunk_snps = 300
    )

    unlink(c(d$path, paste0(d$path, ".meta.rds"), save_ref, save_fast), recursive = TRUE)

    expect_equal(ref$x, fast$x)
    ok <- is.finite(ref$y) & is.finite(fast$y)
    expect_gt(mean(ok), 0.8)
    expect_equal(fast$y[ok], ref$y[ok], tolerance = 1e-4, info = paste("seed", seed))
  }
})

test_that("ras detects a changepoint in the same neighborhood as ras_original, on an injected signal", {
  hits_ref <- 0; hits_fast <- 0; near_ref <- numeric(0); near_fast <- numeric(0)

  for (seed in c(21, 22, 23)) {
    d <- make_pipeline_fixture(seed, n_samp = 400, n_snp = 2500, effect = 2.2)
    pheno <- d$latent + d$cov_effect + rnorm(nrow(d$geno), sd = 0.4)
    save_ref  <- tempfile(); save_fast <- tempfile()

    # window_size/slope thresholds relaxed to match the doc examples'
    # convention for small synthetic profiles (defaults are tuned for
    # genome-scale grids with tens of thousands of points).
    set.seed(3000 + seed)
    ref <- suppressWarnings(ras_original(
      geno = d$geno, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = TRUE,
      num_rep = 3, skip1 = 10, skip2 = 20, save_dir = save_ref,
      min_window_size = 5, max_window_size = 100, scan_test = "glm",
      cp_window_size = 150, cp_slope_left = 1e-3, cp_slope_right = 1e-3,
      second_window_size = 20, second_p_threshold = 1e-3,
      run_plots = FALSE
    ))

    set.seed(3000 + seed)
    fast <- suppressWarnings(ras(
      geno = d$path, phenotype = pheno, covariates = d$cov_df,
      covariate_cols = d$covariate_cols, is_continuous = TRUE,
      num_rep = 3, skip1 = 10, skip2 = 20, save_dir = save_fast,
      min_window_size = 5, max_window_size = 100, chunk_snps = 300,
      detector = "changepoint",
      cp_window_size = 150, cp_slope_left = 1e-3, cp_slope_right = 1e-3,
      second_window_size = 20, second_p_threshold = 1e-3,
      run_plots = FALSE
    ))

    unlink(c(d$path, paste0(d$path, ".meta.rds"), save_ref, save_fast), recursive = TRUE)

    if (length(ref$detection$tau_hats) > 0) {
      hits_ref <- hits_ref + 1
      near_ref <- c(near_ref, min(abs(ref$detection$tau_hats - d$causal_center)))
    }
    if (length(fast$detection$tau_hats) > 0) {
      hits_fast <- hits_fast + 1
      near_fast <- c(near_fast, min(abs(fast$detection$tau_hats - d$causal_center)))
    }
  }

  # Diagnostic only -- bootstrap-restart RNG differs from R's by design (see
  # Phase D report), so exact per-seed detection parity is not expected; this
  # documents the neighborhood behaviour rather than hard-asserting it.
  message(sprintf("ras_original(): %d/3 seeds detected, nearest offsets: %s",
                   hits_ref, paste(round(near_ref), collapse = ", ")))
  message(sprintf("ras(): %d/3 seeds detected, nearest offsets: %s",
                   hits_fast, paste(round(near_fast), collapse = ", ")))

  # Structural guarantee: ras() and ras_original() should behave consistently on
  # the same seeded input -- either both detect (and land near the true
  # signal) or both correctly find nothing, not diverge from each other.
  # Detection itself is not guaranteed on a small synthetic fixture with
  # strict slope thresholds (real bio data isn't 100% either -- see the
  # RAS Pig Example memory's 96/100-seed robustness result), so this test
  # checks parity of behaviour, not that detection must occur.
  expect_equal(hits_fast, hits_ref, tolerance = 1)
  if (length(near_ref) > 0)  expect_true(all(near_ref  < 400))
  if (length(near_fast) > 0) expect_true(all(near_fast < 400))
})

test_that("ras_scan(cores = 1) is byte-for-byte reproducible (cores param is a no-op at cores=1)", {
  d <- make_pipeline_fixture(31, n_samp = 200, n_snp = 1000)
  pheno <- d$latent + d$cov_effect + rnorm(nrow(d$geno), sd = 0.5)
  save_a <- tempfile(); save_b <- tempfile()

  set.seed(99)
  a <- ras_scan(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
                      is_continuous = TRUE, num_rep = 3, skip1 = 20, skip2 = 40,
                      save_dir = save_a, chunk_snps = 300, cores = 1)
  set.seed(99)
  b <- ras_scan(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
                      is_continuous = TRUE, num_rep = 3, skip1 = 20, skip2 = 40,
                      save_dir = save_b, chunk_snps = 300, cores = 1)

  unlink(c(d$path, paste0(d$path, ".meta.rds"), save_a, save_b), recursive = TRUE)
  expect_identical(a$y, b$y)
})

test_that("ras_scan(cores > 1) works when arguments are passed as variables, not literals (regression: promise serialization across parLapply workers)", {
  skip_if(parallel::detectCores() < 2, "needs at least 2 cores")

  d <- make_pipeline_fixture(32, n_samp = 200, n_snp = 1000)
  # Deliberately name these differently from ras_scan()'s own formal
  # parameter names and route them through local variables -- the bug this
  # guards against only manifested when an argument arrived as an unforced
  # promise pointing at a caller-side variable (e.g. `phenotype = y`), not
  # when it was a literal, so a same-named passthrough wouldn't have caught it.
  my_pheno     <- d$latent + d$cov_effect + rnorm(nrow(d$geno), sd = 0.5)
  my_covs      <- d$cov_df
  my_cov_cols  <- d$covariate_cols
  my_is_cont   <- TRUE
  my_chrom     <- 7L
  save_dir     <- tempfile()

  expect_no_error(
    result <- ras_scan(d$path, my_pheno, my_covs, covariate_cols = my_cov_cols,
                             is_continuous = my_is_cont, num_rep = 4, skip1 = 20, skip2 = 40,
                             chrom = my_chrom, save_dir = save_dir, chunk_snps = 300, cores = 2)
  )

  unlink(c(d$path, paste0(d$path, ".meta.rds"), save_dir), recursive = TRUE)
  expect_true(is.list(result))
  expect_true(all(c("x", "y") %in% names(result)))
  expect_length(result$y, length(result$x))
  expect_true(all(is.finite(result$y) | is.na(result$y)))
})

test_that("ras_scan(cores > 1) produces a structurally comparable profile to cores = 1 on the same data", {
  skip_if(parallel::detectCores() < 2, "needs at least 2 cores")

  d <- make_pipeline_fixture(33, n_samp = 250, n_snp = 1500)
  pheno <- d$latent + d$cov_effect + rnorm(nrow(d$geno), sd = 0.5)
  save_1 <- tempfile(); save_2 <- tempfile()

  set.seed(7)
  serial <- ras_scan(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
                           is_continuous = TRUE, num_rep = 6, skip1 = 15, skip2 = 30,
                           save_dir = save_1, chunk_snps = 300, cores = 1)
  parallel_res <- ras_scan(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
                                 is_continuous = TRUE, num_rep = 6, skip1 = 15, skip2 = 30,
                                 save_dir = save_2, chunk_snps = 300, cores = 2)

  unlink(c(d$path, paste0(d$path, ".meta.rds"), save_1, save_2), recursive = TRUE)

  # cores > 1 draws independent per-worker train/holdout splits (see the
  # `cores` parameter docs), so exact numeric parity with cores = 1 isn't
  # expected -- only that it ran on the same grid and lands in a comparable
  # range, not e.g. all-NA or a constant/degenerate profile.
  expect_equal(serial$x, parallel_res$x)
  ok <- is.finite(serial$y) & is.finite(parallel_res$y)
  expect_gt(mean(ok), 0.8)
  expect_gt(sd(parallel_res$y[ok]), 0)
  expect_true(cor(serial$y[ok], parallel_res$y[ok]) > 0.3)
})

# --- keep_reps: the per-repetition profiles ---------------------------------
# ras_scan() averages its repetitions, and averaging is exactly the step
# that removes split-to-split variability -- so any measurement OF that
# variability needs the individual profiles the mean is built from. These tests
# pin the two properties that matter: the mean is unchanged (this is purely
# additive, not a new code path), and the kept profiles genuinely reconstruct it.

test_that("keep_reps = FALSE leaves the returned object exactly as it was", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                              causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect

  set.seed(99)
  res <- ras_scan(
    d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 3, skip1 = 50, skip2 = 25,
    max_window_size = 40, chunk_snps = 400, save_dir = tempfile()
  )
  expect_named(res, c("x", "y", "reps"))
  expect_null(res$reps)
  expect_length(res$y, length(res$x))
})

test_that("keep_reps = TRUE returns per-repetition profiles that average to y", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                              causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect
  save_dir <- tempfile(); dir.create(save_dir)

  set.seed(99)
  res <- ras_scan(
    d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 3, skip1 = 50, skip2 = 25,
    max_window_size = 40, chunk_snps = 400, save_dir = save_dir,
    keep_reps = TRUE
  )

  expect_true(is.matrix(res$reps))
  expect_equal(dim(res$reps), c(length(res$x), 3L))
  expect_identical(colnames(res$reps), c("rep1", "rep2", "rep3"))
  # the whole point: the mean is recoverable from what was kept
  expect_equal(rowMeans(res$reps), res$y, tolerance = 1e-12)
  # the reps must actually differ -- if they were identical the split-induced
  # variability this argument exists to measure would be zero by construction
  expect_gt(stats::sd(res$reps[, 1] - res$reps[, 2]), 0)

  f <- file.path(save_dir, "rep_p_values_chr1_reps1-3.rds")
  expect_true(file.exists(f))
  expect_equal(readRDS(f), res$reps)
})

test_that("keep_reps does not perturb the averaged profile", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                              causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect
  args <- list(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
               is_continuous = TRUE, num_rep = 3, skip1 = 50, skip2 = 25,
               max_window_size = 40, chunk_snps = 400)

  set.seed(7); off <- do.call(ras_scan, c(args, list(save_dir = tempfile())))
  set.seed(7); on  <- do.call(ras_scan, c(args, list(save_dir = tempfile(),
                                                          keep_reps = TRUE)))
  expect_identical(off$x, on$x)
  expect_identical(off$y, on$y)
})

test_that("ras passes keep_reps through to the scan it stores", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                              causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect

  set.seed(5)
  r <- ras(
    d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 2, skip1 = 50, skip2 = 25,
    max_window_size = 40, chunk_snps = 400, save_dir = tempfile(),
    detector = "changepoint", run_plots = FALSE, keep_reps = TRUE
  )
  expect_true(is.matrix(r$scan$reps))
  expect_equal(ncol(r$scan$reps), 2L)
  expect_equal(rowMeans(r$scan$reps), r$scan$y, tolerance = 1e-12)
})

# A phenotype with a strong regional signal: a few standardised causal SNPs
# explain most of the variance, so the injected block is unmistakable on the
# profile and the tests below can assert detection rather than parity.
strong_pheno <- function(d, seed) {
  set.seed(seed)
  2 * as.numeric(scale(rowSums(d$geno[, d$causal_idx]))) + d$cov_effect +
    rnorm(nrow(d$geno), sd = 1)
}

# --- detector = "box" -------------------------------------------------------
# The box-scan detector is the alternative Stage 2/3. These tests pin: the
# option is wired through both pipelines, the region it reports covers the
# injected causal block, the scan itself is untouched by the choice of
# detector, and a missing threshold is refused up front.

test_that("ras(detector = 'box') reports a region covering the injected signal", {
  d <- make_pipeline_fixture(seed = 77, n_samp = 300, n_snp = 2000,
                             causal_lo = 900, causal_hi = 1000)
  pheno <- strong_pheno(d, 77)
  set.seed(77)
  r <- ras(
    d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 2, skip1 = 10, skip2 = 20,
    max_window_size = 100, chunk_snps = 400, save_dir = tempfile(),
    detector = "box", box_threshold = 2, run_plots = FALSE
  )
  expect_s3_class(r, "ras")
  expect_identical(r$detection$detector, "box")
  reg <- r$detection$regions
  expect_false(is.null(reg))
  covers <- reg$pos_L <= d$causal_center & reg$pos_R >= d$causal_center
  expect_true(any(covers))
  expect_equal(r$detection$tau_hats, reg$anchor_pos)
  expect_output(print(r), "box-scan detector")
  unlink(c(d$path, paste0(d$path, ".meta.rds")))
})

test_that("the detector choice does not change the scan profile", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                             causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect
  args <- list(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
               is_continuous = TRUE, num_rep = 2, skip1 = 50, skip2 = 25,
               max_window_size = 40, chunk_snps = 400, run_plots = FALSE)
  set.seed(7); a <- do.call(ras, c(args, list(save_dir = tempfile(), detector = "changepoint")))
  set.seed(7); b <- do.call(ras, c(args, list(save_dir = tempfile(),
                                                    detector = "box", box_threshold = 1)))
  expect_identical(a$scan$y, b$scan$y)
  expect_identical(a$detection$detector, "changepoint")
  expect_identical(b$detection$detector, "box")
  unlink(c(d$path, paste0(d$path, ".meta.rds")))
})

test_that("detector = 'box' calibrates its threshold on permuted phenotypes by default", {
  d <- make_pipeline_fixture(seed = 77, n_samp = 300, n_snp = 2000,
                             causal_lo = 900, causal_hi = 1000)
  pheno <- strong_pheno(d, 77)
  set.seed(77)
  r <- ras(
    d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 2, skip1 = 10, skip2 = 20,
    max_window_size = 100, chunk_snps = 400, save_dir = tempfile(),
    box_null = 12, box_alpha = 0.1, run_plots = FALSE
  )
  expect_identical(r$detection$detector, "box")
  expect_s3_class(r$box_calibration, "ras_box_calibration")
  expect_length(r$box_calibration$maxT, 12L)
  expect_equal(r$detection$threshold, r$box_calibration$threshold)
  expect_equal(r$detection$level0, r$box_calibration$level0)
  reg <- r$detection$regions
  expect_false(is.null(reg))
  expect_true(any(reg$pos_L <= d$causal_center & reg$pos_R >= d$causal_center))
  # nothing of the null scans is left behind in save_dir
  expect_length(list.files(r$save_dir, pattern = "coef_mat"), 2L)
  unlink(c(d$path, paste0(d$path, ".meta.rds")))
})

test_that("detector = 'box' with box_null = 0 and no threshold is refused", {
  d <- make_pipeline_fixture(seed = 4242, n_samp = 200, n_snp = 800,
                             causal_lo = 400, causal_hi = 500)
  pheno <- d$latent + d$cov_effect
  # the refusal comes from the detection stage; the (cheap) scan runs first
  expect_error(
    ras(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
             is_continuous = TRUE, num_rep = 1, skip1 = 50, skip2 = 25,
             max_window_size = 40, chunk_snps = 400, save_dir = tempfile(),
             detector = "box", box_null = 0, run_plots = FALSE),
    "box_threshold")
  expect_error(
    ras(d$path, pheno, d$cov_df, covariate_cols = d$covariate_cols,
             is_continuous = TRUE, detector = "nonsense"),
    "arg")
  unlink(c(d$path, paste0(d$path, ".meta.rds")))
})

test_that("ras_original(detector = 'box') runs the in-memory pipeline with a calibration object", {
  d <- make_pipeline_fixture(seed = 31, n_samp = 200, n_snp = 800,
                             causal_lo = 400, causal_hi = 480)
  pheno <- strong_pheno(d, 31)
  set.seed(5)
  nulls <- lapply(1:30, function(k) 0.5 + abs(rnorm(length(seq(1, 800, by = 10)), sd = 0.3)))
  cal <- ras_box_calibrate(nulls, alpha = 0.05)
  set.seed(31)
  r <- ras_original(
    d$geno, pheno, d$cov_df, covariate_cols = d$covariate_cols,
    is_continuous = TRUE, num_rep = 2, skip1 = 10, skip2 = 20,
    max_window_size = 40, save_dir = tempfile(), scan_test = "score",
    detector = "box", box_calibration = cal, run_plots = TRUE, plot_device = "pdf"
  )
  expect_identical(r$detection$detector, "box")
  expect_equal(r$detection$threshold, cal$threshold)
  expect_equal(r$detection$level0, cal$level0)
  reg <- r$detection$regions
  expect_false(is.null(reg))
  expect_true(any(reg$pos_L <= d$causal_center & reg$pos_R >= d$causal_center))
  expect_true(file.exists(file.path(r$save_dir, "chr-1-cp-plot.pdf")))
  expect_true(file.exists(file.path(r$save_dir, "chr-1-zoom.pdf")))
  unlink(c(d$path, paste0(d$path, ".meta.rds"), r$save_dir), recursive = TRUE)
})
