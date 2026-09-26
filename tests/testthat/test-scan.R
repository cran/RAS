make_fast_scan_fixture <- function(seed, n_samp = 90, n_snp = 150, na_frac = 0.03) {
  set.seed(seed)
  geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
                 nrow = n_samp, ncol = n_snp)
  storage.mode(geno) <- "double"
  na_idx <- sample(seq_len(n_samp * n_snp), size = floor(n_samp * n_snp * na_frac))
  geno[na_idx] <- NA_real_

  this.leftout <- sort(sample(seq_len(n_samp), n_samp %/% 2))
  n_holdout <- length(this.leftout)
  weights <- rnorm(n_snp)

  df <- data.frame(
    sex = sample(c("Male", "Female"), n_holdout, replace = TRUE),
    age = sample(20:60, n_holdout, replace = TRUE)
  )
  df$age_squared <- df$age^2
  df$age_sex     <- df$age * as.numeric(df$sex == "Male")
  for (i in 1:6) df[[paste0("pc", i)]] <- rnorm(n_holdout)

  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(geno, path, chunk_write_size = 37)  # awkward chunk size on write, too

  list(geno = geno, this.leftout = this.leftout, weights = weights, df = df,
       path = path, n_snp = n_snp)
}

compare_fast_to_reference <- function(d, is_continuous, chunk_snps,
                                       covariate_formula = "sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6",
                                       skip1 = 7, skip2 = 6,
                                       min_window_size = 3, max_window_size = 18) {
  pgs.mat <- compute_pgs_matrix(d$geno, d$this.leftout, d$weights)

  ref <- screen_forward_max_region_original(
    geno = d$geno, pgs.mat = pgs.mat, this.df = d$df, num_signals = 0,
    isPlot = FALSE, skip1 = skip1, skip2 = skip2,
    min_window_size = min_window_size, max_window_size = max_window_size,
    is_continuous = is_continuous, covariate_formula = covariate_formula,
    scan_test = if (is_continuous) "glm" else "score"
  )

  fast <- screen_forward_max_region(
    rasbin_path = d$path, weights = d$weights, this.leftout = d$this.leftout,
    this.df = d$df, is_continuous = is_continuous,
    covariate_formula = covariate_formula,
    skip1 = skip1, skip2 = skip2,
    min_window_size = min_window_size, max_window_size = max_window_size,
    chunk_snps = chunk_snps
  )

  list(ref = ref, fast = fast)
}

test_that("fused C scan matches R reference for continuous trait across seeds", {
  for (seed in 1:5) {
    d <- make_fast_scan_fixture(seed)
    d$df$phenotype2 <- rnorm(length(d$this.leftout))
    res <- compare_fast_to_reference(d, is_continuous = TRUE, chunk_snps = 11)  # small, forces reloads
    unlink(c(d$path, paste0(d$path, ".meta.rds")))

    expect_equal(length(res$fast), length(res$ref))
    ok <- is.finite(res$ref) & is.finite(res$fast)
    expect_gt(mean(ok), 0.8)
    expect_equal(res$fast[ok], res$ref[ok], tolerance = 1e-6,
                 info = paste("seed", seed))
  }
})

test_that("fused C scan matches R reference (score test) for binary trait across seeds", {
  for (seed in 1:5) {
    d <- make_fast_scan_fixture(seed + 100)
    d$df$phenotype2 <- rbinom(length(d$this.leftout), 1, 0.5)
    res <- compare_fast_to_reference(d, is_continuous = FALSE, chunk_snps = 9)
    unlink(c(d$path, paste0(d$path, ".meta.rds")))

    expect_equal(length(res$fast), length(res$ref))
    ok <- is.finite(res$ref) & is.finite(res$fast)
    expect_gt(mean(ok), 0.8)
    expect_equal(res$fast[ok], res$ref[ok], tolerance = 1e-6,
                 info = paste("seed", seed))
  }
})

test_that("fused C scan is insensitive to chunk_snps (chunking is an implementation detail)", {
  d <- make_fast_scan_fixture(999)
  d$df$phenotype2 <- rnorm(length(d$this.leftout))

  common <- list(rasbin_path = d$path, weights = d$weights,
                 this.leftout = d$this.leftout, this.df = d$df,
                 is_continuous = TRUE,
                 covariate_formula = "sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6",
                 skip1 = 7, skip2 = 6, min_window_size = 3, max_window_size = 18)

  small_chunk <- do.call(screen_forward_max_region, c(common, chunk_snps = 5))
  big_chunk   <- do.call(screen_forward_max_region, c(common, chunk_snps = 10000))
  unlink(c(d$path, paste0(d$path, ".meta.rds")))

  expect_identical(small_chunk, big_chunk)
})

test_that("mode must be 'continuous' or a supported binary path", {
  d <- make_fast_scan_fixture(42, n_samp = 30, n_snp = 20)
  d$df$phenotype2 <- rnorm(length(d$this.leftout))
  on.exit(unlink(c(d$path, paste0(d$path, ".meta.rds"))))

  expect_no_error(screen_forward_max_region(
    rasbin_path = d$path, weights = d$weights, this.leftout = d$this.leftout,
    this.df = d$df, is_continuous = TRUE,
    covariate_formula = "sex + age",
    skip1 = 5, skip2 = 5, min_window_size = 2, max_window_size = 8
  ))
})
