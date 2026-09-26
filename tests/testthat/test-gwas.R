make_gwas_fast_fixture <- function(seed, n_samp = 120, n_snp = 90, na_frac = 0.04) {
  set.seed(seed)
  geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
                 nrow = n_samp, ncol = n_snp)
  storage.mode(geno) <- "double"
  na_idx <- sample(seq_len(n_samp * n_snp), size = floor(n_samp * n_snp * na_frac))
  geno[na_idx] <- NA_real_

  this.sample <- sort(sample(seq_len(n_samp), n_samp %/% 2))

  df <- data.frame(
    sex = sample(c("Male", "Female"), n_samp, replace = TRUE),
    age = sample(20:60, n_samp, replace = TRUE)
  )
  df$age_squared <- df$age^2
  df$age_sex     <- df$age * as.numeric(df$sex == "Male")
  for (i in 1:5) df[[paste0("pc", i)]] <- rnorm(n_samp)

  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(geno, path, chunk_write_size = 13)

  list(geno = geno, this.sample = this.sample, df = df, path = path, n_snp = n_snp)
}

expect_coefmat_close <- function(fast, ref, tolerance = 1e-6) {
  expect_equal(dim(fast), dim(ref))
  na_ref  <- is.na(ref)
  na_fast <- is.na(fast)
  expect_identical(na_fast, na_ref, info = "NA pattern must match exactly")
  expect_equal(fast[!na_ref], ref[!na_ref], tolerance = tolerance)
}

test_that("compute_gwas_weights matches R reference for continuous trait, incl. NA fallback", {
  for (seed in 1:5) {
    d <- make_gwas_fast_fixture(seed)
    phenotype1 <- rnorm(length(d$this.sample))

    ref <- compute_gwas_weights_original(d$geno, phenotype1, d$this.sample, d$df, is_continuous = TRUE)
    fast <- compute_gwas_weights(d$path, phenotype1, d$this.sample, d$df,
                                       is_continuous = TRUE, chunk_snps = 7)
    unlink(c(d$path, paste0(d$path, ".meta.rds")))

    expect_coefmat_close(fast, ref, tolerance = 1e-8)
  }
})

test_that("compute_gwas_weights matches R reference for binary trait, incl. NA fallback", {
  for (seed in 1:5) {
    d <- make_gwas_fast_fixture(seed + 200)
    phenotype1 <- rbinom(length(d$this.sample), 1, 0.4)
    covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5"

    ref <- compute_gwas_weights_original(d$geno, phenotype1, d$this.sample, d$df,
                                 is_continuous = FALSE, covariate_formula = covariate_formula)
    fast <- compute_gwas_weights(d$path, phenotype1, d$this.sample, d$df,
                                       is_continuous = FALSE, covariate_formula = covariate_formula,
                                       chunk_snps = 6)
    unlink(c(d$path, paste0(d$path, ".meta.rds")))

    expect_coefmat_close(fast, ref, tolerance = 1e-6)
  }
})

test_that("compute_gwas_weights is insensitive to chunk_snps", {
  d <- make_gwas_fast_fixture(777)
  phenotype1 <- rbinom(length(d$this.sample), 1, 0.5)
  covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5"

  small <- compute_gwas_weights(d$path, phenotype1, d$this.sample, d$df,
                                      is_continuous = FALSE, covariate_formula = covariate_formula,
                                      chunk_snps = 3)
  big   <- compute_gwas_weights(d$path, phenotype1, d$this.sample, d$df,
                                      is_continuous = FALSE, covariate_formula = covariate_formula,
                                      chunk_snps = 10000)
  unlink(c(d$path, paste0(d$path, ".meta.rds")))

  expect_identical(small, big)
})

test_that("compute_gwas_weights handles a genotype-NA-heavy training subset (forces the fallback path)", {
  d <- make_gwas_fast_fixture(55, n_samp = 80, n_snp = 40, na_frac = 0.15)
  phenotype1 <- rbinom(length(d$this.sample), 1, 0.5)
  covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5"

  ref <- compute_gwas_weights_original(d$geno, phenotype1, d$this.sample, d$df,
                               is_continuous = FALSE, covariate_formula = covariate_formula)
  fast <- compute_gwas_weights(d$path, phenotype1, d$this.sample, d$df,
                                     is_continuous = FALSE, covariate_formula = covariate_formula,
                                     chunk_snps = 5)
  unlink(c(d$path, paste0(d$path, ".meta.rds")))

  expect_coefmat_close(fast, ref, tolerance = 1e-6)
})
