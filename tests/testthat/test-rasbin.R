make_rasbin_matrix <- function(n = 37, N = 53, na_frac = 0.05) {
  set.seed(7)
  m <- matrix(sample(0:2, n * N, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
              nrow = n, ncol = N)
  storage.mode(m) <- "double"
  m <- m + matrix(rnorm(n * N, sd = 0.01), n, N)  # non-integer dosages too
  na_idx <- sample(seq_len(n * N), size = floor(n * N * na_frac))
  m[na_idx] <- NA_real_
  rownames(m) <- paste0("S", seq_len(n))
  colnames(m) <- paste0("rs", seq_len(N))
  m
}

test_that("geno_to_rasbin + rasbin_header round-trips dimensions", {
  m <- make_rasbin_matrix()
  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(m, path)
  on.exit(unlink(c(path, paste0(path, ".meta.rds"))))

  hdr <- rasbin_header(path)
  expect_equal(unname(hdr["n_samples"]), nrow(m))
  expect_equal(unname(hdr["n_snps"]), ncol(m))
  expect_equal(unname(hdr["dtype"]), 0)
})

test_that("rasbin_read_chunk reproduces the full matrix exactly, including NA", {
  m <- make_rasbin_matrix()
  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(m, path)
  on.exit(unlink(c(path, paste0(path, ".meta.rds"))))

  full <- rasbin_read_chunk(path, 1, ncol(m))
  dimnames(full) <- NULL
  m_unnamed <- m
  dimnames(m_unnamed) <- NULL
  expect_identical(full, m_unnamed)
})

test_that("rasbin_read_chunk reproduces arbitrary column ranges, including boundaries", {
  m <- make_rasbin_matrix(n = 41, N = 97)
  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(m, path, chunk_write_size = 17)  # deliberately not a divisor of N
  on.exit(unlink(c(path, paste0(path, ".meta.rds"))))

  ranges <- list(c(1, 1), c(1, 5), c(50, 60), c(97, 97), c(90, 97), c(1, 97))
  for (r in ranges) {
    got <- rasbin_read_chunk(path, r[1], r[2])
    want <- m[, r[1]:r[2], drop = FALSE]
    dimnames(want) <- NULL
    expect_identical(got, want, info = paste("range", r[1], r[2]))
  }
})

test_that("rasbin_read_chunk rejects out-of-range column requests", {
  m <- make_rasbin_matrix(n = 10, N = 20)
  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(m, path)
  on.exit(unlink(c(path, paste0(path, ".meta.rds"))))

  expect_error(rasbin_read_chunk(path, 0, 5))
  expect_error(rasbin_read_chunk(path, 15, 25))
  expect_error(rasbin_read_chunk(path, 10, 5))
})

test_that("geno_to_rasbin writes a usable metadata sidecar", {
  m <- make_rasbin_matrix(n = 5, N = 8)
  path <- tempfile(fileext = ".rasbin")
  geno_to_rasbin(m, path)
  on.exit(unlink(c(path, paste0(path, ".meta.rds"))))

  meta <- readRDS(paste0(path, ".meta.rds"))
  expect_equal(meta$sample_ids, rownames(m))
  expect_equal(meta$snp_ids, colnames(m))
})
