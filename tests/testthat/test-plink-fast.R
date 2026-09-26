# Builds a synthetic PLINK 1 .bed/.bim/.fam fileset by hand (no PLINK binary
# needed) so bed_to_rasbin()'s byte-level unpacking can be validated directly
# against a known dosage matrix.
encode_bed_body <- function(m) {
  n <- nrow(m); N <- ncol(m)
  bytes_per_snp <- ceiling(n / 4)
  dosage_to_code <- c(`2` = 0L, `1` = 2L, `0` = 3L)  # NA -> 1L (missing)
  body <- raw(bytes_per_snp * N)
  idx <- 1L
  for (j in seq_len(N)) {
    col <- m[, j]
    codes <- vapply(col, function(v) if (is.na(v)) 1L else dosage_to_code[[as.character(v)]], integer(1))
    for (b in seq_len(bytes_per_snp)) {
      base <- (b - 1L) * 4L
      byte_val <- 0L
      for (k in 0:3) {
        i <- base + k + 1L
        code <- if (i <= n) codes[i] else 0L
        byte_val <- bitwOr(byte_val, bitwShiftL(code, 2L * k))
      }
      body[idx] <- as.raw(byte_val)
      idx <- idx + 1L
    }
  }
  body
}

make_test_plink <- function(dir, n = 41, N = 97, na_frac = 0.05, stem = "toy") {
  set.seed(11)
  m <- matrix(sample(0:2, n * N, replace = TRUE, prob = c(0.6, 0.3, 0.1)), nrow = n, ncol = N)
  storage.mode(m) <- "double"
  na_idx <- sample(seq_len(n * N), size = floor(n * N * na_frac))
  m[na_idx] <- NA_real_

  bed_path <- file.path(dir, paste0(stem, ".bed"))
  bim_path <- file.path(dir, paste0(stem, ".bim"))
  fam_path <- file.path(dir, paste0(stem, ".fam"))

  con <- file(bed_path, open = "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  writeBin(encode_bed_body(m), con)
  close(con)

  fam <- data.frame(FID = paste0("FAM", seq_len(n)), IID = paste0("IND", seq_len(n)),
                     PID = 0, MID = 0, SEX = 1, PHENO = -9)
  utils::write.table(fam, fam_path, quote = FALSE, row.names = FALSE, col.names = FALSE)

  bim <- data.frame(CHR = 1, SNP = paste0("rs", seq_len(N)), CM = 0,
                     POS = seq_len(N) * 1000, A1 = "A", A2 = "G")
  utils::write.table(bim, bim_path, quote = FALSE, row.names = FALSE, col.names = FALSE)

  list(bed = bed_path, bim = bim_path, fam = fam_path, m = m, fam_df = fam, bim_df = bim)
}

test_that("bed_to_rasbin round-trips dimensions and dtype", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir)
  out <- file.path(dir, "toy.rasbin")

  bed_to_rasbin(fx$bed, out)
  hdr <- rasbin_header(out)
  expect_equal(unname(hdr["n_samples"]), nrow(fx$m))
  expect_equal(unname(hdr["n_snps"]), ncol(fx$m))
  expect_equal(unname(hdr["dtype"]), 1)
})

test_that("bed_to_rasbin reproduces the full dosage matrix exactly, including NA", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir)
  out <- file.path(dir, "toy.rasbin")

  bed_to_rasbin(fx$bed, out)
  full <- rasbin_read_chunk(out, 1, ncol(fx$m))
  dimnames(full) <- NULL
  m_unnamed <- fx$m; dimnames(m_unnamed) <- NULL
  expect_identical(full, m_unnamed)
})

test_that("bed_to_rasbin reproduces arbitrary column ranges via rasbin_read_chunk, including boundaries", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 37, N = 53)
  out <- file.path(dir, "toy.rasbin")

  bed_to_rasbin(fx$bed, out, chunk_snps = 17)  # deliberately not a divisor of N
  ranges <- list(c(1, 1), c(1, 5), c(20, 30), c(53, 53), c(45, 53), c(1, 53))
  for (r in ranges) {
    got <- rasbin_read_chunk(out, r[1], r[2])
    want <- fx$m[, r[1]:r[2], drop = FALSE]
    dimnames(want) <- NULL
    expect_identical(got, want, info = paste("range", r[1], r[2]))
  }
})

test_that("bed_to_rasbin is insensitive to chunk_snps (small vs large vs single-chunk)", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 29, N = 40)

  outs <- lapply(c(3, 11, 10000), function(cs) {
    out <- file.path(dir, paste0("toy_", cs, ".rasbin"))
    bed_to_rasbin(fx$bed, out, chunk_snps = cs)
    full <- rasbin_read_chunk(out, 1, ncol(fx$m))
    dimnames(full) <- NULL
    full
  })
  expect_identical(outs[[1]], outs[[2]])
  expect_identical(outs[[1]], outs[[3]])
})

test_that("bed_to_rasbin writes a usable metadata sidecar from .fam/.bim", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 5, N = 8)
  out <- file.path(dir, "toy.rasbin")

  bed_to_rasbin(fx$bed, out)
  meta <- readRDS(paste0(out, ".meta.rds"))
  expect_equal(meta$sample_ids, paste(fx$fam_df$FID, fx$fam_df$IID, sep = "_"))
  expect_equal(meta$snp_ids, fx$bim_df$SNP)
})

test_that("bed_to_rasbin rejects a .bed whose size doesn't match .fam/.bim counts", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 10, N = 12)
  out <- file.path(dir, "toy.rasbin")

  # Truncate the .bed by one byte.
  raw_bytes <- readBin(fx$bed, "raw", n = file.size(fx$bed))
  writeBin(raw_bytes[-length(raw_bytes)], fx$bed)

  expect_error(bed_to_rasbin(fx$bed, out), "does not match")
})

test_that("bed_to_rasbin rejects a file with bad .bed magic bytes", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 10, N = 12)
  out <- file.path(dir, "toy.rasbin")

  con <- file(fx$bed, "r+b")
  writeBin(as.raw(0x00), con)
  close(con)

  expect_error(bed_to_rasbin(fx$bed, out))
})

test_that("bed_to_rasbin errors clearly when a companion .bim/.fam is missing", {
  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  fx <- make_test_plink(dir, n = 6, N = 6)
  unlink(fx$bim)

  expect_error(bed_to_rasbin(fx$bed, file.path(dir, "toy.rasbin")), "bim")
})
