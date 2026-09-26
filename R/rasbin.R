#' Write a Genotype Matrix to the Chunked .rasbin Format
#'
#' Converts an in-memory genotype dosage matrix (or a matrix stored in an
#' \code{.rds}/\code{.csv}/\code{.tsv} file) into \code{.rasbin}: a flat,
#' column-major binary file that supports reading an arbitrary SNP-column
#' range with a single seek, without loading the whole matrix into memory.
#'
#' @param geno Either a numeric matrix (\eqn{n} samples \eqn{\times} \eqn{N}
#'   variants), or a character path to a \code{.rds}, \code{.csv}, or
#'   \code{.tsv}/\code{.txt} file containing one (first column treated as
#'   sample ID if present as rownames-style data). The source is read into
#'   memory exactly once, regardless of downstream chunked access.
#' @param out_path Character. Destination path for the \code{.rasbin} file.
#'   A companion \code{<out_path>.meta.rds} is written alongside it holding
#'   sample and SNP identifiers.
#' @param chunk_write_size Integer. Number of SNP columns written per
#'   \code{writeBin} call. Only affects write throughput, not the resulting
#'   file. Default \code{2000}.
#'
#' @details
#' File layout (see \code{src/rasbin.h} for the authoritative spec):
#' \itemize{
#'   \item Header (32 bytes): 8-byte magic \code{"RASBIN01"}, 8-byte
#'     \code{n_samples} (int64), 8-byte \code{n_snps} (int64), 4-byte
#'     \code{dtype} (int32; \code{0} = double), 4 reserved bytes.
#'   \item Body: \code{n_samples * n_snps} doubles, column-major, so column
#'     \code{j} (0-based) occupies bytes
#'     \code{[32 + j*n_samples*8, 32 + (j+1)*n_samples*8)}.
#' }
#' \code{NA} dosages round-trip exactly: they are ordinary IEEE-754 doubles
#' at the bit level, so a raw byte copy preserves R's \code{NA_real_}
#' sentinel without any special-casing.
#'
#' @return Invisibly, the \code{out_path} that was written.
#' @export
geno_to_rasbin <- function(geno, out_path, chunk_write_size = 2000) {
  if (is.character(geno)) {
    ext <- tolower(tools::file_ext(geno))
    raw <- switch(ext,
      rds = readRDS(geno),
      csv = utils::read.csv(geno, stringsAsFactors = FALSE, check.names = FALSE),
      tsv = ,
      txt = utils::read.delim(geno, stringsAsFactors = FALSE, check.names = FALSE),
      stop("geno_to_rasbin: unsupported source extension '", ext, "'"))
    if (is.data.frame(raw)) {
      rownames(raw) <- raw[[1]]
      raw <- as.matrix(raw[, -1, drop = FALSE])
      mode(raw) <- "numeric"
    }
    geno <- raw
  }
  if (!is.matrix(geno)) stop("geno_to_rasbin: geno must resolve to a matrix")
  storage.mode(geno) <- "double"

  n <- nrow(geno)
  N <- ncol(geno)
  sample_ids <- rownames(geno)
  snp_ids    <- colnames(geno)

  con <- file(out_path, open = "wb")
  on.exit(close(con), add = TRUE)

  writeBin(charToRaw("RASBIN01"), con, useBytes = TRUE)
  writeBin(as.integer(n), con, size = 8L)   # int64 n_samples
  writeBin(as.integer(N), con, size = 8L)   # int64 n_snps
  writeBin(0L,            con, size = 4L)   # int32 dtype = 0 (double)
  writeBin(0L,            con, size = 4L)   # reserved

  starts <- seq(1L, N, by = chunk_write_size)
  for (s in starts) {
    e <- min(s + chunk_write_size - 1L, N)
    writeBin(as.vector(geno[, s:e, drop = FALSE]), con, size = 8L)
  }

  saveRDS(list(sample_ids = sample_ids, snp_ids = snp_ids),
          paste0(out_path, ".meta.rds"))

  invisible(out_path)
}

#' Read the .rasbin File Header
#'
#' @param path Character. Path to a \code{.rasbin} file.
#' @return Named numeric vector with elements \code{n_samples}, \code{n_snps},
#'   \code{dtype}.
#' @export
rasbin_header <- function(path) {
  out <- .Call("RAS_rasbin_header", path.expand(path), PACKAGE = "RAS")
  stats::setNames(out, c("n_samples", "n_snps", "dtype"))
}

#' Read a Column Range from a .rasbin File
#'
#' @param path Character. Path to a \code{.rasbin} file.
#' @param col_start Integer. First column to read (1-based, inclusive).
#' @param col_end Integer. Last column to read (1-based, inclusive).
#' @return An \code{n_samples x (col_end - col_start + 1)} numeric matrix.
#' @export
rasbin_read_chunk <- function(path, col_start, col_end) {
  .Call("RAS_rasbin_read_chunk", path.expand(path),
        as.double(col_start), as.double(col_end), PACKAGE = "RAS")
}

# Resolve the `geno` argument of ras(), ras_scan() and ras_scan_external():
# a path to an existing .rasbin file is used as is; an in-memory genotype
# matrix (or data frame) is written to a temporary .rasbin file that the
# caller removes on exit through the returned cleanup function. This keeps
# the 1.0.x calling convention `ras(geno = <matrix>, ...)` working on the
# disk-backed engine.
.ras_resolve_geno <- function(geno, what = "geno") {
  if (is.character(geno)) {
    if (length(geno) != 1L || !file.exists(geno))
      stop(what, ": file not found: ", paste(geno, collapse = ", "), call. = FALSE)
    return(list(path = geno, cleanup = function() invisible(NULL)))
  }
  if (is.data.frame(geno)) geno <- as.matrix(geno)
  if (!is.matrix(geno))
    stop(what, " must be the path to a .rasbin file (see geno_to_rasbin() / ",
         "bed_to_rasbin()) or a numeric genotype matrix", call. = FALSE)
  tmp <- tempfile(fileext = ".rasbin")
  message(sprintf(paste0(
    "Converting the %s x %s in-memory genotype matrix to a temporary .rasbin file. ",
    "For repeated runs convert once with geno_to_rasbin() and pass the path."),
    format(nrow(geno), big.mark = ","), format(ncol(geno), big.mark = ",")))
  geno_to_rasbin(geno, tmp)
  list(path = tmp, cleanup = function() unlink(c(tmp, paste0(tmp, ".meta.rds"))))
}
