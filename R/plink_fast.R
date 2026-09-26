#' Convert a PLINK 1 .bed/.bim/.fam Fileset to .rasbin
#'
#' Converts a PLINK 1 binary fileset directly into the chunked \code{.rasbin}
#' format, without ever materializing a dense genotype matrix in R. Both the
#' read from \code{.bed} and the write to \code{.rasbin} are done in fixed-size
#' SNP chunks in C (see \code{src/plink_io.c}), so peak memory is
#' \eqn{O(\text{n\_samples} \times \text{chunk\_snps})} regardless of how many
#' samples or SNPs the fileset contains.
#'
#' @param bed_path Character. Path to the \code{.bed} file. The companion
#'   \code{.bim} and \code{.fam} files are expected alongside it with the same
#'   stem (standard PLINK convention, e.g. \code{--bfile} prefix).
#' @param out_path Character. Destination path for the \code{.rasbin} file. A
#'   companion \code{<out_path>.meta.rds} is written alongside it holding
#'   sample and SNP identifiers, parsed from \code{.fam}/\code{.bim}.
#' @param chunk_snps Integer. Number of SNPs read from \code{.bed} and written
#'   to \code{.rasbin} per chunk. Only affects peak memory and throughput, not
#'   the resulting file. Default \code{2000}, matching
#'   \code{\link{geno_to_rasbin}}'s \code{chunk_write_size} default.
#'
#' @details
#' Only SNP-major \code{.bed} files are supported (the only mode any current
#' PLINK version writes). Genotype dosage counts copies of the A1 allele (the
#' fourth column of \code{.bim}), matching PLINK's own \code{--recode A} /
#' plink2's \code{--export A} convention: \code{2} = homozygous A1, \code{1} =
#' heterozygous, \code{0} = homozygous A2, \code{NA} = missing.
#'
#' Unlike \code{\link{geno_to_rasbin}} (dtype 0, double, 8 bytes/genotype),
#' this writes dtype 1 (int8, 1 byte/genotype) -- 8x smaller on disk, though
#' still ~4x larger than \code{.bed}'s 2-bit packing (see \code{src/rasbin.h}
#' for why byte-per-genotype was chosen over bit-packing). Both dtypes are
#' read transparently by \code{\link{rasbin_read_chunk}} and every
#' \verb{_fast} pipeline function.
#'
#' @return Invisibly, the \code{out_path} that was written.
#' @export
bed_to_rasbin <- function(bed_path, out_path, chunk_snps = 2000) {
  stem <- sub("\\.bed$", "", bed_path)
  bim_path <- paste0(stem, ".bim")
  fam_path <- paste0(stem, ".fam")
  if (!file.exists(bed_path)) stop("bed_to_rasbin: '", bed_path, "' not found")
  if (!file.exists(bim_path)) stop("bed_to_rasbin: companion '", bim_path, "' not found")
  if (!file.exists(fam_path)) stop("bed_to_rasbin: companion '", fam_path, "' not found")

  fam <- utils::read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
  bim <- utils::read.table(bim_path, header = FALSE, stringsAsFactors = FALSE)
  n_samples <- nrow(fam)
  n_snps <- nrow(bim)
  if (n_samples < 1) stop("bed_to_rasbin: '", fam_path, "' has no samples")
  if (n_snps < 1) stop("bed_to_rasbin: '", bim_path, "' has no variants")

  bytes_per_snp <- ceiling(n_samples / 4)
  expected_size <- 3 + n_snps * bytes_per_snp
  actual_size <- file.size(bed_path)
  if (is.na(actual_size) || actual_size != expected_size) {
    stop("bed_to_rasbin: '", bed_path, "' size (", actual_size, " bytes) does not match ",
         "n_samples=", n_samples, " x n_snps=", n_snps, " (expected ", expected_size,
         " bytes) -- .fam/.bim counts don't match this .bed")
  }

  .Call("RAS_bed_to_rasbin", path.expand(bed_path), as.double(n_samples),
        as.double(n_snps), path.expand(out_path), as.double(chunk_snps),
        PACKAGE = "RAS")

  sample_ids <- paste(fam[[1]], fam[[2]], sep = "_")
  snp_ids <- bim[[2]]
  saveRDS(list(sample_ids = sample_ids, snp_ids = snp_ids),
          paste0(out_path, ".meta.rds"))

  invisible(out_path)
}
