#' Forward Scan of the RAS Profile
#'
#' The Stage-1 forward scan on the hold-out half of the sample. At every grid
#' position (every \code{skip1} SNPs) windows of \code{min_window_size},
#' \code{min_window_size + skip2}, ... up to \code{max_window_size} SNPs are
#' formed, each window's weighted dosage score is tested against the
#' phenotype, and the smallest p-value is recorded as the profile value. SNP
#' columns are streamed from the \code{.rasbin} file and the weighted scores
#' are accumulated on the fly, so neither the genotype matrix nor a score
#' matrix is held in memory. The pure-R, in-memory implementation is
#' \code{\link{screen_forward_max_region_original}} together with
#' \code{\link{compute_pgs_matrix}}.
#'
#' @param rasbin_path Character. Path to the \code{.rasbin} genotype file.
#' @param weights Numeric vector of length \eqn{N} (one weight per variant),
#'   typically the \code{Estimate} column of
#'   \code{\link{compute_gwas_weights}} or aligned external weights from
#'   \code{\link{ras_harmonize_sumstats}}.
#' @param this.leftout Integer vector. 1-based row indices of the samples to
#'   scan (the hold-out half, or all samples for external weights).
#' @param this.df Data frame with one row per element of \code{this.leftout},
#'   holding the phenotype in a column named \code{phenotype2} and the
#'   covariates named in \code{covariate_formula}.
#' @param is_continuous Logical. \code{TRUE}: exact Frisch-Waugh-Lovell
#'   regression of the window score on the phenotype after the covariates.
#'   \code{FALSE}: Rao score test of the window score in a logistic model
#'   with the covariates.
#' @param covariate_formula Character. Right-hand side listing the covariates,
#'   for example \code{"age + sex"}.
#' @param skip1 Integer. Grid stride in SNPs. Default \code{10}.
#' @param skip2 Integer. Window growth step in SNPs. Default \code{20}.
#' @param min_window_size,max_window_size Integer. Smallest and largest
#'   window in SNPs. Default \code{5} and \code{100}.
#' @param chunk_snps Integer. SNP columns read per disk chunk. Default
#'   \code{5000}.
#'
#' @return Numeric vector with one \eqn{-\log_{10}(p)} value per grid
#'   position, \code{ceiling(N / skip1)} values in total.
#'
#' @seealso \code{\link{ras_scan}}, which averages this scan over
#'   repetitions; \code{\link{ras_scan_external}};
#'   \code{\link{screen_forward_max_region_original}}.
#' @export
screen_forward_max_region <- function(rasbin_path, weights, this.leftout,
                                            this.df, is_continuous,
                                            covariate_formula,
                                            skip1 = 10, skip2 = 20,
                                            min_window_size = 5,
                                            max_window_size = 100,
                                            chunk_snps = 5000) {
  cov_form <- stats::as.formula(paste("~", covariate_formula))
  n_holdout <- length(this.leftout)

  # Completeness must be resolved on the data frame BEFORE model.matrix(), same
  # rationale as screen_forward_max_region_original(): model.matrix()'s na.action would
  # otherwise silently shorten the design, desynchronising it from the
  # full-length accumulator built from this.leftout.
  keep <- stats::complete.cases(
    this.df[, unique(c(all.vars(cov_form), "phenotype2")), drop = FALSE])
  Zmat <- stats::model.matrix(cov_form, data = this.df[keep, , drop = FALSE])
  yvec <- as.numeric(this.df$phenotype2[keep])

  if (is_continuous) {
    qrZ  <- qr(Zmat)
    resid_vec <- qr.resid(qrZ, yvec)
    SyyC <- sum(resid_vec * resid_vec)
    dfC  <- length(yvec) - (qrZ$rank + 1L)
    Mmat <- solve(crossprod(Zmat), t(Zmat))
    w0   <- rep(0, length(yvec))     # unused in continuous mode, passed for a uniform C signature
    mode <- "continuous"
  } else {
    fit0 <- stats::glm.fit(Zmat, yvec, family = stats::binomial())
    p0   <- fit0$fitted.values
    w0   <- p0 * (1 - p0)
    resid_vec <- yvec - p0
    Mmat <- solve(crossprod(Zmat, Zmat * w0), t(Zmat * w0))
    SyyC <- 0; dfC <- 0                # unused in score mode
    mode <- "score"
  }

  .Call("RAS_fused_pgs_scan", path.expand(rasbin_path),
        as.double(weights), as.integer(this.leftout), mode,
        as.integer(keep), Zmat, Mmat, as.double(resid_vec), as.double(w0),
        as.double(SyyC), as.double(dfC),
        as.integer(skip1), as.integer(skip2),
        as.integer(min_window_size), as.integer(max_window_size),
        as.integer(chunk_snps),
        PACKAGE = "RAS")
}
