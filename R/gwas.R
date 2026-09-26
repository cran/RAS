#' Per-SNP GWAS Effect Size Weights
#'
#' Fits the per-SNP regressions of RAS Stage 1 on the training half of the
#' sample: for each variant, \code{phenotype1 ~ dosage} for a continuous trait
#' or \code{phenotype1 ~ dosage + covariates} for a binary trait. SNP columns
#' are streamed from the \code{.rasbin} file in chunks and only the training
#' rows are kept, so the full genotype matrix is never held in memory. The
#' pure-R, in-memory implementation is
#' \code{\link{compute_gwas_weights_original}}.
#'
#' @param rasbin_path Character. Path to the \code{.rasbin} genotype file.
#' @param phenotype1 Numeric vector, one value per training sample: the
#'   covariate-adjusted residuals for a continuous trait, the raw 0/1 outcome
#'   for a binary trait.
#' @param this.sample Integer vector. 1-based row indices of the training
#'   samples in the genotype file.
#' @param this.df Data frame of covariates for all samples (indexed by
#'   \code{this.sample}); used for binary traits only.
#' @param is_continuous Logical. \code{TRUE} for a quantitative trait.
#' @param covariate_formula Character. Right-hand side of the binary-trait
#'   model, including the SNP term \code{this.x}, for example
#'   \code{"this.x + age + sex"}. Ignored for continuous traits.
#' @param chunk_snps Integer. SNP columns read per disk chunk. Default
#'   \code{5000}.
#'
#' @details
#' Binary traits are fitted by the Frisch-Waugh-Lovell identity: the covariate
#' projection is computed once and each SNP costs one residual regression.
#' SNPs with a missing dosage in the training rows are refitted on their
#' complete rows by downdating the shared covariate cross-products in blocks,
#' so post-QC sequencing data with sparse missingness stays fast.
#'
#' @return Numeric matrix with \eqn{N} rows (variants) and columns
#'   \code{Estimate}, \code{Std. Error}, \code{t value} and \code{Pr(>|t|)}.
#'   The \code{Estimate} column is the weight vector used by
#'   \code{\link{screen_forward_max_region}}.
#'
#' @seealso \code{\link{ras_scan}}, which calls this function once per
#'   repetition; \code{\link{compute_gwas_weights_original}}.
#' @export
compute_gwas_weights <- function(rasbin_path, phenotype1, this.sample,
                                       this.df, is_continuous,
                                       covariate_formula = NULL,
                                       chunk_snps = 5000) {
  if (is.null(covariate_formula)) {
    covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  }
  message("Starting GWAS ...")

  if (is_continuous) {
    coef.mat <- .Call("RAS_chunked_gwas_continuous", path.expand(rasbin_path),
                       as.integer(this.sample), as.double(phenotype1),
                       as.integer(chunk_snps), PACKAGE = "RAS")
  } else {
    cov_rhs  <- .drop_term(covariate_formula, "this.x")
    cov_form <- stats::as.formula(paste("~", cov_rhs))

    base_df <- this.df[this.sample, , drop = FALSE]
    base_df$phenotype1 <- as.numeric(phenotype1)
    base_ok <- stats::complete.cases(base_df)
    base_df <- base_df[base_ok, , drop = FALSE]

    Wb   <- stats::model.matrix(cov_form, data = base_df)
    qrWb <- qr(Wb)
    yb   <- base_df$phenotype1
    ryb  <- qr.resid(qrWb, yb)
    pW   <- qrWb$rank
    nB   <- length(yb)
    dfree <- nB - (pW + 1L)
    Mmat  <- solve(crossprod(Wb), t(Wb))
    Syy   <- sum(ryb * ryb)

    train_ok_idx <- this.sample[base_ok]

    res <- .Call("RAS_chunked_gwas_binary_clean", path.expand(rasbin_path),
                 as.integer(train_ok_idx), Mmat, Wb, as.double(ryb),
                 as.double(Syy), as.double(dfree), as.integer(chunk_snps),
                 PACKAGE = "RAS")
    coef.mat <- res$coef_mat

    # SNP columns carrying an NA within the training/base_ok rows: each one is
    # fitted on its own complete rows. This is NOT a rare fallback on real data
    # -- with post-QC dosages nearly every column has at least one missing call
    # (All of Us chromosome 16: 99.9% of columns, and this path was ~95% of a binary
    # replicate), so it is batched, not looped one column at a time:
    #   * columns are re-read in blocks instead of one seek per column;
    #   * the per-column fit DOWNDATES the shared covariate cross-products --
    #     for a column missing the row set S, (W'W)_S = W'W - W_S'W_S and
    #     likewise for W'y, W'x, x'x, x'y, y'y -- so a column costs one
    #     (pW x pW) solve instead of a model.matrix() + qr() over all rows.
    # Same linear algebra as the per-column refit, so the coefficients match it
    # exactly; a column whose reduced design loses rank (a dropped factor level)
    # makes solve() fail and falls back to the original per-column path.
    bad <- res$na_cols
    if (length(bad) > 0L) {
      A  <- crossprod(Wb)
      bY <- as.numeric(crossprod(Wb, yb))
      yy <- sum(yb * yb)

      fit_one_refit <- function(jcol, xi) {          # original path, used on rank loss
        keep <- !is.na(xi)
        sub  <- base_df[keep, , drop = FALSE]
        ni   <- nrow(sub)
        if (ni < 4L) return(NULL)
        Wi   <- stats::model.matrix(cov_form, data = sub)
        qrWi <- qr(Wi)
        dfi  <- ni - (qrWi$rank + 1L)
        if (dfi <= 0L) return(NULL)
        ryi  <- qr.resid(qrWi, sub$phenotype1)
        rxi  <- qr.resid(qrWi, xi[keep])
        srxx <- sum(rxi * rxi)
        if (!is.finite(srxx) || srxx <= 0) return(NULL)
        srxy <- sum(rxi * ryi)
        b    <- srxy / srxx
        rss  <- max(sum(ryi * ryi) - b * srxy, 0)
        se   <- sqrt((rss / dfi) / srxx)
        c(b, se, b / se, 2 * stats::pt(-abs(b / se), df = dfi))
      }

      # Block size is capped by BYTES, not by chunk_snps: rasbin_read_chunk()
      # returns all n rows of the block as doubles, and with cores > 1 every
      # worker holds one such block at the same time. At full-cohort scale
      # (453,698 samples) chunk_snps = 2,000 would be 7.3 GB per worker, which
      # is what killed a worker before this cap existed. ~1 GB per block keeps
      # five workers inside a few GB, and the block only feeds one matmul per
      # column, so smaller blocks cost nothing but a few extra reads.
      n_total  <- as.integer(rasbin_header(rasbin_path)[["n_samples"]])
      per_blk  <- max(1L, min(as.integer(chunk_snps), as.integer(1.25e8 %/% max(1L, n_total))))
      blocks   <- split(bad, ceiling(seq_along(bad) / per_blk))
      for (blk in blocks) {
        Xr <- rasbin_read_chunk(rasbin_path, min(blk), max(blk))
        for (j in seq_along(blk)) {
          xj  <- Xr[train_ok_idx, blk[j] - min(blk) + 1L]   # one column, base rows
          S   <- which(is.na(xj))
          ni  <- nB - length(S)
          dfi <- ni - (pW + 1L)
          if (ni < 4L || dfi <= 0L) next
          xj0 <- xj; xj0[S] <- 0          # zeroed rows drop out of every sum below
          bXs <- as.numeric(crossprod(Wb, xj0))
          Ws  <- Wb[S, , drop = FALSE]
          ys  <- yb[S]
          bYs <- bY - as.numeric(crossprod(Ws, ys))
          sol <- try(solve(A - crossprod(Ws), cbind(bXs, bYs)), silent = TRUE)
          if (inherits(sol, "try-error")) {
            v <- fit_one_refit(blk[j], xj)
            if (!is.null(v)) coef.mat[blk[j], ] <- v
            next
          }
          srxx <- sum(xj0 * xj0) - sum(bXs * sol[, 1])
          if (!is.finite(srxx) || srxx <= 0) next
          srxy <- sum(xj0 * yb) - sum(bYs * sol[, 1])
          sryy <- (yy - sum(ys * ys)) - sum(bYs * sol[, 2])
          b    <- srxy / srxx
          rss  <- max(sryy - b * srxy, 0)
          se   <- sqrt((rss / dfi) / srxx)
          tt   <- b / se
          coef.mat[blk[j], ] <- c(b, se, tt, 2 * stats::pt(-abs(tt), df = dfi))
        }
        rm(Xr)
      }
    }
  }

  colnames(coef.mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coef.mat
}
