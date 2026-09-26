# Internal helper: which columns of a matrix contain at least one NA.
.col_any_na <- function(M) {
  (colSums(is.na(M)) > 0L)
}

# Internal helper: drop a term (e.g. "this.x") from a formula RHS string,
# robust to spacing/order. Returns "1" if nothing remains.
.drop_term <- function(rhs, drop) {
  tl <- attr(stats::terms(stats::as.formula(paste("~", rhs))), "term.labels")
  tl <- setdiff(tl, drop)
  if (length(tl) == 0L) return("1")
  paste(tl, collapse = " + ")
}

#' Per-SNP GWAS Effect Size Weights (Pure-R, In-Memory)
#'
#' This is the pure-R, in-memory implementation that RAS 1.0.x shipped under
#' the plain name. Since 1.1.0 the plain name (\code{\link{compute_gwas_weights}}) is the
#' compiled, disk-backed implementation; this function is kept for reference
#' and gives the same results.
#' Runs a per-SNP regression on the training split to obtain effect-size
#' estimates used as polygenic score (PGS) weights in the forward scan.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants). Only the rows indexed by \code{this.sample} are used.
#' @param phenotype1 Numeric vector of length \eqn{|\texttt{this.sample}|}.
#'   Training-split phenotype values.  For continuous traits these should be
#'   covariate-residualised values; for binary traits, the raw 0/1 indicator.
#' @param this.sample Integer vector. Row indices of the training split in
#'   \code{geno}.
#' @param this.df Data frame of covariates with rows aligned to \code{geno}.
#'   Required for binary traits; used to subset and build the regression data
#'   frame.  Ignored for continuous traits.
#' @param is_continuous Logical. \code{TRUE} for quantitative traits,
#'   \code{FALSE} for binary (case/control) traits.
#' @param covariate_formula Character. The right-hand side of the regression
#'   formula including the SNP term \code{this.x}.  For example,
#'   \code{"this.x + age + sex + pc1"}.  If \code{NULL}, a default formula
#'   with sex, age, age-squared, age-sex interaction, and the first ten
#'   ancestry PCs is used.
#'
#' @details
#' For \strong{continuous traits} (\code{is_continuous = TRUE}) the function
#' fits \code{lm(phenotype1 ~ this.x)} for each SNP, excluding samples with
#' missing dosage via \code{na.omit}.  The covariates are assumed to have
#' already been residualised from \code{phenotype1} before this function is
#' called (as done in \code{\link{ras_scan_original}}).
#'
#' For \strong{binary traits} (\code{is_continuous = FALSE}) the function
#' builds a data frame combining the training-split rows of \code{this.df}
#' with the dosage column (\code{this.x}) and phenotype (\code{phenotype1}),
#' removes incomplete cases with \code{na.omit}, and fits
#' \code{lm(phenotype1 ~ covariate_formula)}.  Note that \code{lm} rather
#' than \code{glm(family = binomial)} is used at the GWAS weight stage
#' to obtain a linear probability approximation to the effect size; this is
#' consistent with standard GWAS practice for weight derivation.
#'
#' Both branches are vectorised: instead of fitting one \code{lm} per SNP, all
#' NA-free SNPs are solved in closed form with a single BLAS cross-product
#' (continuous) or a batched Frisch-Waugh-Lovell residualisation (binary).
#' SNP columns that contain missing dosages are handled individually so the
#' result matches \code{lm} + \code{na.omit} exactly.  A single
#' \dQuote{Starting GWAS ...} message is printed.
#'
#' @return
#' A numeric matrix with \eqn{N} rows and 4 columns:
#' \describe{
#'   \item{\code{Estimate}}{Per-SNP effect-size estimate (slope of the dosage
#'     term).  The first column is used as PGS weights by
#'     \code{\link{compute_pgs_matrix}}.}
#'   \item{\code{Std. Error}}{Standard error of the estimate.}
#'   \item{\code{t value}}{t-statistic.}
#'   \item{\code{Pr(>|t|)}}{Two-sided p-value.}
#' }
#'
#' @seealso
#' \code{\link{compute_pgs_matrix}} for the next step, which multiplies
#' these weights by the hold-out genotype matrix.
#' \code{\link{ras_scan_original}} for the recommended high-level entry point that
#' calls this function automatically.
#'
#' @examples
#' set.seed(1)
#' n_samp <- 60; n_snp <- 10
#' geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE,
#'                       prob = c(0.6, 0.3, 0.1)),
#'                nrow = n_samp, ncol = n_snp)
#' pheno_train <- rnorm(n_samp)
#'
#' coef_mat <- compute_gwas_weights_original(
#'   geno              = geno,
#'   phenotype1        = pheno_train,
#'   this.sample       = seq_len(n_samp),
#'   this.df           = data.frame(row.names = seq_len(n_samp)),
#'   is_continuous     = TRUE
#' )
#' dim(coef_mat)        # n_snp x 4
#' head(coef_mat, 3)
#' @export
compute_gwas_weights_original <- function(geno, phenotype1, this.sample, this.df,
                                 is_continuous,
                                 covariate_formula = NULL) {
  if (is.null(covariate_formula)) {
    covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  }

  N <- ncol(geno)
  coef.mat <- matrix(NA_real_, nrow = N, ncol = 4)
  message("Starting GWAS ...")

  # Training-split genotype submatrix (m x N), shared by both branches.
  Xg <- geno[this.sample, , drop = FALSE]
  m  <- length(this.sample)

  if (is_continuous) {
    # ── Vectorized simple regression  y ~ x_j  (with intercept), NA-dropped ──
    # phenotype1 is already covariate-residualised, so each SNP model is the
    # closed-form OLS of y on x_j. All NA-free SNPs are solved with one BLAS
    # matmul (crossprod), reproducing lm()'s Estimate/SE/t/p exactly (df = n-2).
    y      <- as.numeric(phenotype1)
    na_col <- .col_any_na(Xg)

    # (a) NA-free columns: a single BLAS matmul covers all of them at once.
    clean <- which(!na_col)
    if (length(clean) > 0L) {
      Xc  <- Xg[, clean, drop = FALSE]
      yc  <- y - mean(y)
      Xc  <- sweep(Xc, 2L, colMeans(Xc), "-")
      Sxx <- colSums(Xc * Xc)
      Sxy <- as.numeric(crossprod(Xc, yc))   # t(Xc) %*% yc  — the heavy step
      Syy <- sum(yc * yc)
      df  <- m - 2L
      beta <- Sxy / Sxx
      RSS  <- Syy - beta * Sxy
      RSS[RSS < 0] <- 0
      SE   <- sqrt((RSS / df) / Sxx)
      tval <- beta / SE
      pval <- 2 * stats::pt(-abs(tval), df = df)
      ok   <- is.finite(Sxx) & Sxx > 0 & df > 0L   # guard monomorphic columns
      block <- cbind(beta, SE, tval, pval)
      block[!ok, ] <- NA_real_
      coef.mat[clean, ] <- block
    }

    # (b) columns with NA: handled individually to match lm + na.omit exactly.
    for (i in which(na_col)) {
      xi   <- Xg[, i]
      keep <- !is.na(xi)
      ni   <- sum(keep)
      if (ni < 3L) next
      xk <- xi[keep]; yk <- y[keep]
      xk <- xk - mean(xk); yk <- yk - mean(yk)
      sxx <- sum(xk * xk)
      if (!is.finite(sxx) || sxx <= 0) next
      sxy <- sum(xk * yk); syy <- sum(yk * yk)
      dfi <- ni - 2L
      b   <- sxy / sxx
      rss <- max(syy - b * sxy, 0)
      se  <- sqrt((rss / dfi) / sxx)
      tt  <- b / se
      coef.mat[i, ] <- c(b, se, tt, 2 * stats::pt(-abs(tt), df = dfi))
    }

  } else {
    # ── Binary trait: lm(y ~ this.x + covariates) via Frisch–Waugh–Lovell ────
    # Residualise y and each SNP against the covariate design W, then do a
    # simple regression on the residuals. FWL guarantees the SNP coefficient,
    # its SE, t and p equal the `this.x` row of the full lm fit, provided the
    # full-model residual df = n - (rank(W) + 1) is used.
    cov_rhs  <- .drop_term(covariate_formula, "this.x")
    cov_form <- stats::as.formula(paste("~", cov_rhs))

    # Base completeness over the training split: drop rows with NA in the
    # phenotype or ANY column of this.df — matches the original
    # na.omit(this.df.sub) contract (this.df.sub = this.df cols + this.x + pheno).
    base_df <- this.df[this.sample, , drop = FALSE]
    base_df$phenotype1 <- as.numeric(phenotype1)
    base_ok <- stats::complete.cases(base_df)
    base_df <- base_df[base_ok, , drop = FALSE]

    Wb   <- stats::model.matrix(cov_form, data = base_df)  # incl. intercept
    qrWb <- qr(Wb)
    yb   <- base_df$phenotype1
    ryb  <- qr.resid(qrWb, yb)
    pW   <- qrWb$rank                  # rank of covariate design (incl intercept)
    nB   <- length(yb)

    Xb         <- Xg[base_ok, , drop = FALSE]
    na_in_base <- .col_any_na(Xb)

    # (a) SNP columns complete on the base rows: shared QR + batched residualize.
    clean <- which(!na_in_base)
    if (length(clean) > 0L) {
      rX   <- qr.resid(qrWb, Xb[, clean, drop = FALSE])
      Srxx <- colSums(rX * rX)
      Srxy <- as.numeric(crossprod(rX, ryb))
      Syy  <- sum(ryb * ryb)
      df   <- nB - (pW + 1L)           # n - rank(full design)
      beta <- Srxy / Srxx
      RSS  <- Syy - beta * Srxy
      RSS[RSS < 0] <- 0
      SE   <- sqrt((RSS / df) / Srxx)
      tval <- beta / SE
      pval <- 2 * stats::pt(-abs(tval), df = df)
      ok   <- is.finite(Srxx) & Srxx > 0 & df > 0L
      block <- cbind(beta, SE, tval, pval)
      block[!ok, ] <- NA_real_
      coef.mat[clean, ] <- block
    }

    # (b) SNP columns with NA dosage: same estimator as a per-column re-fit on
    #     that column's own rows, but obtained by DOWNDATING the shared
    #     covariate cross-products instead of rebuilding model.matrix() + qr()
    #     for every column. For a column missing the row set S,
    #     (W'W)_S = W'W - W_S'W_S and likewise for W'y, W'x, x'x, x'y, y'y, so
    #     each column costs one (pW x pW) solve instead of an n x pW QR. With
    #     post-QC genotypes almost every column has at least one missing dosage,
    #     so this path -- not the batched one above -- sets the cost of the
    #     binary GWAS stage: measured on an All of Us chromosome 16 run it was ~95% of a replicate.
    #     Exactness: identical linear algebra, so results match the per-column
    #     refit to floating-point error (see tests/testthat/test-gwas-binary.R).
    bad <- which(na_in_base)
    if (length(bad) > 0L) {
      A    <- crossprod(Wb)                          # W'W over all base rows
      bY   <- as.numeric(crossprod(Wb, yb))          # W'y
      yy   <- sum(yb * yb)
      X0   <- Xb[, bad, drop = FALSE]
      miss <- is.na(X0)
      X0[miss] <- 0                    # zeroed rows contribute nothing to any sum
      BX   <- crossprod(Wb, X0)        # W'x for every column with NA, one matmul
      xx   <- colSums(X0 * X0)
      xy   <- as.numeric(crossprod(X0, yb))
      for (j in seq_along(bad)) {
        S   <- which(miss[, j])
        ni  <- nB - length(S)
        dfi <- ni - (pW + 1L)
        if (ni < 4L || dfi <= 0L) next
        Ws  <- Wb[S, , drop = FALSE]
        ys  <- yb[S]
        As  <- A - crossprod(Ws)
        bYs <- bY - as.numeric(crossprod(Ws, ys))
        bXs <- BX[, j]
        sol <- try(solve(As, cbind(bXs, bYs)), silent = TRUE)
        if (inherits(sol, "try-error")) next        # rank dropped with those rows
        srxx <- xx[j] - sum(bXs * sol[, 1])
        if (!is.finite(srxx) || srxx <= 0) next
        srxy <- xy[j]  - sum(bYs * sol[, 1])
        sryy <- (yy - sum(ys * ys)) - sum(bYs * sol[, 2])
        b    <- srxy / srxx
        rss  <- max(sryy - b * srxy, 0)
        se   <- sqrt((rss / dfi) / srxx)
        tt   <- b / se
        coef.mat[bad[j], ] <- c(b, se, tt, 2 * stats::pt(-abs(tt), df = dfi))
      }
    }
  }

  colnames(coef.mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  return(coef.mat)
}


#' Build the Per-Individual PGS Contribution Matrix
#'
#' Multiplies each hold-out individual's genotype dosages by the per-SNP PGS
#' weights to produce a matrix of weighted dosage contributions used in the
#' forward scan.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants).  \code{NA} values are replaced with \code{0}
#'   in-place before multiplication.
#' @param this.leftout Integer vector. Row indices of the hold-out split in
#'   \code{geno}.
#' @param pgs.weights Numeric vector of length \eqn{N}. Per-SNP effect-size
#'   weights, typically the \code{Estimate} column from
#'   \code{\link{compute_gwas_weights_original}}.
#'
#' @details
#' Missing genotype values (\code{NA}) in \code{geno} are imputed to zero
#' before computing the product.  This is a simple mean-imputation equivalent
#' under a centred dosage scale and is sufficient for the forward scan, where
#' the primary goal is to aggregate regional signals rather than obtain
#' individual-level accuracy.
#'
#' The caller's \code{geno} object is \strong{not} modified: the
#' \code{NA}-replacement assignment triggers R's copy-on-modify semantics, so a
#' local copy of \code{geno} is made inside the function.  Callers working with
#' very large genotype matrices should be aware that this temporarily doubles
#' the memory footprint of \code{geno} during the call; use
#' \code{\link{ras_memory}} to check feasibility before running.
#'
#' @return
#' A numeric matrix of dimensions
#' \eqn{|\texttt{this.leftout}| \times N}, where entry
#' \eqn{[i, j]} is the weighted dosage contribution of SNP \eqn{j} for
#' hold-out individual \eqn{i} (i.e., \code{geno[this.leftout[i], j] *
#' pgs.weights[j]}).
#'
#' @seealso
#' \code{\link{compute_gwas_weights_original}} for the step that produces
#' \code{pgs.weights}.
#' \code{\link{screen_forward_max_region_original}} for the step that consumes this
#' matrix.
#' \code{\link{ras_memory}} for pre-flight memory estimation.
#'
#' @examples
#' set.seed(2)
#' geno    <- matrix(sample(0:2, 50 * 20, replace = TRUE), nrow = 50, ncol = 20)
#' weights <- rnorm(20)
#' leftout <- 1:15
#'
#' pgs_mat <- compute_pgs_matrix(geno, leftout, weights)
#' dim(pgs_mat)   # 15 x 20
#'
#' ## Entry [i, j] equals geno[leftout[i], j] * weights[j]
#' stopifnot(pgs_mat[1, 1] == geno[leftout[1], 1] * weights[1])
#' @export
compute_pgs_matrix <- function(geno, this.leftout, pgs.weights) {
  geno[is.na(geno)] <- 0L
  # Each hold-out row is its dosage vector scaled column-wise by pgs.weights;
  # sweep does this in one vectorized pass (identical to the per-row product).
  pgs.mat <- sweep(geno[this.leftout, , drop = FALSE], 2L, pgs.weights, "*")
  return(pgs.mat)
}
