# One train/holdout repetition of the Stage-1 scan: draw a random 50/50
# split, fit GWAS weights on the training half (compute_gwas_weights),
# then scan the holdout half (screen_forward_max_region). Factored out
# of ras_scan() so the exact same rep logic can run either in a serial
# for-loop (cores = 1) or dispatched to a parallel::parLapply() worker
# (cores > 1) -- see the `cores` parameter docs on ras_scan() for why
# this loop was worth parallelizing and what changes when it is.
.ras_scan_one_rep <- function(this.rep, rasbin_path, phenotype, cov.df,
                               covariate_cols, is_continuous, gwas_formula,
                               scan_formula, rows, chrom, save_dir, skip1, skip2,
                               min_window_size, max_window_size, chunk_snps,
                               verbose) {
  this.sample  <- sort(sample(rows, length(rows) %/% 2))
  this.leftout <- setdiff(rows, this.sample)

  if (is_continuous) {
    train.df            <- cov.df[this.sample, , drop = FALSE]
    train.df$phenotype1 <- phenotype[this.sample]
    lm0 <- lm(
      as.formula(paste("phenotype1 ~", paste(covariate_cols, collapse = " + "))),
      data = train.df
    )
    phenotype1 <- as.numeric(lm0$residuals)
  } else {
    phenotype1 <- phenotype[this.sample]
  }
  phenotype2 <- phenotype[this.leftout]

  coef.mat <- compute_gwas_weights(
    rasbin_path       = rasbin_path,
    phenotype1        = phenotype1,
    this.sample       = this.sample,
    this.df           = cov.df,
    is_continuous     = is_continuous,
    covariate_formula = gwas_formula,
    chunk_snps        = chunk_snps
  )
  pgs.weights <- coef.mat[, 1]
  if (!is.null(save_dir))
    saveRDS(coef.mat, file.path(save_dir,
      paste0("chr-", chrom, "_coef_mat-", this.rep, ".rds")))

  scan.df            <- cov.df[this.leftout, , drop = FALSE]
  scan.df$phenotype2 <- phenotype2

  return.p.values <- screen_forward_max_region(
    rasbin_path       = rasbin_path,
    weights           = pgs.weights,
    this.leftout      = this.leftout,
    this.df           = scan.df,
    is_continuous     = is_continuous,
    covariate_formula = scan_formula,
    skip1             = skip1,
    skip2             = skip2,
    min_window_size   = min_window_size,
    max_window_size   = max_window_size,
    chunk_snps        = chunk_snps
  )

  rm(coef.mat, pgs.weights, scan.df)
  release_memory(verbose = verbose)

  return.p.values
}

#' RAS Stage 1: Averaged Regional Association Profile
#'
#' Computes the \eqn{-\log_{10}(p)} profile that the RAS detectors work on.
#' The sample is split at random into a training half and a hold-out half
#' \code{num_rep} times; in each repetition per-SNP regression weights are
#' fitted on the training half (\code{\link{compute_gwas_weights}}) and a
#' forward regional scan is run on the hold-out half
#' (\code{\link{screen_forward_max_region}}); the profiles are averaged over
#' repetitions. Genotypes are streamed from a \code{.rasbin} file in chunks of
#' \code{chunk_snps} SNP columns, so peak memory does not grow with the
#' chromosome. For a genotype matrix small enough to sit in memory the
#' pure-R implementation \code{\link{ras_scan_original}} gives the same
#' profile.
#'
#' @param geno Either the path to a \code{.rasbin} genotype file (written by
#'   \code{\link{geno_to_rasbin}} or \code{\link{bed_to_rasbin}}), or an
#'   in-memory numeric genotype matrix (\eqn{n} samples by \eqn{N} variants).
#'   A matrix is converted to a temporary \code{.rasbin} file for the run;
#'   convert once with \code{\link{geno_to_rasbin}} when running repeatedly.
#' @param phenotype Numeric vector of length \eqn{n}, in the sample order of
#'   \code{geno}.
#' @param covariates Data frame with \eqn{n} rows, in the sample order of
#'   \code{geno}.
#' @param covariate_cols Character vector. Names of the columns of
#'   \code{covariates} to adjust for.
#' @param is_continuous Logical. \code{TRUE} for a quantitative trait,
#'   \code{FALSE} for a binary (0/1) trait.
#' @param num_rep Integer. Number of train/hold-out repetitions to average
#'   over. Default \code{5}.
#' @param skip1 Integer. Stride of the profile grid in SNPs: the profile has
#'   one value every \code{skip1} SNPs. Default \code{10}.
#' @param skip2 Integer. Step, in SNPs, by which the scan window grows from
#'   \code{min_window_size} to \code{max_window_size} at each grid position.
#'   Default \code{20}.
#' @param min_window_size,max_window_size Integer. Smallest and largest scan
#'   window, in SNPs. Default \code{5} and \code{100}.
#' @param chrom Integer or character. Chromosome label used in the output file
#'   names. Default \code{1}.
#' @param save_dir Character or \code{NULL}. Directory that receives the
#'   per-repetition coefficient matrices and the averaged profile; \code{NULL}
#'   writes nothing. Default \code{file.path(tempdir(), "RAS")}.
#' @param scan_test Character. Per-window test for a binary trait:
#'   \code{"score"} (default, Rao score test in closed form) or \code{"glm"}
#'   (per-window logistic regression, the default of RAS 1.0.x).
#'   \code{"glm"} runs the in-memory implementation
#'   \code{\link{ras_scan_original}} and therefore needs a genotype matrix,
#'   not a \code{.rasbin} path. Ignored for continuous traits.
#' @param chunk_snps Integer. SNP columns held in memory at a time. Bounds
#'   peak memory at roughly \code{n * chunk_snps * 8} bytes per worker.
#'   Default \code{5000}.
#' @param cores Integer. Number of worker processes the repetitions are spread
#'   over (a PSOCK cluster). Default \code{1}. Each worker draws its own
#'   random-number stream, so results with \code{cores > 1} are not
#'   bit-identical to a serial run with the same seed, and peak memory grows
#'   roughly in proportion to \code{cores}.
#' @param keep_reps Logical. Also return (as \code{$reps}) and save the
#'   individual per-repetition profiles, which are needed to study the
#'   split-to-split variability that the average removes. Default
#'   \code{FALSE}.
#' @param rows Integer vector or \code{NULL}. Row indices of the samples to
#'   use; the train/hold-out splits are drawn within this set. \code{NULL}
#'   (default) uses all samples. Used by \code{\link{ras}} to calibrate the
#'   box-scan threshold on a random subset of a large cohort.
#'
#' @return Invisibly, a list with \code{x} (the SNP index of each grid point,
#'   \code{seq(1, N, by = skip1)}), \code{y} (the averaged
#'   \eqn{-\log_{10}(p)} profile) and \code{reps} (a \code{length(x)} by
#'   \code{num_rep} matrix of per-repetition profiles when
#'   \code{keep_reps = TRUE}, otherwise \code{NULL}). The averaged profile is
#'   also written to \code{save_dir} as
#'   \code{mean_p_values_chr<chrom>_reps1-<num_rep>.rds}.
#'
#' @seealso \code{\link{ras}} for the full pipeline; \code{\link{ras_detect}}
#'   and \code{\link{ras_box_detect}} for the detectors that consume the
#'   profile; \code{\link{ras_scan_external}} for a scan with external
#'   weights; \code{\link{ras_scan_original}} for the pure-R implementation.
#' @export
ras_scan <- function(geno, phenotype, covariates, covariate_cols,
                           is_continuous,
                           num_rep         = 5,
                           skip1           = 10,
                           skip2           = 20,
                           chrom           = 1,
                           save_dir        = file.path(tempdir(), "RAS"),
                           min_window_size = 5,
                           max_window_size = 100,
                           scan_test       = c("score", "glm"),
                           chunk_snps      = 5000,
                           cores           = 1,
                           keep_reps       = FALSE,
                           rows            = NULL) {

  scan_test <- match.arg(scan_test)
  if (scan_test == "glm" && !isTRUE(is_continuous)) {
    if (is.character(geno))
      stop("scan_test = \"glm\" is only available on the in-memory route: ",
           "call ras_scan_original() with a genotype matrix, or use scan_test = \"score\".",
           call. = FALSE)
    message("scan_test = \"glm\": running the in-memory reference scan ras_scan_original().")
    out <- ras_scan_original(
      geno = geno, phenotype = phenotype, covariates = covariates,
      covariate_cols = covariate_cols, is_continuous = is_continuous,
      num_rep = num_rep, skip1 = skip1, skip2 = skip2, chrom = chrom,
      save_dir = save_dir, min_window_size = min_window_size,
      max_window_size = max_window_size, scan_test = "glm")
    out$reps <- NULL
    return(invisible(out))
  }
  rb <- .ras_resolve_geno(geno)
  on.exit(rb$cleanup(), add = TRUE)
  rasbin_path <- rb$path

  # Force every argument that .ras_scan_one_rep() reads but that isn't
  # otherwise touched before the cores > 1 branch below. R arguments are
  # lazy promises tied to the *caller's* environment (e.g. a caller writing
  # `ras_scan(..., phenotype = y, ...)` leaves `phenotype` as an
  # unevaluated reference to the symbol `y` in the caller's frame). The
  # cores > 1 path hands a closure over these arguments to
  # parallel::parLapply(), which serializes it to a freshly-spawned worker
  # process -- if any argument is still an unforced promise at that point,
  # the worker inherits a promise pointing at an environment that doesn't
  # exist there, and forcing it on first use fails with "object 'y' not
  # found" instead of quietly evaluating to the caller's value the way it
  # would in the serial (cores = 1) path. Forcing here, in the master
  # process, avoids that regardless of which arguments a given call site
  # passes as literals vs. variables.
  force(phenotype); force(is_continuous); force(chrom)
  force(skip1); force(skip2); force(min_window_size); force(max_window_size)
  force(chunk_snps)

  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  hdr <- rasbin_header(rasbin_path)
  n <- as.integer(hdr["n_samples"])
  N <- as.integer(hdr["n_snps"])
  rows <- if (is.null(rows)) seq_len(n) else sort(unique(as.integer(rows)))
  if (length(rows) < 4L || any(rows < 1L | rows > n))
    stop("rows must index at least 4 samples within 1..", n, call. = FALSE)
  force(rows)

  cov.df <- covariates[, covariate_cols, drop = FALSE]

  gwas_formula <- paste("this.x +", paste(covariate_cols, collapse = " + "))
  scan_formula <- paste(covariate_cols, collapse = " + ")

  this.seq <- seq(1, N, by = skip1)
  cores    <- max(1L, as.integer(cores))

  if (cores <= 1L) {
    rep_p_values <- vector("list", num_rep)
    for (this.rep in seq_len(num_rep)) {
      message("============================== Repetition ", this.rep, " / ", num_rep)
      rep_p_values[[this.rep]] <- .ras_scan_one_rep(
        this.rep        = this.rep,
        rasbin_path     = rasbin_path,
        phenotype       = phenotype,
        cov.df          = cov.df,
        covariate_cols  = covariate_cols,
        is_continuous   = is_continuous,
        gwas_formula    = gwas_formula,
        scan_formula    = scan_formula,
        rows            = rows,
        chrom           = chrom,
        save_dir        = save_dir,
        skip1           = skip1,
        skip2           = skip2,
        min_window_size = min_window_size,
        max_window_size = max_window_size,
        chunk_snps      = chunk_snps,
        verbose         = TRUE
      )
    }
  } else {
    message("============================== Running ", num_rep,
            " repetition(s) across ", cores, " parallel worker(s)")
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressMessages(library(RAS)))
    parallel::clusterSetRNGStream(cl)
    rep_p_values <- parallel::parLapply(cl, seq_len(num_rep), function(this.rep) {
      .ras_scan_one_rep(
        this.rep        = this.rep,
        rasbin_path     = rasbin_path,
        phenotype       = phenotype,
        cov.df          = cov.df,
        covariate_cols  = covariate_cols,
        is_continuous   = is_continuous,
        gwas_formula    = gwas_formula,
        scan_formula    = scan_formula,
        rows            = rows,
        chrom           = chrom,
        save_dir        = save_dir,
        skip1           = skip1,
        skip2           = skip2,
        min_window_size = min_window_size,
        max_window_size = max_window_size,
        chunk_snps      = chunk_snps,
        verbose         = FALSE
      )
    })
  }

  full.p.values <- Reduce(`+`, rep_p_values)
  mean.p.values <- full.p.values / num_rep
  if (!is.null(save_dir))
    saveRDS(mean.p.values,
      file.path(save_dir,
        paste0("mean_p_values_chr", chrom, "_reps1-", num_rep, ".rds")))

  # Averaging is exactly the step that destroys the split-to-split variability,
  # so anything measuring that variability -- pairwise profile correlation,
  # pointwise SD, boundary dispersion, region-recovery frequency -- needs the
  # per-repetition profiles rather than their mean. They are already in memory
  # here (`rep_p_values`, built above and consumed by the Reduce), so keeping
  # them costs no extra computation and no extra peak memory: it only declines
  # to throw them away. Off by default so the returned object and the files on
  # disk stay exactly what callers already expect.
  reps <- NULL
  if (isTRUE(keep_reps)) {
    reps <- do.call(cbind, rep_p_values)
    colnames(reps) <- paste0("rep", seq_len(num_rep))
    if (!is.null(save_dir))
      saveRDS(reps,
        file.path(save_dir,
          paste0("rep_p_values_chr", chrom, "_reps1-", num_rep, ".rds")))
    message(sprintf("Kept %d per-repetition profiles (%d grid points each).",
                    num_rep, nrow(reps)))
  }
  message("Scanning complete!")

  invisible(list(x = this.seq, y = mean.p.values, reps = reps))
}


#' RAS: Regional Association Score Analysis
#'
#' One-call entry point of the package. Runs the RAS pipeline on one
#' chromosome: the Stage-1 scan (\code{\link{ras_scan}}), region detection on
#' the resulting profile with the box-scan detector
#' (\code{\link{ras_box_detect}}, default) or the changepoint detector
#' (\code{\link{ras_detect}} followed by \code{\link{ras_validate}}), and the
#' diagnostic plots (\code{\link{plot.ras}}). Genotypes are streamed from a
#' \code{.rasbin} file; an in-memory genotype matrix is accepted as well.
#'
#' The box-scan threshold is calibrated on the data: unless
#' \code{box_calibration} or \code{box_threshold} is supplied, the scan is
#' repeated on \code{box_null} permuted phenotypes and the threshold is the
#' order statistic of their maximum scores at level \code{box_alpha}
#' (\code{\link{ras_box_calibrate}}), which controls the family-wise error
#' rate without distributional assumptions. This multiplies the run time by
#' about \code{box_null}; for large cohorts calibrate once with
#' \code{\link{ras_box_calibrate}} and pass the result.
#'
#' @inheritParams ras_scan
#' @param cp_p_threshold Numeric. Davies test p-value threshold for a
#'   first-pass candidate changepoint. Default \code{0.01}.
#' @param cp_window_size Integer. Width, in grid points, of the sliding window
#'   of the first pass. Default \code{3000}.
#' @param cp_min_length Integer. Minimum number of grid points on each side of
#'   a candidate. Default \code{10}.
#' @param cp_slope_check_window Integer. Half-width, in grid points, of the
#'   local window in which the slopes on each side of a candidate are tested.
#'   Default \code{30}.
#' @param cp_slope_left,cp_slope_right Numeric. One-tailed p-value thresholds
#'   for the rising left slope and the falling right slope. Default
#'   \code{1e-10} and \code{1e-20}.
#' @param second_window_size Integer. Half-width, in grid points, of the
#'   second-pass local Davies tests. Default \code{50}.
#' @param second_p_threshold Numeric. Davies test p-value threshold of the
#'   second pass. Default \code{1e-10}.
#' @param min_signal Numeric. Minimum profile value \eqn{-\log_{10}(p)} at a
#'   changepoint for it to be kept; also the colour boundary in the plots.
#'   Default \code{2.5}.
#' @param detector Character. \code{"box"} (default): the box-scan detector
#'   \code{\link{ras_box_detect}}, which reports intervals
#'   \eqn{[\tau_L, \tau_R]} and also delimits broad plateau-shaped regions.
#'   \code{"changepoint"}: the two-pass changepoint detector of RAS 1.0.x,
#'   \code{\link{ras_detect}} then \code{\link{ras_validate}}, which reports
#'   single positions. The \code{box_*} arguments belong to \code{"box"}, the
#'   \code{cp_*}, \code{second_*} and \code{min_signal} arguments to
#'   \code{"changepoint"}.
#' @param box_threshold Numeric or \code{NULL}. Fixed score threshold for
#'   \code{detector = "box"}. \code{NULL} (default) uses
#'   \code{box_calibration}, or calibrates on \code{box_null} permutations.
#' @param box_calibration A \code{\link{ras_box_calibrate}} result or
#'   \code{NULL} (default). Supplies the threshold and \code{level0} for
#'   \code{detector = "box"} and skips the permutation calibration.
#' @param box_null Integer. Number of permuted-phenotype scans used to
#'   calibrate the box-scan threshold when neither \code{box_threshold} nor
#'   \code{box_calibration} is given. Default \code{100}. Continuous
#'   phenotypes are permuted as residuals after the covariates
#'   (Freedman-Lane); binary phenotypes are permuted directly. \code{0}
#'   disables the calibration, in which case a threshold must be supplied.
#' @param box_alpha Numeric. Family-wise error rate of the calibrated
#'   threshold. Default \code{0.05}.
#' @param box_null_n Integer or \code{NULL}. Number of samples, drawn at
#'   random once, on which the permuted-phenotype scans are run. \code{NULL}
#'   (default) uses all samples. The null distribution of the profile maximum
#'   depends on the linkage structure, not on the sample size, so a subset of
#'   a few tens of thousands of individuals calibrates a biobank-scale cohort
#'   at a fraction of the cost (on a 352-pig chromosome, subsets of a third
#'   and a sixth of the animals gave thresholds within 3 percent of the
#'   full-sample one).
#' @param box_method Character. How the calibrated threshold is read off the
#'   null maxima: \code{"order"} (default, exact order statistic, use with
#'   \code{box_null} of 100 or more) or \code{"gumbel"} (Gumbel fit, usable
#'   from \code{box_null = 20} with a slightly liberal threshold). See
#'   \code{\link{ras_box_calibrate}}.
#' @param box_level0 Numeric. Null-profile median for the box scan's level
#'   factor (see \code{\link{ras_box_stat}}). \code{NA} (default) takes it
#'   from \code{box_calibration} when given and otherwise disables the factor.
#' @param run_plots Logical. Save the diagnostic plots to \code{save_dir}.
#'   Default \code{TRUE}.
#' @param plot_device Character. \code{"pdf"} (default), \code{"png"} or
#'   \code{"screen"}.
#' @param plot_p_threshold Numeric. Significance reference line drawn on the
#'   plots, on the \eqn{-\log_{10}} scale. Default \code{8}.
#' @param plot_y_cap Numeric or \code{NULL}. Cap of the plotted y-axis.
#'   Default \code{NULL}.
#'
#' @details
#' The pure-R, in-memory implementation that RAS 1.0.x shipped as
#' \code{ras()} is available as \code{\link{ras_original}}; it defaults to
#' the changepoint detector, and with \code{detector = "changepoint"} the two
#' functions give the same scan profile to machine precision and the same
#' detected positions.
#'
#' @return Invisibly, an object of class \code{"ras"}: a list with
#'   \code{scan} (the \code{\link{ras_scan}} result), \code{detection} (for
#'   \code{detector = "box"} the \code{\link{ras_box_detect}} result, whose
#'   \code{regions} table has one row per interval; for
#'   \code{detector = "changepoint"} the \code{\link{ras_validate}} result,
#'   whose \code{tau_hats} are the detected positions),
#'   \code{box_calibration} (the calibration used, or \code{NULL}),
#'   \code{chrom} and \code{save_dir}. Use \code{print()} and \code{plot()}
#'   on it.
#'
#' @seealso \code{\link{ras_scan}}, \code{\link{ras_detect}},
#'   \code{\link{ras_validate}}, \code{\link{ras_box_detect}},
#'   \code{\link{plot.ras}}; \code{\link{geno_to_rasbin}} and
#'   \code{\link{bed_to_rasbin}} to prepare the genotype file;
#'   \code{\link{ras_original}} for the pure-R implementation.
#'
#' @examples
#' \donttest{
#' set.seed(3)
#' n_samp <- 120; n_snp <- 400
#' geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE,
#'                       prob = c(0.6, 0.3, 0.1)), n_samp, n_snp)
#' causal <- 181:220                       # 40 causal SNPs, together explaining
#' pheno  <- 2 * as.numeric(scale(rowSums(geno[, causal]))) + rnorm(n_samp)  # ~80% of the variance
#' cov_df <- data.frame(age = rnorm(n_samp), sex = rbinom(n_samp, 1, 0.5))
#'
#' ## default: box-scan detector, threshold calibrated on permuted phenotypes
#' ## (20 permutations with the Gumbel fit keep the example short; the
#' ## default is the exact order statistic on 100 permutations)
#' res <- ras(geno, pheno, cov_df, covariate_cols = c("age", "sex"),
#'            is_continuous = TRUE, num_rep = 2, skip1 = 2, skip2 = 5,
#'            min_window_size = 2, max_window_size = 20,
#'            box_null = 20, box_method = "gumbel",
#'            save_dir = tempdir(), run_plots = FALSE)
#' res$detection$regions
#' res$box_calibration
#'
#' ## large data: convert once and pass the file; here with the changepoint
#' ## detector of RAS 1.0.x (settings scaled down to the 200-point profile)
#' rb <- tempfile(fileext = ".rasbin")
#' geno_to_rasbin(geno, rb)
#' res_cp <- ras(rb, pheno, cov_df, covariate_cols = c("age", "sex"),
#'               is_continuous = TRUE, num_rep = 2, skip1 = 2, skip2 = 5,
#'               min_window_size = 2, max_window_size = 20,
#'               detector = "changepoint",
#'               cp_window_size = 100, cp_slope_check_window = 10,
#'               cp_slope_left = 1e-2, cp_slope_right = 1e-2,
#'               second_window_size = 20, second_p_threshold = 1e-2,
#'               save_dir = tempdir(), run_plots = FALSE)
#' res_cp$detection$tau_hats
#' unlink(c(rb, paste0(rb, ".meta.rds")))
#' }
#' @export
ras <- function(geno, phenotype, covariates, covariate_cols,
                      is_continuous,
                      num_rep               = 5,
                      skip1                 = 10,
                      skip2                 = 20,
                      chrom                 = 1,
                      save_dir              = file.path(tempdir(), "RAS"),
                      min_window_size       = 5,
                      max_window_size       = 100,
                      scan_test             = c("score", "glm"),
                      chunk_snps            = 5000,
                      cores                 = 1,
                      keep_reps             = FALSE,
                      cp_p_threshold        = 0.01,
                      cp_window_size        = 3000,
                      cp_min_length         = 10,
                      cp_slope_check_window = 30,
                      cp_slope_left         = 1e-10,
                      cp_slope_right        = 1e-20,
                      second_window_size    = 50,
                      second_p_threshold    = 1e-10,
                      min_signal            = 2.5,
                      detector              = c("box", "changepoint"),
                      box_threshold         = NULL,
                      box_calibration       = NULL,
                      box_level0            = NA_real_,
                      box_null              = 100,
                      box_alpha             = 0.05,
                      box_null_n            = NULL,
                      box_method            = c("order", "gumbel"),
                      run_plots             = TRUE,
                      plot_device           = "pdf",
                      plot_p_threshold      = 8,
                      plot_y_cap            = NULL) {

  detector   <- match.arg(detector)
  scan_test  <- match.arg(scan_test)
  box_method <- match.arg(box_method)

  # Legacy binary glm scan: only the in-memory reference pipeline offers it.
  if (scan_test == "glm" && !isTRUE(is_continuous)) {
    if (is.character(geno))
      stop("scan_test = \"glm\" is only available on the in-memory route: ",
           "call ras_original() with a genotype matrix, or use scan_test = \"score\".",
           call. = FALSE)
    message("scan_test = \"glm\": running the in-memory reference pipeline ras_original().")
    return(invisible(ras_original(
      geno = geno, phenotype = phenotype, covariates = covariates,
      covariate_cols = covariate_cols, is_continuous = is_continuous,
      num_rep = num_rep, skip1 = skip1, skip2 = skip2, chrom = chrom,
      save_dir = save_dir, min_window_size = min_window_size,
      max_window_size = max_window_size, scan_test = "glm",
      cp_p_threshold = cp_p_threshold, cp_window_size = cp_window_size,
      cp_min_length = cp_min_length, cp_slope_check_window = cp_slope_check_window,
      cp_slope_left = cp_slope_left, cp_slope_right = cp_slope_right,
      second_window_size = second_window_size, second_p_threshold = second_p_threshold,
      min_signal = min_signal, detector = detector, box_threshold = box_threshold,
      box_calibration = box_calibration, box_level0 = box_level0,
      run_plots = run_plots, plot_device = plot_device,
      plot_p_threshold = plot_p_threshold, plot_y_cap = plot_y_cap)))
  }

  # Convert a matrix once here so that ras_scan() does not do it again.
  rb <- .ras_resolve_geno(geno)
  on.exit(rb$cleanup(), add = TRUE)

  # ── Stage 1: scan ─────────────────────────────────────────────────────────
  scan <- ras_scan(
    geno            = rb$path,
    phenotype       = phenotype,
    covariates      = covariates,
    covariate_cols  = covariate_cols,
    is_continuous   = is_continuous,
    num_rep         = num_rep,
    skip1           = skip1,
    skip2           = skip2,
    chrom           = chrom,
    save_dir        = save_dir,
    min_window_size = min_window_size,
    max_window_size = max_window_size,
    scan_test       = scan_test,
    chunk_snps      = chunk_snps,
    cores           = cores,
    keep_reps       = keep_reps
  )

  x <- scan$x
  y <- scan$y

  # -- Box-scan threshold: calibrate on permuted phenotypes unless supplied --
  if (detector == "box" && is.null(box_threshold) && is.null(box_calibration) &&
      box_null >= 2) {
    box_calibration <- .ras_box_autocalibrate(
      rasbin_path = rb$path, phenotype = phenotype, covariates = covariates,
      covariate_cols = covariate_cols, is_continuous = is_continuous,
      box_null = box_null, box_alpha = box_alpha, box_method = box_method,
      box_null_n = box_null_n, cores = cores,
      scan_args = list(num_rep = num_rep, skip1 = skip1, skip2 = skip2,
                       chrom = chrom, min_window_size = min_window_size,
                       max_window_size = max_window_size,
                       scan_test = scan_test, chunk_snps = chunk_snps))
  }

  # -- Stage 2 + 3: region detection on the profile -------------------------
  detection <- .ras_run_detection(
    x = x, y = y, detector = detector, fast = TRUE, skip1 = skip1,
    cp_p_threshold = cp_p_threshold, cp_window_size = cp_window_size,
    cp_min_length = cp_min_length, cp_slope_check_window = cp_slope_check_window,
    cp_slope_left = cp_slope_left, cp_slope_right = cp_slope_right,
    second_window_size = second_window_size, second_p_threshold = second_p_threshold,
    min_signal = min_signal, box_threshold = box_threshold,
    box_calibration = box_calibration, box_level0 = box_level0)

  result <- structure(
    list(scan = scan, detection = detection,
         box_calibration = if (detector == "box") box_calibration else NULL,
         chrom = chrom, save_dir = save_dir),
    class = "ras"
  )

  # ── Stage 4: plots (unchanged, already cheap) ───────────────────────────────
  if (run_plots) {
    message(" Generating plots...")
    plot.ras(result, zoom = FALSE, device = plot_device,
             p.threshold = plot_p_threshold, y_cap = plot_y_cap,
             min_signal  = min_signal)
    plot.ras(result, zoom = TRUE,  device = plot_device,
             p.threshold = plot_p_threshold, min_signal = min_signal)
    message(sprintf(" Plots saved to: %s", save_dir))
  }

  message("============================== Pipeline complete.")

  invisible(result)
}
