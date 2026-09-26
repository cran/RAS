# Stage 2/3 dispatch shared by ras_original() and ras(): run whichever region
# detector the caller asked for on the scan profile and return the
# `detection` element of the "ras" object.
#
#   detector = "changepoint"  the original two-pass detector: ras_detect_original()
#                             (or its C port ras_detect() when fast =
#                             TRUE) followed by ras_validate().
#   detector = "box"          the box-scan region detector ras_box_detect(),
#                             which needs a threshold (box_threshold) or a
#                             calibration object (box_calibration).
.ras_run_detection <- function(x, y, detector, fast, skip1,
                               cp_p_threshold, cp_window_size, cp_min_length,
                               cp_slope_check_window, cp_slope_left, cp_slope_right,
                               second_window_size, second_p_threshold, min_signal,
                               box_threshold, box_calibration, box_level0) {
  if (detector == "changepoint") {
    message("============================== Changepoint Detection (Pass 1)")
    detect_fun <- if (isTRUE(fast)) ras_detect else ras_detect_original
    cp_result <- detect_fun(
      x                              = x,
      y                              = y,
      p.values.threshold             = cp_p_threshold,
      min.length                     = cp_min_length,
      window_size                    = cp_window_size,
      slope_check_window_size        = cp_slope_check_window,
      slope.p.values.threshold.left  = cp_slope_left,
      slope.p.values.threshold.right = cp_slope_right
    )

    message("============================== Changepoint Validation (Pass 2)")
    detection <- ras_validate(
      this.result        = cp_result,
      x                  = x,
      y                  = y,
      this.start         = 1,
      this.skip          = skip1,
      second_window_size = second_window_size,
      p.value.threshold  = second_p_threshold,
      min_signal         = min_signal
    )
    detection$detector <- "changepoint"

    n_detected <- length(detection$tau_hats)
    if (n_detected == 0) {
      message("============================== No changepoints detected.")
    } else {
      message(sprintf(
        "============================== Detected %d changepoint(s) at position(s): %s",
        n_detected, paste(detection$tau_hats, collapse = ", ")))
    }
    return(detection)
  }

  # detector == "box"
  if (is.null(box_threshold) && is.null(box_calibration))
    stop("detector = \"box\" needs a threshold: pass `box_threshold` (numeric), ",
         "`box_calibration` (a ras_box_calibrate() result), or leave box_null > 0 ",
         "so that ras() calibrates on permuted phenotypes. See ?ras_box_calibrate.",
         call. = FALSE)
  message("============================== Box-Scan Region Detection")
  detection <- ras_box_detect(
    x           = x,
    y           = y,
    threshold   = box_threshold,
    calibration = box_calibration,
    level0      = if (is.null(box_calibration) || !is.na(box_level0)) box_level0 else NULL
  )
  n_detected <- if (is.null(detection$regions)) 0L else nrow(detection$regions)
  if (n_detected == 0) {
    message(sprintf("============================== No regions reach T >= %.3f.",
                    detection$threshold))
  } else {
    r <- detection$regions
    message(sprintf(
      "============================== Detected %d region(s) (T >= %.3f): %s",
      n_detected, detection$threshold,
      paste(sprintf("[%s, %s] T=%.2f", format(r$pos_L), format(r$pos_R), r$T_box),
            collapse = "; ")))
  }
  detection
}

# Calibrate the box-scan threshold for ras() on permuted phenotypes: the same
# Stage-1 scan (same settings, same cohort) is run on `box_null` null
# phenotypes and ras_box_calibrate() turns the maxima into the threshold.
#
# Null phenotypes keep the covariate structure. Continuous: Freedman-Lane, the
# residuals of phenotype ~ covariates are permuted and the fitted values added
# back. Binary: the 0/1 outcome is permuted among the non-missing samples.
.ras_box_autocalibrate <- function(rasbin_path, phenotype, covariates, covariate_cols,
                                   is_continuous, box_null, box_alpha, box_method, box_null_n,
                                   cores, scan_args) {
  n <- length(phenotype)
  rows <- NULL
  if (!is.null(box_null_n) && box_null_n < n) {
    rows <- sort(sample.int(n, box_null_n))
    message(sprintf(
      "============================== Box-Scan Calibration: %d permuted-phenotype scans on %s of %s samples",
      box_null, format(box_null_n, big.mark = ","), format(n, big.mark = ",")))
  } else {
    message(sprintf(
      "============================== Box-Scan Calibration: %d permuted-phenotype scans", box_null))
  }
  message("  (each one repeats the Stage-1 scan; supply box_calibration or box_threshold to skip)")
  ok <- !is.na(phenotype)
  if (isTRUE(is_continuous)) {
    df <- covariates[, covariate_cols, drop = FALSE]
    df$.y <- phenotype
    fit <- stats::lm(.y ~ ., data = df, na.action = stats::na.exclude)
    fitted_y <- as.numeric(stats::fitted(fit))
    resid_y  <- as.numeric(stats::residuals(fit))
    ok <- ok & !is.na(resid_y)
    draw <- function() { y <- phenotype; y[ok] <- fitted_y[ok] + sample(resid_y[ok]); y }
  } else {
    draw <- function() { y <- phenotype; y[ok] <- sample(phenotype[ok]); y }
  }
  nulls <- vector("list", box_null)
  for (b in seq_len(box_null)) {
    nulls[[b]] <- suppressMessages(do.call(ras_scan, c(
      list(geno = rasbin_path, phenotype = draw(), covariates = covariates,
           covariate_cols = covariate_cols, is_continuous = is_continuous,
           save_dir = NULL, cores = cores, rows = rows), scan_args)))$y
    if (b %% 10 == 0 || b == box_null)
      message(sprintf("  null scan %d / %d", b, box_null))
  }
  cal <- ras_box_calibrate(nulls, alpha = box_alpha, method = box_method)
  message(sprintf("  threshold = %.3f (%s, alpha %.3g, %d nulls), level0 = %.3f",
                  cal$threshold, box_method, box_alpha, box_null, cal$level0))
  cal
}
