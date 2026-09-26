# =============================================================================
# Box-scan region detector: a sliding-window contrast in RAW profile units.
#
# WHY A SECOND DETECTOR
#   The original RAS detector (ras_detect_original + ras_validate) looks for one
#   slope-reversal changepoint: background -> rising -> peak -> falling ->
#   background. That shape test structurally cannot place a boundary on a
#   plateau (a broad region whose interior is flat), and its first-pass
#   window can miss a wide signal altogether. The box scan instead scores
#   every window of a set of widths by how far its mean stands above BOTH
#   adjacent flanks and above the profile background, in the raw
#   -log10(p) units of the profile. No noise scale is estimated: on the
#   published benchmark about half of the alternative profiles have a
#   globally raised, wandering background, and every scale estimate tried
#   either over-penalised obvious regions or could not be estimated
#   symmetrically under the null and the alternative. Height in raw units
#   separates background bumps (~2 above their surroundings) from true
#   regions (7-15); height divided by a noise scale does not.
#
# THE STATISTIC
#   For a window of w grid points starting at i, with flanks of
#   f = max(5, round(w / 2)) points on each side:
#       contrast(i, w) = mean(window) - max(mean(left flank), mean(right flank))
#       height(i, w)   = mean(window) - median(profile)
#       T(i, w)        = min(contrast, height) / max(1, median(profile) / level0)^gamma
#   The window must stand above EACH flank, so a window inside a plateau, on
#   the shoulder of a peak, or spanning a cluster of neighbouring peaks scores
#   about 0. level0 is the typical median of a null profile; the level
#   factor raises the bar on a profile whose whole background is raised
#   (-log10 p noise grows with its level). Widths are NOT normalised against
#   each other: dividing each width by its own null maximum favours wide
#   windows and merges neighbouring regions into one.
#
# REGIONS
#   Windows are taken greedily by T; a window is dropped when it touches an
#   accepted one expanded by 10 + max(5, width / 2) grid points (the response
#   skirt plus the flank), which is what removes side lobes.
#
# CALIBRATION
#   The threshold on T is an order statistic of the maximum T over null
#   profiles (ras_box_calibrate). Null profiles come from the same scan run
#   on permuted phenotypes, or on the same genotypes with a simulated null
#   phenotype.
# =============================================================================

.ras_box_default_widths <- function() c(5L, 8L, 11L, 14L, 17L, 20L, 24L, 28L, 34L, 42L, 55L, 75L, 100L)

# flank width used with a window of w grid points
.ras_box_flank <- function(w) pmax(5L, as.integer(round(w / 2)))

#' Box-Scan Statistic for Every Window Start and Width
#'
#' Computes the box-scan score \eqn{T(i, w)} of every window of every width
#' in \code{widths} along a RAS scan profile. This is the building block of
#' \code{\link{ras_box_detect}}; most users will call that function instead.
#'
#' @param y Numeric vector. The scan profile (\eqn{-\log_{10}(p)} values from
#'   \code{\link{ras_scan}} or \code{\link{ras_scan_external}}). Must be
#'   free of \code{NA}.
#' @param widths Integer vector. Window widths in grid points. Default
#'   \code{c(5, 8, 11, 14, 17, 20, 24, 28, 34, 42, 55, 75, 100)}. Widths
#'   that do not fit the profile are dropped.
#' @param level0 Numeric. Typical median of a null profile (element
#'   \code{level0} of a \code{\link{ras_box_calibrate}} result). \code{NA}
#'   (default) disables the level factor.
#' @param gamma Numeric. Exponent of the level factor. Default \code{0.5}.
#' @param edge \code{"open"} (default) or \code{"strict"}. With
#'   \code{"strict"} both flanks must lie inside the profile, so a region at
#'   a chromosome end can never be reported. With \code{"open"} a flank that
#'   runs off the profile is truncated, and dropped altogether when fewer
#'   than 3 points remain; the window then only has to stand above the side
#'   that exists (and above the background).
#'
#' @details
#' For a window of \eqn{w} grid points starting at \eqn{i}, with flanks of
#' \eqn{f = \max(5, w/2)} points on each side,
#' \deqn{contrast = mean(window) - max(mean(left flank), mean(right flank))}
#' \deqn{height = mean(window) - median(y)}
#' \deqn{T = min(contrast, height) / max(1, median(y)/level0)^{gamma}.}
#' All quantities are in the raw units of the profile; no noise scale is
#' estimated. A window inside a plateau, on the shoulder of a peak, or
#' spanning a cluster of neighbouring peaks scores about 0 because at least
#' one flank is as high as the window.
#'
#' @return A data frame with one row per (start, width) pair and columns
#'   \code{i} (window start, grid index), \code{w} (width), \code{contrast},
#'   \code{height} and \code{T}.
#'
#' @seealso \code{\link{ras_box_detect}}, \code{\link{ras_box_calibrate}}.
#'
#' @examples
#' y <- rep(1, 400); y[201:230] <- 8
#' S <- ras_box_stat(y)
#' S[which.max(S$T), ]          # the 30-point window starting at 201
#' @export
ras_box_stat <- function(y, widths = .ras_box_default_widths(), level0 = NA_real_,
                         gamma = 0.5, edge = c("open", "strict")) {
  edge <- match.arg(edge)
  y <- as.numeric(y); n <- length(y)
  if (any(!is.finite(y)))
    stop("ras_box_stat: 'y' contains NA / non-finite values (",
         sum(!is.finite(y)), " of ", n, ")", call. = FALSE)
  widths <- sort(unique(as.integer(widths)))
  if (any(widths < 1L)) stop("ras_box_stat: window widths must be >= 1", call. = FALSE)
  # a window needs both flanks and at least two start positions inside the profile
  fits <- widths + 2L * .ras_box_flank(widths) < n
  if (!any(fits))
    stop("ras_box_stat: a profile of ", n, " points is too short for the smallest window (",
         min(widths), " + 2 x ", .ras_box_flank(min(widths)), " flank points)", call. = FALSE)
  widths <- widths[fits]
  cs <- cumsum(c(0, y)); bg <- stats::median(y)
  lf <- if (is.finite(level0) && gamma > 0) max(1, bg / level0)^gamma else 1
  out <- lapply(widths, function(w) {
    f <- .ras_box_flank(w)
    if (edge == "strict") {
      lo <- f + 1L; hi <- n - w - f + 1L
      if (hi <= lo) return(NULL)
      i <- lo:hi
      win <- (cs[i + w] - cs[i]) / w
      fl <- (cs[i] - cs[i - f]) / f
      fr <- (cs[i + w + f] - cs[i + w]) / f
    } else {
      i <- 1L:(n - w + 1L)
      win <- (cs[i + w] - cs[i]) / w
      a <- pmax(0L, i - 1L - f); nl <- i - 1L - a               # left flank: points a+1 .. i-1
      b <- pmin(n, i + w - 1L + f); nr <- b - (i + w - 1L)      # right flank: points i+w .. b
      fl <- ifelse(nl >= 3L, (cs[i] - cs[a + 1L]) / pmax(nl, 1L), -Inf)
      fr <- ifelse(nr >= 3L, (cs[b + 1L] - cs[i + w]) / pmax(nr, 1L), -Inf)
      ok <- is.finite(fl) | is.finite(fr)
      if (!any(ok)) return(NULL)
      i <- i[ok]; win <- win[ok]; fl <- fl[ok]; fr <- fr[ok]
    }
    contrast <- win - pmax(fl, fr); height <- win - bg
    data.frame(i = i, w = w, contrast = contrast, height = height,
               T = pmin(contrast, height) / lf)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Greedy, side-lobe-suppressed selection of windows from a ras_box_stat() table.
.ras_box_pick <- function(S, n, floor = 0.2, max_regions = Inf) {
  empty <- data.frame(tau_L = integer(0), tau_R = integer(0), width = integer(0),
                      T_box = numeric(0), contrast = numeric(0), height = numeric(0))
  S <- S[is.finite(S$T) & S$T >= floor, , drop = FALSE]
  if (!nrow(S)) return(empty)
  S <- S[order(-S$T), , drop = FALSE]
  L <- S$i; R <- S$i + S$w - 1L
  blocked <- rep(FALSE, n); acc <- integer(0)
  for (m in seq_len(nrow(S))) {
    if (any(blocked[L[m]:R[m]])) next
    marg <- 10L + max(5L, as.integer(round(S$w[m] / 2)))
    blocked[max(1L, L[m] - marg):min(n, R[m] + marg)] <- TRUE
    acc <- c(acc, m)
    if (length(acc) >= max_regions) break
  }
  data.frame(tau_L = L[acc], tau_R = R[acc], width = S$w[acc], T_box = S$T[acc],
             contrast = S$contrast[acc], height = S$height[acc])
}

#' Detect Elevated Regions with the Box Scan
#'
#' The box-scan region detector, the alternative to the changepoint detector
#' (\code{\link{ras_detect}} + \code{\link{ras_validate}}) selected with
#' \code{detector = "box"} in \code{\link{ras}}.
#' It reports intervals \eqn{[\tau_L, \tau_R]} of the scan profile whose
#' mean stands above both adjacent flanks and above the profile background,
#' which lets it delimit broad plateau-shaped association regions as well as
#' sharp peaks.
#'
#' @param x Numeric vector. Grid positions of the profile (\code{scan$x},
#'   the SNP index of each grid point).
#' @param y Numeric vector. The scan profile (\code{scan$y}), same length as
#'   \code{x}.
#' @param threshold Numeric. Keep regions with \code{T_box >= threshold}.
#'   Either \code{threshold} or \code{calibration} must be supplied.
#' @param calibration A \code{\link{ras_box_calibrate}} result. Supplies
#'   \code{threshold} (or \code{threshold_scaled} when \code{scaled = TRUE}),
#'   \code{level0}, \code{widths}, \code{gamma} and \code{edge}; an
#'   explicitly supplied argument overrides the calibrated value.
#' @param scaled Logical. With a \code{calibration} whose
#'   \code{length_ratio > 1}, use its \code{threshold_scaled} (the threshold
#'   extrapolated to a profile that many times longer than the null
#'   profiles) instead of \code{threshold}. Default \code{FALSE}.
#' @param level0,widths,gamma,edge As in \code{\link{ras_box_stat}};
#'   \code{NULL} (default) takes the value from \code{calibration} when one
#'   is given and the \code{\link{ras_box_stat}} default otherwise.
#' @param floor Numeric. Windows with \code{T} below this value are never
#'   considered, whatever the threshold. Default \code{0.2}.
#' @param max_regions Integer. Stop after this many regions. Default
#'   \code{Inf}.
#' @param keep_stat Logical. Also return the full \code{\link{ras_box_stat}}
#'   table as element \code{stat}. Default \code{FALSE} (the table has one row
#'   per window start and width, so it is large for a chromosome-length
#'   profile).
#'
#' @details
#' Windows are scored with \code{\link{ras_box_stat}} and taken greedily by
#' score; a window is dropped when it touches an already accepted one
#' expanded by \code{10 + max(5, width/2)} grid points, which removes side
#' lobes. Each accepted window becomes one region; its anchor is the highest
#' point of the profile inside it.
#'
#' The threshold is a property of the null distribution of the maximum
#' score along a profile and therefore depends on the profile's length and on
#' the scan settings. Obtain it with \code{\link{ras_box_calibrate}} from null
#' profiles produced by the same scan on permuted phenotypes; the order
#' statistic it returns controls the family-wise error rate at the requested
#' level without any distributional assumption.
#'
#' @return A list that \code{\link{plot.ras}} and \code{\link{print.ras}}
#'   accept as a \code{detection} element:
#' \describe{
#'   \item{\code{detector}}{\code{"box"}.}
#'   \item{\code{regions}}{Data frame with one row per region, ordered by
#'     decreasing score: \code{tau_L}, \code{tau_R} (grid indices),
#'     \code{pos_L}, \code{pos_R} (the same boundaries in \code{x} units),
#'     \code{width}, \code{T_box}, \code{contrast}, \code{height},
#'     \code{anchor_idx}, \code{anchor_pos}, \code{anchor_y}. \code{NULL}
#'     when nothing reaches the threshold.}
#'   \item{\code{tau_hats}}{Numeric vector. Anchor positions in \code{x}
#'     units (the same role as in \code{\link{ras_validate}}'s result).}
#'   \item{\code{all.changepoints}, \code{all.p.values}}{Anchor positions and
#'     their \code{T_box} scores, so that the candidate overlay of
#'     \code{\link{plot.ras}} works.}
#'   \item{\code{left.slopes}, \code{right.slopes}}{\code{NULL}; the box scan
#'     estimates no slopes.}
#'   \item{\code{threshold}, \code{level0}, \code{level},
#'     \code{level_factor}}{The threshold used, the null level, the median of
#'     \code{y} and the resulting level factor.}
#'   \item{\code{stat}}{The \code{\link{ras_box_stat}} table when
#'     \code{keep_stat = TRUE}, otherwise \code{NULL}.}
#' }
#'
#' @seealso \code{\link{ras_box_calibrate}} for the threshold,
#'   \code{\link{ras_box_stat}} for the statistic, \code{\link{ras_detect}}
#'   for the changepoint detector, \code{\link{ras}} which calls this
#'   function when \code{detector = "box"}.
#'
#' @examples
#' set.seed(1)
#' x <- seq(1, 6000, by = 10)                 # 600 grid points
#' y <- 0.7 + abs(rnorm(600, sd = 0.3))       # null-like background
#' y[301:342] <- y[301:342] + 8               # a 42-point plateau
#' det <- ras_box_detect(x, y, threshold = 2)
#' det$regions[, c("tau_L", "tau_R", "pos_L", "pos_R", "T_box")]
#'
#' ## a window inside the plateau scores about 0: only its edges are reported
#' S <- ras_box_stat(y)
#' max(S$T[S$w == 20 & S$i >= 305 & S$i + 19 <= 338])
#' @export
ras_box_detect <- function(x, y, threshold = NULL, calibration = NULL, scaled = FALSE,
                           level0 = NULL, widths = NULL, gamma = NULL, edge = NULL,
                           floor = 0.2, max_regions = Inf, keep_stat = FALSE) {
  y <- as.numeric(y); n <- length(y)
  if (length(x) != n)
    stop("ras_box_detect: 'x' and 'y' differ in length (", length(x), " vs ", n, ")", call. = FALSE)
  if (!is.null(calibration)) {
    if (!inherits(calibration, "ras_box_calibration"))
      stop("ras_box_detect: 'calibration' must come from ras_box_calibrate()", call. = FALSE)
    if (is.null(threshold))
      threshold <- if (isTRUE(scaled)) calibration$threshold_scaled else calibration$threshold
    if (is.null(level0)) level0 <- calibration$level0
    if (is.null(widths)) widths <- calibration$widths
    if (is.null(gamma))  gamma  <- calibration$gamma
    if (is.null(edge))   edge   <- calibration$edge
  }
  if (is.null(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("ras_box_detect: supply a numeric 'threshold', or a 'calibration' from ",
         "ras_box_calibrate(); see ?ras_box_calibrate", call. = FALSE)
  if (is.null(level0)) level0 <- NA_real_
  if (is.null(widths)) widths <- .ras_box_default_widths()
  if (is.null(gamma))  gamma  <- 0.5
  if (is.null(edge))   edge   <- "open"

  S <- ras_box_stat(y, widths = widths, level0 = level0, gamma = gamma, edge = edge)
  reg <- .ras_box_pick(S, n, floor = max(floor, threshold), max_regions = max_regions)
  bg <- stats::median(y)
  lf <- if (is.finite(level0) && gamma > 0) max(1, bg / level0)^gamma else 1

  if (nrow(reg)) {
    reg$anchor_idx <- vapply(seq_len(nrow(reg)), function(k)
      reg$tau_L[k] - 1L + which.max(y[reg$tau_L[k]:reg$tau_R[k]]), integer(1))
    reg$anchor_pos <- x[reg$anchor_idx]
    reg$anchor_y   <- y[reg$anchor_idx]
    reg$pos_L <- x[reg$tau_L]
    reg$pos_R <- x[reg$tau_R]
    reg <- reg[, c("tau_L", "tau_R", "pos_L", "pos_R", "width", "T_box", "contrast",
                   "height", "anchor_idx", "anchor_pos", "anchor_y")]
    rownames(reg) <- NULL
  }

  list(detector         = "box",
       regions          = if (nrow(reg)) reg else NULL,
       tau_hats         = if (nrow(reg)) reg$anchor_pos else numeric(0),
       all.changepoints = if (nrow(reg)) reg$anchor_pos else numeric(0),
       all.p.values     = if (nrow(reg)) reg$T_box else numeric(0),
       left.slopes      = NULL,
       right.slopes     = NULL,
       threshold        = threshold,
       level0           = level0,
       level            = bg,
       level_factor     = lf,
       stat             = if (isTRUE(keep_stat)) S else NULL)
}

#' Calibrate the Box-Scan Threshold on Null Profiles
#'
#' Turns a set of null scan profiles into the threshold that
#' \code{\link{ras_box_detect}} needs. Null profiles are obtained by running
#' the same scan (\code{\link{ras_scan}} or
#' \code{\link{ras_scan_external}}, with the same settings) on the same
#' genotypes with the phenotype permuted, or with a phenotype simulated
#' under the null.
#'
#' @param null_profiles A list of numeric vectors, or a matrix with one null
#'   profile per row. All profiles should have the length of the profile
#'   that will be analysed, or the length that \code{length_ratio} refers to.
#' @param alpha Numeric. Target family-wise error rate. Default \code{0.05}.
#' @inheritParams ras_box_stat
#' @param length_ratio Numeric. If the profile to be analysed is
#'   \code{length_ratio} times longer than the null profiles, the null maxima
#'   are fitted by a Gumbel law (method of moments) and the
#'   \eqn{(1-\alpha)^{1/length\_ratio}} quantile is returned as
#'   \code{threshold_scaled}. That number is an extrapolation past the data;
#'   \code{null_exceed_scaled} reports how many null maxima reach it, and a
#'   null of the full length is always preferable. Default \code{1}.
#' @param calib Integer vector. Indices of the profiles used to set the
#'   threshold; the remaining profiles, if any, give a held-out error rate
#'   (\code{fwer_holdout}). Default: all profiles.
#' @param method Character. \code{"order"} (default): the threshold is the
#'   order statistic \eqn{k = \lceil (B+1)(1-\alpha) \rceil} of the \eqn{B}
#'   null maxima, which bounds the family-wise error rate by \eqn{\alpha}
#'   exactly but needs \eqn{B \ge 1/\alpha} profiles and is noisy for small
#'   \eqn{B}. \code{"gumbel"}: the \eqn{1-\alpha} quantile of a Gumbel law
#'   fitted to the null maxima by the method of moments; stable from about 20
#'   profiles on, at the price of a slightly liberal threshold (in a
#'   permutation study on a 224-point profile, 20 to 50 profiles gave 97 to
#'   98 percent of the 200-profile order-statistic threshold and a realised
#'   error rate of about 0.07 for \eqn{\alpha = 0.05}).
#'
#' @details
#' \code{level0} is the median of the null profiles' medians. Both thresholds
#' are always computed and returned (\code{threshold_order},
#' \code{threshold_gumbel}); \code{threshold} is the one selected by
#' \code{method}.
#'
#' @return An object of class \code{"ras_box_calibration"}: a list with
#'   \code{level0}, \code{threshold}, \code{threshold_order},
#'   \code{threshold_gumbel}, \code{threshold_scaled}, \code{method},
#'   \code{length_ratio}, \code{maxT} (the maximum score of every null
#'   profile), \code{calib}, \code{alpha}, \code{gumbel} (fitted location and
#'   scale), \code{null_exceed_scaled}, \code{fwer_holdout} (\code{NA} when no
#'   profile was held out), \code{widths}, \code{gamma}, \code{edge}.
#'
#' @seealso \code{\link{ras_box_detect}}, which consumes the result.
#'
#' @examples
#' ## 60 short synthetic null profiles (a real calibration uses several
#' ## hundred profiles from the actual scan on permuted phenotypes)
#' set.seed(2)
#' nulls <- lapply(1:60, function(k) 0.7 + abs(rnorm(200, sd = 0.3)))
#' cal <- ras_box_calibrate(nulls, alpha = 0.05, calib = 1:40)
#' cal$threshold
#' cal$fwer_holdout           # error rate on the 20 held-out profiles
#'
#' y <- nulls[[41]]; y[91:110] <- y[91:110] + 6
#' ras_box_detect(seq_along(y), y, calibration = cal)$regions
#' @export
ras_box_calibrate <- function(null_profiles, alpha = 0.05, widths = .ras_box_default_widths(),
                              gamma = 0.5, length_ratio = 1, calib = NULL,
                              edge = c("open", "strict"), method = c("order", "gumbel")) {
  edge <- match.arg(edge); method <- match.arg(method)
  NY <- if (is.list(null_profiles)) null_profiles
        else lapply(seq_len(nrow(null_profiles)), function(i) null_profiles[i, ])
  if (length(NY) < 2L) stop("ras_box_calibrate: at least two null profiles are needed", call. = FALSE)
  if (is.null(calib)) calib <- seq_along(NY)
  level0 <- stats::median(vapply(NY[calib], stats::median, numeric(1)))
  m <- vapply(NY, function(y) {
    S <- ras_box_stat(y, widths = widths, level0 = level0, gamma = gamma, edge = edge)
    if (is.null(S) || !nrow(S)) 0 else max(S$T)
  }, numeric(1))
  mc <- sort(m[calib]); B <- length(mc)
  thr_order <- mc[min(B, ceiling((B + 1) * (1 - alpha)))]
  b <- stats::sd(mc) * sqrt(6) / pi; mu <- mean(mc) - 0.5772156649 * b
  thr_gumbel <- mu - b * log(-log(1 - alpha))
  thr <- if (method == "gumbel") thr_gumbel else thr_order
  thr_scaled <- if (length_ratio > 1) mu - b * log(-log((1 - alpha)^(1 / length_ratio))) else thr
  structure(list(
    level0 = level0, threshold = thr, threshold_order = thr_order, threshold_gumbel = thr_gumbel,
    threshold_scaled = thr_scaled, method = method,
    length_ratio = length_ratio, maxT = m, calib = calib, alpha = alpha,
    gumbel = c(mu = mu, beta = b),
    null_exceed_scaled = sum(mc >= thr_scaled),
    fwer_holdout = if (length(calib) < length(NY)) mean(m[-calib] >= thr) else NA_real_,
    widths = widths, gamma = gamma, edge = edge), class = "ras_box_calibration")
}

#' @export
print.ras_box_calibration <- function(x, ...) {
  cat(sprintf("Box-scan calibration: %d null profile(s), alpha = %g, method = %s\n",
              length(x$maxT), x$alpha, if (is.null(x$method)) "order" else x$method))
  cat(sprintf("  level0 = %.3f | threshold = %.3f (order %.3f, gumbel %.3f)", x$level0, x$threshold,
              x$threshold_order, x$threshold_gumbel))
  if (x$length_ratio > 1)
    cat(sprintf(" | threshold_scaled (x%.1f length) = %.3f", x$length_ratio, x$threshold_scaled))
  cat("\n")
  if (is.finite(x$fwer_holdout))
    cat(sprintf("  held-out error rate = %.3f (%d profile(s))\n",
                x$fwer_holdout, length(x$maxT) - length(x$calib)))
  invisible(x)
}
