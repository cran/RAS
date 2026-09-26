# Box-scan region detector (R/box_scan.R). Synthetic profiles only; the
# regression checks against the real chr16 / AoU results live outside the
# package (paired_boundary_detector/tests/test_box_scan.R) because they need
# data that is not shipped.

# a RAS-like null: positive, right-skewed, autocorrelated, median near 0.7
null_profile <- function(n, seed) {
  set.seed(seed)
  0.7 + as.numeric(stats::filter(stats::rexp(n + 20, 3) - 1 / 3,
                                 rep(1 / 5, 5), sides = 2))[11:(n + 10)]
}

test_that("noiseless plateau: exact edges, score = height, interior windows score 0", {
  y <- rep(1, 600); y[301:342] <- 9
  S <- ras_box_stat(y)
  r <- ras_box_detect(seq_along(y), y, threshold = 1)$regions
  expect_equal(nrow(r), 1L)
  expect_equal(r$tau_L, 301L)
  expect_equal(r$tau_R, 342L)
  expect_equal(r$T_box, 8)
  # a 20-point window inside the plateau: one flank is as high as the window
  expect_lte(max(S$T[S$w == 20 & S$i >= 305 & S$i + 19 <= 338]), 1e-12)
})

test_that("three neighbouring peaks stay three regions", {
  y <- rep(1, 800); for (a in c(300, 360, 420)) y[a:(a + 19)] <- 7
  r <- ras_box_detect(seq_along(y), y, threshold = 1)$regions
  expect_equal(nrow(r), 3L)
  expect_equal(sort(r$tau_L), c(300L, 360L, 420L))
  expect_true(all(r$width == 20L))
})

test_that("the shoulders of a tall Gaussian peak are not separate regions", {
  y <- rep(1, 600)
  y[281:320] <- 1 + 30 * exp(-((281:320) - 300.5)^2 / (2 * 6^2))
  r <- ras_box_detect(seq_along(y), y, threshold = 1)$regions
  expect_equal(nrow(r), 1L)
  expect_lte(r$tau_L, 300L)
  expect_gte(r$tau_R, 301L)
})

test_that("raised background: level factor applied, small bump rejected, tall region kept", {
  set.seed(4); n <- 572
  y <- 6 + 2 * sin(seq(0, 6 * pi, length.out = n)) + rnorm(n, sd = 0.15)
  y[100:120] <- y[100:120] + 2.0
  y[400:424] <- y[400:424] + 14
  d <- ras_box_detect(seq_len(n), y, threshold = 2.5, level0 = 0.7)
  expect_equal(d$level_factor, sqrt(stats::median(y) / 0.7))
  expect_gt(d$level_factor, 2.5)
  expect_equal(nrow(d$regions), 1L)
  expect_gte(d$regions$tau_L, 396L)
  expect_lte(d$regions$tau_R, 428L)
  # level0 = NA disables the adjustment without error
  d0 <- ras_box_detect(seq_len(n), y, threshold = 2.5, level0 = NA)
  expect_equal(d0$level_factor, 1)
  expect_gte(nrow(d0$regions), 1L)
})

test_that("edge = 'open' can report a region at the profile end; 'strict' cannot", {
  y <- rep(1, 400); y[1:20] <- 9
  r_open   <- ras_box_detect(seq_along(y), y, threshold = 1, edge = "open")$regions
  r_strict <- ras_box_detect(seq_along(y), y, threshold = 1, edge = "strict")$regions
  expect_equal(nrow(r_open), 1L)
  expect_equal(r_open$tau_L, 1L)
  expect_equal(r_open$tau_R, 20L)
  expect_null(r_strict)
})

test_that("calibration: order-statistic threshold, held-out error rate, length scaling", {
  NY <- lapply(1:200, function(k) null_profile(400, k))
  cal <- ras_box_calibrate(NY, alpha = 0.05, calib = 1:120, length_ratio = 20)
  expect_s3_class(cal, "ras_box_calibration")
  expect_lte(sum(cal$maxT[1:120] >= cal$threshold), floor(0.05 * 121))
  expect_lte(cal$fwer_holdout, 0.15)
  expect_gt(cal$threshold_scaled, cal$threshold)
  expect_output(print(cal), "threshold")
  # both read-offs are returned; method selects which one is `threshold`
  expect_equal(cal$threshold, cal$threshold_order)
  calg <- ras_box_calibrate(NY, alpha = 0.05, calib = 1:120, method = "gumbel")
  expect_equal(calg$threshold, calg$threshold_gumbel)
  expect_equal(calg$threshold_order, cal$threshold_order)
  expect_lt(abs(calg$threshold_gumbel / cal$threshold_order - 1), 0.25)
  # a fresh null profile reports (almost) nothing at the scaled threshold
  d <- ras_box_detect(1:400, NY[[121]], calibration = cal, scaled = TRUE)
  expect_true(is.null(d$regions) || nrow(d$regions) <= 1)
  # the calibration's settings are picked up, explicit arguments override
  expect_equal(d$threshold, cal$threshold_scaled)
  expect_equal(d$level0, cal$level0)
  d2 <- ras_box_detect(1:400, NY[[121]], calibration = cal, threshold = 99)
  expect_equal(d2$threshold, 99)
  expect_null(d2$regions)
})

test_that("positions, anchors, tau_hats alias and the plot-compatible fields", {
  y <- rep(1, 600)
  y[301:342] <- c(seq(5, 9, length.out = 21), seq(9, 5, length.out = 21))
  x <- 1000L + 10L * seq_along(y)
  d <- ras_box_detect(x, y, threshold = 1)
  r <- d$regions
  expect_equal(r$anchor_pos, x[r$anchor_idx])
  expect_equal(r$anchor_y, max(y))
  expect_equal(r$pos_L, x[r$tau_L])
  expect_equal(r$pos_R, x[r$tau_R])
  expect_identical(d$tau_hats, r$anchor_pos)
  expect_identical(d$all.p.values, r$T_box)
  expect_null(d$left.slopes)
  expect_identical(d$detector, "box")
  expect_null(d$stat)
  expect_true(is.data.frame(ras_box_detect(x, y, threshold = 1, keep_stat = TRUE)$stat))
})

test_that("input checks", {
  expect_error(ras_box_stat(c(1, NA, 3, rep(1, 50))), "NA")
  expect_error(ras_box_stat(rep(1, 12)), "too short")
  expect_error(ras_box_detect(1:10, rep(1, 600), threshold = 1), "differ in length")
  expect_error(ras_box_detect(1:600, rep(1, 600)), "threshold")
  expect_error(ras_box_detect(1:600, rep(1, 600), calibration = list(threshold = 1)),
               "ras_box_calibrate")
  expect_error(ras_box_calibrate(list(rep(1, 100))), "at least two")
  expect_true(all(ras_box_stat(rep(1, 60))$w <= 28))
})

test_that("plot and print methods accept a box-scan detection", {
  set.seed(9)
  x <- seq(1, 6000, by = 10)
  y <- 0.7 + abs(rnorm(600, sd = 0.3)); y[301:342] <- y[301:342] + 8
  det <- ras_box_detect(x, y, threshold = 2)
  res <- structure(list(scan = list(x = x, y = y), detection = det,
                        chrom = 1, save_dir = tempdir()), class = "ras")
  expect_output(print(res), "box-scan detector")
  expect_output(print(res), "region")
  out_dir <- file.path(tempdir(), "ras_box_plot_test"); dir.create(out_dir, showWarnings = FALSE)
  res$save_dir <- out_dir
  # pdf() exists in every R build; png() needs X11/cairo and is absent from
  # minimal builds (e.g. rhub/r-minimal on Alpine), so it is only tried when available
  expect_no_error(plot(res, device = "pdf"))
  expect_no_error(plot(res, zoom = TRUE, device = "pdf"))
  expect_true(file.exists(file.path(out_dir, "chr-1-cp-plot.pdf")))
  expect_true(file.exists(file.path(out_dir, "chr-1-cp-p-values-plot.pdf")))
  expect_true(file.exists(file.path(out_dir, "chr-1-zoom.pdf")))
  if (isTRUE(capabilities("png"))) {
    expect_no_error(plot(res, device = "png"))
    expect_true(file.exists(file.path(out_dir, "chr-1-cp-plot.png")))
  }
  unlink(out_dir, recursive = TRUE)
})
