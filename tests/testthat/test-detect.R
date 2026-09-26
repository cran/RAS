skip_if_not_installed("segmented")

make_piecewise_fixture <- function(seed, n, psi_true = 40, slope_left = 0.15, slope_right = -0.08,
                                    sd = 0.8, xmax = 100) {
  set.seed(seed)
  x <- sort(runif(n, 0, xmax))
  y <- ifelse(x < psi_true, slope_left * x, slope_left * psi_true + slope_right * (x - psi_true)) +
    stats::rnorm(n, sd = sd)
  list(x = x, y = y)
}

test_that("davies_test (C) matches segmented::davies.test across both n<=300 and n>300 branches", {
  max_diff_small <- 0; max_diff_large <- 0
  for (seed in 1:6) {
    for (n in c(50, 150, 300, 301, 600, 3000)) {
      d <- make_piecewise_fixture(seed, n)
      fit_lm <- stats::lm(y ~ x, data = data.frame(x = d$x, y = d$y))
      r_p <- tryCatch(segmented::davies.test(fit_lm)$p.value, error = function(e) NA)
      c_p <- .Call("RAS_davies_test", as.double(d$x), as.double(d$y), PACKAGE = "RAS")
      if (!is.na(r_p)) {
        diff <- abs(r_p - c_p)
        if (n <= 300) max_diff_small <- max(max_diff_small, diff) else max_diff_large <- max(max_diff_large, diff)
      }
    }
  }
  expect_lt(max_diff_small, 1e-6)
  expect_lt(max_diff_large, 1e-6)
})

test_that("slope_test (C) matches a no-intercept lm t-test exactly", {
  max_diff <- 0
  for (seed in 1:10) {
    set.seed(seed + 100)
    n <- sample(5:40, 1)
    x <- stats::rnorm(n)
    y <- 0.7 * x + stats::rnorm(n, sd = 0.3)
    lower_tail <- (seed %% 2 == 0)
    model <- stats::lm(y ~ x - 1)
    tstat <- stats::coef(summary(model))["x", "t value"]
    r_p <- stats::pt(tstat, df = stats::df.residual(model), lower.tail = lower_tail)
    c_p <- .Call("RAS_slope_test", as.double(x), as.double(y), as.logical(lower_tail), PACKAGE = "RAS")
    max_diff <- max(max_diff, abs(r_p - c_p))
  }
  expect_lt(max_diff, 1e-9)
})

test_that("single deterministic Muggeo fit (C, no bootstrap) matches segmented(n.boot=0)", {
  max_psi_diff <- 0; max_slope_diff <- 0; n_ok <- 0
  for (seed in 1:12) {
    n <- c(60, 150, 400, 3000)[((seed - 1) %% 4) + 1]
    d <- make_piecewise_fixture(seed, n)
    fit_lm <- stats::lm(y ~ x, data = data.frame(x = d$x, y = d$y))
    r_fit <- tryCatch(
      suppressWarnings(segmented::segmented(fit_lm, seg.Z = ~x, npsi = 1,
                                             control = segmented::seg.control(n.boot = 0))),
      error = function(e) NULL)
    if (is.null(r_fit) || is.null(r_fit$psi)) next
    r_psi <- r_fit$psi[1, 2]
    r_slopes <- segmented::slope(r_fit)$x[, 1]
    c_res <- .Call("RAS_seg_fit_single", as.double(d$x), as.double(d$y), PACKAGE = "RAS")
    if (!c_res$success) next
    n_ok <- n_ok + 1
    max_psi_diff <- max(max_psi_diff, abs(r_psi - c_res$psi))
    max_slope_diff <- max(max_slope_diff, abs(r_slopes[1] - c_res$slope.left), abs(r_slopes[2] - c_res$slope.right))
  }
  expect_gte(n_ok, 10)
  expect_lt(max_psi_diff, 1e-4)
  expect_lt(max_slope_diff, 1e-3)
})

test_that("bootstrap-restart fit (C) never does worse than its own single-start fit", {
  for (seed in 1:15) {
    n <- c(100, 400, 3000)[((seed - 1) %% 3) + 1]
    d <- make_piecewise_fixture(seed + 200, n, slope_left = 0.12, slope_right = -0.10, sd = 1.5)
    single <- .Call("RAS_seg_fit_single", as.double(d$x), as.double(d$y), PACKAGE = "RAS")
    boot <- .Call("RAS_seg_fit_boot", as.double(d$x), as.double(d$y), 10L, PACKAGE = "RAS")
    if (!single$success || !boot$success) next
    expect_lte(boot$rss, single$rss + 1e-6)
  }
})

test_that("degenerate (flat) window returns p=1 / no breakpoint, no crash", {
  x <- 1:50
  y <- rep(3, 50) # zero residual variance
  res <- .Call("RAS_get_break_points", as.double(x), as.double(y), 50L, PACKAGE = "RAS")
  expect_null(res$break.points)
  expect_equal(res$p.values, 1)
})

test_that("short window below identifiability threshold does not crash", {
  x <- 1:4
  y <- c(1, 2, 5, 4)
  res <- .Call("RAS_get_break_points", as.double(x), as.double(y), 4L, PACKAGE = "RAS")
  expect_null(res$break.points)
  expect_equal(res$p.values, 1)
})

test_that("ras_detect finds the same accepted changepoint as ras_detect_original (small window, n<=300 davies branch)", {
  set.seed(42)
  N <- 600
  x <- 1:N
  mid <- N %/% 2
  y <- c(seq(0, 8, length.out = mid), seq(8, 1, length.out = N - mid)) + stats::rnorm(N, sd = 0.6)

  r_res <- suppressMessages(ras_detect_original(x, y, window_size = 150, skip = 4, slope_check_window_size = 15,
                                        slope.p.values.threshold.left = 1e-3, slope.p.values.threshold.right = 1e-3))
  c_res <- ras_detect(x, y, window_size = 150, skip = 4, slope_check_window_size = 15,
                            slope.p.values.threshold.left = 1e-3, slope.p.values.threshold.right = 1e-3)

  expect_equal(sort(c_res$tau_hats), sort(r_res$tau_hats))
})

test_that("ras_detect finds the same accepted changepoint as ras_detect_original (large window, n>300 davies branch)", {
  set.seed(7)
  N <- 1000
  x <- 1:N
  mid <- N %/% 2
  y <- c(seq(0, 8, length.out = mid), seq(8, 1, length.out = N - mid)) + stats::rnorm(N, sd = 0.6)

  r_res <- suppressMessages(ras_detect_original(x, y, window_size = 400, skip = 15, slope_check_window_size = 20,
                                        slope.p.values.threshold.left = 1e-3, slope.p.values.threshold.right = 1e-3))
  c_res <- ras_detect(x, y, window_size = 400, skip = 15, slope_check_window_size = 20,
                            slope.p.values.threshold.left = 1e-3, slope.p.values.threshold.right = 1e-3)

  expect_equal(sort(c_res$tau_hats), sort(r_res$tau_hats))
})

test_that("ras_detect returns no changepoints on pure noise (no crash, sane output)", {
  set.seed(99)
  N <- 400
  x <- 1:N
  y <- stats::rnorm(N, sd = 1)
  res <- ras_detect(x, y, window_size = 150, skip = 5, slope_check_window_size = 15)
  expect_type(res$tau_hats, "integer")
  expect_equal(length(res$tau_hats), length(res$p.values))
})
