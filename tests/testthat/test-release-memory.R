test_that("release_memory() runs without error", {
  expect_no_error(release_memory(verbose = FALSE))
})

test_that("release_memory() returns an integer (or NA_integer_)", {
  out <- release_memory(verbose = FALSE)
  expect_true(is.integer(out) || is.na(out))
})

# malloc_trim() is a glibc extension. The build records whether it was compiled
# in (Linux/glibc: yes; Windows, macOS, musl/Alpine Linux, ...: no), and the
# contract differs accordingly: 0/1 where it exists, NA_integer_ elsewhere.
test_that("release_memory() returns 0 or 1 where malloc_trim() was compiled in", {
  skip_if_not(RAS:::.ras_malloc_trim_available(), "malloc_trim() not available in this build")
  out <- release_memory(verbose = FALSE)
  expect_true(out %in% c(0L, 1L))
})

test_that("release_memory() returns NA where malloc_trim() was not compiled in", {
  skip_if(RAS:::.ras_malloc_trim_available(), "malloc_trim() available in this build")
  out <- release_memory(verbose = FALSE)
  expect_true(is.na(out))
})

test_that("malloc_trim() availability matches the platform", {
  avail <- RAS:::.ras_malloc_trim_available()
  is_linux <- identical(.Platform$OS.type, "unix") && grepl("linux", R.version$os)
  if (!is_linux) expect_false(avail)     # never on Windows / macOS / *BSD
  expect_type(avail, "logical")
})

test_that("release_memory() emits a message when verbose = TRUE", {
  expect_message(release_memory(verbose = TRUE), "malloc_trim returned")
})

test_that("release_memory() emits no message when verbose = FALSE", {
  expect_no_message(release_memory(verbose = FALSE))
})

test_that("release_memory() return value is invisible", {
  # withVisible() reports whether the return was auto-printed
  rv <- withVisible(release_memory(verbose = FALSE))
  expect_false(rv$visible)
})
