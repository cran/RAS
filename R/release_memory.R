#' Release Memory Back to the OS
#'
#' Runs full garbage collection and, on Linux with the GNU C library (glibc),
#' calls \code{malloc_trim(0)} to return free heap pages to the operating
#' system. Useful after large temporary matrices are removed in memory-heavy
#' RAS pipeline stages.
#'
#' \code{malloc_trim()} is a glibc extension. On every other platform,
#' including Linux systems with another C library such as musl (Alpine
#' Linux), Windows and macOS, the call is not compiled in and \code{NA} is
#' returned silently; no error is raised.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a one-line message
#'   with the \code{malloc_trim} return value so RSS changes can be monitored
#'   in pipeline logs.
#'
#' @return Invisibly returns the \code{malloc_trim(0)} result:
#'   \code{1} if heap pages were returned to the OS,
#'   \code{0} if nothing was returned,
#'   \code{NA_integer_} where \code{malloc_trim()} is unavailable (all
#'   platforms other than glibc Linux).
#' @export
release_memory <- function(verbose = TRUE) {
  invisible(gc(full = TRUE))

  out <- tryCatch(
    .Call("RAS_malloc_trim", PACKAGE = "RAS"),
    error = function(e) NA_integer_
  )

  if (isTRUE(verbose)) {
    message("gc(full = TRUE) complete; malloc_trim returned: ", out)
  }

  invisible(out)
}

# TRUE when this build of the package was compiled with glibc's malloc_trim()
# (Linux/glibc), FALSE elsewhere. Internal; used by the unit tests to decide
# what release_memory() must return on the current platform.
.ras_malloc_trim_available <- function() {
  isTRUE(.Call("RAS_malloc_trim_available", PACKAGE = "RAS"))
}
