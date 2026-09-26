#' Aliases Kept from the Development Versions
#'
#' In the development versions of this package the disk-backed C
#' implementations carried a \code{_fast} suffix. Since RAS 1.1.0 they are
#' the default implementations and carry the plain names, while the original
#' pure-R, in-memory implementations carry an \code{_original} suffix. The
#' \code{_fast} names are kept as thin aliases so that scripts written against
#' the development versions keep working; new code should use the plain names.
#'
#' @param ... Arguments passed on unchanged to the aliased function.
#'
#' @return Whatever the aliased function returns.
#'
#' @seealso \code{\link{ras}}, \code{\link{ras_scan}}, \code{\link{ras_detect}},
#'   \code{\link{compute_gwas_weights}}, \code{\link{screen_forward_max_region}},
#'   \code{\link{ras_scan_external}}.
#'
#' @examples
#' identical(formals(ras_fast), formals(ras))   # FALSE only because of `...`
#' @name ras-aliases
NULL

#' @rdname ras-aliases
#' @export
ras_fast <- function(...) ras(...)

#' @rdname ras-aliases
#' @export
ras_scan_fast <- function(...) ras_scan(...)

#' @rdname ras-aliases
#' @export
ras_detect_fast <- function(...) ras_detect(...)

#' @rdname ras-aliases
#' @export
compute_gwas_weights_fast <- function(...) compute_gwas_weights(...)

#' @rdname ras-aliases
#' @export
screen_forward_max_region_fast <- function(...) screen_forward_max_region(...)

#' @rdname ras-aliases
#' @export
ras_scan_external_fast <- function(...) ras_scan_external(...)
