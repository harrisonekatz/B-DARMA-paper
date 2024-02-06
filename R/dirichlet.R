#' The Dirichlet Distribution
#'
#' Density function and random number generation for the Dirichlet distribution
#' with shape parameter vector `alpha`.
#'
#' @name Dirichlet
#'
#' @param x Matrix of quantiles. Each row corresponds to one probability vector.
#' @param n Number of draws to sample from the distribution.
#' @param alpha Matrix of positive shape parameters. Each row corresponds to one
#'   probability vector.
#' @param log Logical; If `TRUE`, values are returned on the log scale.
#'
#' @details See \code{vignette("brms_families")} for details on the
#'   parameterization.
#'
#' @export
ddirichlet <- function(x, alpha, log = FALSE) {
  checkmate::assert_flag(log)

  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
  }
  if (nrow(x) == 1L && nrow(alpha) > 1L) {
    x <- replicate(nrow(alpha), x, simplify = FALSE)
    x <- do.call(rbind, x)
  } else if (nrow(x) > 1L && nrow(alpha) == 1L) {
    alpha <- replicate(nrow(x), alpha, simplify = FALSE)
    alpha <- do.call(rbind, alpha)
  }
  if (isTRUE(any(x < 0))) {
    stop("x must be non-negative.")
  }
  if (!isTRUE(all.equal(rowSums(x), rep(1, nrow(x))))) {
    stop("x must sum to 1 per row.")
  }
  if (isTRUE(any(alpha <= 0))) {
    stop("alpha must be positive.")
  }
  out <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) +
    rowSums((alpha - 1) * log(x))
  if (!log) {
    out <- exp(out)
  }
  return(out)
}

#' @rdname Dirichlet
#' @export
rdirichlet <- function(n, alpha) {
  checkmate::assert_count(n)

  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop("alpha should be non-empty.")
  }
  if (isTRUE(any(alpha <= 0))) {
    stop("alpha must be positive.")
  }
  if (n == 1) {
    n <- nrow(alpha)
  }
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, nrow = n, ncol = ncol(alpha), byrow = TRUE)
  }
  x <- matrix(rgamma(ncol(alpha) * n, alpha), ncol = ncol(alpha))
  x / rowSums(x)
}
