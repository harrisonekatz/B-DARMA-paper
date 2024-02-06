#' Additive log ratio transformation
#'
#' @param x numeric vector, must be a simplex
#' @param denom which element of `x` to use as the denominator
#' @param `...` other arguments passed to `is_compositional()`
#'
#' @return the transform of `x`
alr <- function(x, denom = 1, ...) {
  checkmate::assert_true(is_compositional(x, ...))

  log(x[-denom] / x[denom])
}

#' Inverse additive log ratio transformation
#'
#' @param x numeric vector
#' @inheritParams alr
#'
#' @return the inverse of `alr(x)`, where the reference element is equal to 1
alrinv <- function(x, denom = 1) {
  xlen <- length(x)
  checkmate::assert_numeric(x)
  checkmate::assert_number(denom, lower = 1, upper = xlen + 1)

  if (denom == 1) {
    res <- c(0, x)
  } else if (denom == xlen + 1) {
    res <- c(x, 0)
  } else {
    res <- c(x[seq(1, denom - 1)], 0, x[seq(denom, xlen)])
  }

  softmax(res)
}

#' Center log ratio transformation
#'
#' The inverse of the CLR transform is the softmax function.
#'
#' @inheritParams alr
#'
#' @return the transform of `x`
clr <- function(x, ...) {
  checkmate::assert_true(is_compositional(x, ...))

  # we know x > 0 courtesy of `is_compositional()`
  log_x <- log(x)
  log_x - mean(log_x)
}

#' Compress compositional data
#'
#' Compositional data may include 1s or 0s, which can cause problems for certain
#' modeling approaches. Maier (2014) provides a scheme for compressing the data
#' to remove those 1s and 0s.
#'
#' @param x numeric vector in `[0, 1]`
#' @param C the number of components
#' @param N the number of observations, where one observation consists of `C`
#'   components
#'
#' @return a numeric vector compressed as described
compress <- function(x, C, N = length(x) / C) {
  checkmate::assert_numeric(x, lower = 0, upper = 1)
  checkmate::assert_count(C)
  checkmate::assert_count(N)

  (x * (N - 1) + 1 / C) / N
}

#' Create a sequence of cosine terms for a Fourier series
#'
#' @param x a vector of times
#' @param ub upper bound for series in months
#' @param p the seasonal period
#'
#' @return vector of cosine terms
cos_term <- function(x, ub = 10, p = 12) {
  checkmate::assert_numeric(x)
  checkmate::assert_count(ub)
  checkmate::assert_count(p)

  cos(seq(2, ub, by = 2) * pi * x / p)
}

#' Decompress previously compressed compositional data
#'
#' Observe that `compress()` is an affine transformation. Thus, when
#' backtransforming predictions, for example, we should use the same `N` as was
#' used for the original transformation.
#'
#' @inheritParams compress
#'
#' @return a decompressed numeric vector
decompress <- function(x, C, N) {
  checkmate::assert_numeric(x, lower = 0, upper = 1)
  checkmate::assert_count(C)
  checkmate::assert_count(N)

  (N * x - 1 / C) / (N - 1)
}

#' Make a rotation matrix suitable for calculating the ILR transform
#'
#' @param D positive integer giving the number of columns
#'
#' @return a `D - 1` by `D` matrix whose rows each have unit length
ilr_rotation_matrix <- function(D) {
  checkmate::assert_int(D, lower = 2)

  # https://stats.stackexchange.com/a/259223
  H <- stats::contr.helmert(D)
  t(H) / sqrt(2:D * (2:D - 1))
}

#' Isometric log ratio transformation
#'
#' @inheritParams alr
#'
#' @return the transform of `x`
ilr <- function(x, ...) {
  checkmate::assert_true(is_compositional(x, ...))

  as.vector(tcrossprod(clr(x, ...), ilr_rotation_matrix(length(x))))
}

#' Inverse isometric log ratio transformation
#'
#' @inheritParams alrinv
#'
#' @return the inverse of `ilr(x)`
ilrinv <- function(x) {
  checkmate::assert_numeric(x)

  softmax(x %*% ilr_rotation_matrix(length(x) + 1))
}

#' Is a vector compositional?
#'
#' @param x numeric vector
#' @param ... other arguments passed to [base::all.equal()]
#'
#' @return logical indicating whether `x` is compositional
is_compositional <- function(x, ...) {
  checkmate::assert_numeric(x)

  all(x > 0) && isTRUE(all.equal(sum(x), 1, ...))
}

#' Does a tibble contain boundary values?
#'
#' @param .data a tibble
#' @param lower number giving the lower boundary
#' @param upper number giving the upper boundary
#' @inheritParams tidyr::pivot_longer
#'
#' @return logical indicating whether any values are on the boundary
has_boundary_value <- function(.data, cols, lower = 0, upper = 1) {
  checkmate::assert_data_frame(.data)
  checkmate::assert_number(lower)
  checkmate::assert_number(upper)

  # borrowed from pivot_longer()
  cols <- tidyselect::eval_select(
    rlang::enquo(cols),
    .data[unique(names(.data))]
  )
  cols <- names(cols)

  .data |>
    tidyr::pivot_longer(all_of(cols)) |>
    dplyr::filter(near(value, 0) | near(value, 1)) |>
    nrow() |>
    as.logical()
}

#' Create a matrix of Fourier terms
#'
#' @param n number of Fourier terms
#' @param ... additional arguments passed to `sin_term()` and `cos_term()`
#' @inheritParams cos_term
#'
#' @return a matrix of Fourier terms with one row for each element of `x`, where
#'   the terms are ordered, i.e., sines and cosines are interleaved
make_fourier <- function(x, n = 11, ...) {
  checkmate::assert_numeric(x)
  checkmate::assert_count(n)

  # Fourier terms alternate between sine and cosine, so 11 = 6 sines, 5 cosines
  nsin <- ceiling(n / 2)

  if (n %% 2 == 0) {
    ncos <- n / 2
  } else {
    ncos <- floor(n / 2)
  }

  sines <- vapply(x, sin_term, numeric(nsin), ub = nsin * 2, ...)
  cosines <- vapply(x, cos_term, numeric(ncos), ub = ncos * 2, ...)

  stacked <- rbind(sines, cosines)
  t(stacked[order(sequence(c(nsin, ncos))), ])
}

#' Pivot compositional data wider
#'
#' @param d unquoted name of the variable that indexes the compositional data
#' @param y unquoted name of the variable that holds the compositional values
#' @inheritParams prepare_bdarma
#'
#' @return a wider data frame
pivot_compositional <- function(.data, x, d, y) {
  checkmate::assert_data_frame(.data)

  .data |>
    # just in case
    dplyr::ungroup() |>
    dplyr::select({{ x }}, {{ d }}, {{ y }}) |>
    tidyr::pivot_wider(
      names_from = {{ d }},
      values_from = {{ y }},
      values_fill = 0
    ) |>
    dplyr::arrange({{ x }})
}

#' Prepare compositional data model
#'
#' * compress data containing values within a tolerance of 0 or 1
#'
#' * remove observations after the forecast start date
#'
#' @param .data a tibble
#' @param forecast_start_date the `Date` when the forecast should start; the
#'   first of the month will be used regardless of the actual day
#' @param x unquoted name of the `Date` variable that indexes observations
#' @param has_boundary logical indicating whether the data includes values on
#'   the boundary
#'
#' @return a tibble of prepared data
prepare_compositional <- function(.data,
                                  forecast_start_date,
                                  x,
                                  has_boundary = FALSE) {
  checkmate::assert_data_frame(.data)
  checkmate::assert_date(forecast_start_date, len = 1L)
  checkmate::assert_flag(has_boundary)

  fom <- lubridate::floor_date(forecast_start_date, unit = "month")

  # remove data on or after `fom`
  df_training <-
    .data |>
    dplyr::filter({{ x }} < fom) |>
    tidyr::pivot_longer(- {{ x }})

  # compress values if there are any 1s or 0s
  C <- length(unique(df_training[["name"]]))

  if (has_boundary) {
    df_training <- dplyr::mutate(
      df_training,
      across(value, \(x) compress(x, C = C))
    )
  }

  df_training
}

#' Create a sequence of sine terms for a Fourier series
#'
#' @inheritParams cos_term
#'
#' @return vector of sine terms
sin_term <- function(x, ub = 12, p = 12) {
  checkmate::assert_numeric(x)
  checkmate::assert_count(ub)
  checkmate::assert_count(p)

  sin(seq(2, ub, by = 2) * pi * x / p)
}

#' Softmax
#'
#' @param x numeric vector
#'
#' @return the softmax of `x`
softmax <- function(x) {
  checkmate::assert_numeric(x)

  exp(x) / sum(exp(x))
}

#' Tidy the results of a statistical test
#'
#' @param x the result of a statistical test
#' @param ... further arguments passed to methods
#'
#' @return a tibble of tidied results
tidy_test <- function(x, ...) {
  UseMethod("tidy_test")
}

#' Tidy the results of [urca::ur.df()]
#'
#' @param x an object of class `ur.df`
tidy_test.ur.df <- function(x, ...) {
  lhs <-
    x@teststat |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(everything(), names_to = "teststat")

  rhs <- tibble::as_tibble(x@cval, rownames = "teststat")
  dplyr::inner_join(lhs, rhs, by = "teststat")
}
