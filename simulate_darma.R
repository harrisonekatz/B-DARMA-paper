library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(tidyr)

proj_dir <- here::here()
source(file.path(proj_dir, "R", "compositional.R"))
source(file.path(proj_dir, "R", "dirichlet.R"))

C <- 3
P <- 1
Q <- 1
M <- max(P, Q)
T_length <- 500

# VAR coefficients
A <- matrix(c(0.95, 0.3, -0.18, 0.95), nrow = C - 1)

# VMA coefficients
B <- matrix(c(0.65, 0.2, 0.15, 0.65), nrow = C - 1)

# intercepts
beta <- matrix(c(-0.07, 0.1), nrow = C - 1)

set.seed(1234)
y_obs <- matrix(runif(C * M, 0, 0.4), nrow = M)
y_obs <- y_obs / rowSums(y_obs)

y <- matrix(data = rep(0, C * T_length), ncol = C)
y[1:M, ] <- y_obs

alr_y <- matrix(0, nrow = T_length, ncol = C - 1)

alr_y[1:M, ] <-
  y[1:M, , drop = FALSE] |>
  apply(1, alr) |>
  unlist() |>
  matrix(ncol = C - 1, byrow = TRUE)

eta <- matrix(0, nrow = T_length, ncol = C - 1)

# set the first M elements of eta to alr(y)
eta[1:M, ] <- alr_y[1:M, ]

phi <- 1000

for (t in seq(M + 1, T_length)) {
  # this process has only an intercept, so X is just 1s
  eta[t, ] <- A %*% (alr_y[t - 1, ] - beta) +
    B %*% (alr_y[t - 1, ] - eta[t - 1, ]) +
    beta

  alpha <- alrinv(eta[t, ]) * phi
  y_proposed <- rdirichlet(1, alpha)

  # avoid numerical issues
  while(min(y_proposed) < 0.005) {
    y_proposed <- rdirichlet(1, alpha)
  }

  y[t, ] <- y_proposed
  alr_y[t, ] <- alr(y[t, ])
}

# PLOT: the compositional series in the transformed space
alr_y |>
  as.data.frame() |>
  dplyr::mutate(idx = row_number()) |>
  tidyr::pivot_longer(-idx) |>
  dplyr::mutate(across(name, factor)) |>
  ggplot2::ggplot(aes(x = idx, y = value)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(vars(name))

y |>
  as.data.frame() |>
  dplyr::mutate(idx = row_number()) |>
  tidyr::pivot_longer(-idx) |>
  dplyr::mutate(across(name, factor)) |>
  ggplot2::ggplot(aes(x = idx, y = value)) +
  ggplot2::geom_area(aes(fill = name))

readr::write_csv(as.data.frame(y), "darma-1-1.csv")
