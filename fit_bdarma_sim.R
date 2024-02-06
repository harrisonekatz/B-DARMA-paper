library(cmdstanr)
library(dplyr)
library(ggplot2)
library(here)
library(loo)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

proj_dir <- here::here()

df <- readr::read_csv(file.path(proj_dir, "darma-1-1.csv"))

T_new <- 12

Y <-
  df |>
  dplyr::slice_head(n = -T_new) |>
  purrr::transpose() |>
  purrr::map(as.numeric)

C <- ncol(df)

# show the pattern for constructing X in a non-trivial case
X_component <- purrr::map(
  seq_along(Y),
  ~ purrr::map(seq(C - 1), ~ 1)
)

X <- purrr::map(X_component, purrr::flatten_dbl)

X_phi <- matrix(1, nrow = length(Y), ncol = 1)

X_new_component <- purrr::map(
  seq(T_new),
  ~ purrr::map(seq(C - 1), ~ 1)
)

X_new <- purrr::map(X_new_component, purrr::flatten_dbl)

stan_data <- list(
  `T` = length(Y),
  C = C,
  Y = Y,
  ref = 1,
  N = length(X[[1]]),
  K = purrr::map_int(X_component[[1]], length),
  X = X,
  K_phi = ncol(X_phi),
  X_phi = X_phi,
  P = 1,
  Q = 1,
  T_new = T_new,
  X_new = X_new,
  X_phi_new = X_phi[seq(T_new), , drop = FALSE],
  prior_only = 0,
  beta_sd = 0.3
)

cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.33.1/")
stan_dir <- file.path(proj_dir, "stan")
stan_file <- file.path(stan_dir, "bdarma.stan")

mod <- cmdstanr::cmdstan_model(stan_file, include_paths = stan_dir)

set.seed(1234)

fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  init = 0.975,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
)

# true parameters from simulate_darma.R. to recover these values, not only
# should the model reflect the order of the simulation model, e.g., DARMA (1,
# 1), but the index of the reference element for the ALR transformation must
# also be known and used in Stan.
true_param <- tibble::tribble(
  ~variable,     ~value,
  "beta[1]",     -0.07,
  "beta[2]",      0.1,
  "A[1,1,1]",     0.95,
  "A[1,2,1]",     0.3,
  "A[1,1,2]",    -0.18,
  "A[1,2,2]",     0.95,
  "B[1,1,1]",     0.65,
  "B[1,2,1]",     0.2,
  "B[1,1,2]",     0.15,
  "B[1,2,2]",     0.65,
  "beta_phi[1]",  log(1000)
)

fit$summary(c("beta", "beta_phi", "A", "B")) |>
  dplyr::inner_join(true_param, by = "variable") |>
  dplyr::rename(truth = value) |>
  dplyr::relocate(truth, .before = mean)

# check that we can compute PSIS-LOO-CV
loo_res <- fit$loo(cores = 4)
loo_res

Y_true <-
  df |>
  dplyr::slice_tail(n = T_new) |>
  dplyr::mutate(t = row_number()) |>
  tidyr::pivot_longer(starts_with("V")) |>
  dplyr::mutate(
    component = gsub("V", "", name),
    across(component, as.integer)
  ) |>
  dplyr::select(t, component, truth = value)

fit$draws("Y_hat", format = "draws_df") |>
  dplyr::select(starts_with("Y_hat")) |>
  tidyr::pivot_longer(everything()) |>
  dplyr::summarize(estimate = mean(value), .by = name) |>
  dplyr::mutate(idx = stringr::str_extract(name, "(?<=\\[).+?(?=\\])")) |>
  tidyr::separate(idx, c("t", "component"), sep = ",") |>
  dplyr::select(-name) |>
  dplyr::mutate(across(c(t, component), as.integer)) |>
  dplyr::inner_join(Y_true, by = c("t", "component")) |>
  dplyr::mutate(across(component, factor)) |>
  tidyr::pivot_longer(c(estimate, truth)) |>
  ggplot2::ggplot(aes(x = t, y = value)) +
  ggplot2::geom_line(aes(color = name)) +
  ggplot2::facet_wrap(vars(component)) +
  ggplot2::labs(
    title = "Simulated compositional time series",
    x = "time",
    y = element_blank()
  ) +
  ggplot2::theme(legend.title = element_blank())
