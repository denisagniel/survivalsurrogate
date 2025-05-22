eif_delta <- function(y, a, g, e, gamma0, mu0, Q0, gamma1, mu1, Q1) {
  tt <- ncol(y)
  part_j1 <- map(1:tt, ~ eif_delta_partj(., y, a, gamma1, mu1, Q1, tt)) |>
    bind_cols()
  part_j0 <- map(1:tt, ~ eif_delta_partj(., y, a, gamma0, mu0, Q0, tt)) |>
    bind_cols()
  g / e * rowSums(part_j1) - (1 - g) / (1 - e) * rowSums(part_j0) + mu1[, 1] * Q1[, 1] - mu0[, 1] * Q0[, 1]
}

eif_delta_partj <- function(j, y, a, gamma, mu, Q, tt) {
  ay_gamma <- a * y / gamma
  if (j == 1) {
    ay_gamma_jm1 <- 1
  } else {
    ay_gamma_jm1 <- apply(ay_gamma[, 1:(j - 1), drop = FALSE], 1, prod)
  }

  if (j == tt) {
    tibble(!!glue("part{j}") := ay_gamma_jm1 * a[, j] / gamma[, j] * (y[, j] - mu[, j]))
  } else {
    tibble(!!glue("part{j}") := ay_gamma_jm1 * a[, j] / gamma[, j] * (y[, j] * mu[, j + 1] * Q[, j + 1] - mu[, j] * Q[, j]))
  }
}
