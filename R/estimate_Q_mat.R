estimate_Q_j <- function(data, folds, id, x, g, a_jm1, a_j, y_jm1, y_j, sbar_jm1, mu_jp1, Q_jp1, gval, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_jm1) == 1)
  }
  if (!is.null(y_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_jm1) == 1)
  }
  if (!is.null(a_j)) {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval & !!sym(a_j) == 1 & !!sym(y_j) == 1)
  } else {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval & !!sym(y_j) == 1)
  }

  at_risk_data <- mutate(at_risk_data, muQ = !!sym(mu_jp1) * !!sym(Q_jp1))
  at_risk_data <- mutate(at_risk_data, muQ = ifelse(!include_in_training, 0, muQ))
  Q_j <- estimate_cont(at_risk_data, folds, id, c(x, sbar_jm1), "muQ", "include_in_training", lrnr, paste0("Q", gval, "_", j))
  Q_j <- select(Q_j, !!id, !!paste0("Q", gval, "_", j))
  out_Q <- left_join(data, Q_j, by = id)
  out_Q
}

estimate_Q_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, all_mu, gval, lrnr, slim = FALSE) {
  tt <- length(all_y)
  tt_s <- length(all_s)
  Q_dat <- mutate(data, !!paste0("Q", gval, "_", tt) := 1)

  for (t in (tt - 1):1) {
    if (t == 1) {
      Q_dat <- estimate_Q_j(Q_dat, folds, id, x, g, NULL, all_a[1], NULL, all_y[1], NULL, all_mu[t + 1], paste0("Q", gval, "_", t + 1), gval, t, lrnr)
    } else {
      Q_dat <- estimate_Q_j(Q_dat, folds, id, x, g, all_a[t - 1], all_a[t], all_y[t - 1], all_y[t], all_s[1:min(tt_s, t - 1)], all_mu[t + 1], paste0("Q", gval, "_", t + 1), gval, t, lrnr)
    }
  }
  if (slim) {
    return(select(Q_dat, !!id, paste0("Q", gval, "_", 1:tt)))
  } else {
    Q_dat
  }
}
