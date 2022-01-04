estimate_Qstar_j <- function(data, folds, id, x, a_jm1, a_j, y_jm1, y_j, sbar_jm1, mu_jp1, Qstar_jp1, gval, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_jm1) == 1)
  }
  if (!is.null(y_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_jm1) == 1)
  }
  if (!is.null(a_j)) {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(a_j) == 1 & !!sym(y_j) == 1)
  } else {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(y_j) == 1)
    ## do we need this y_j condition?
  }

  at_risk_data <- mutate(at_risk_data, muQstar = !!sym(mu_jp1)*!!sym(Qstar_jp1))
  at_risk_data <- mutate(at_risk_data, muQstar = ifelse(!include_in_training, 0, muQstar))
  Qstar_j <- estimate_cont(at_risk_data, folds, c(x, sbar_jm1), 'muQstar', 'include_in_training', lrnr, paste0('Qstar', gval, '_', j))
  Qstar_j <- select(Qstar_j, -row_id, -include_in_training, -muQstar)
  out_Qstar <- left_join(data, Qstar_j)
  out_Qstar
}

estimate_Qstar_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, all_mu, gval, all_Q, t0 = length(all_s), lrnr, slim = FALSE) {
  tt <- length(all_y)

  tt_s <- length(all_s)
  Qstar_dat <- mutate(data, !!paste0('Qstar', gval, '_', tt) := 1)
  if (t0 < tt - 2) {
    for (ttt in (t0+2):tt) {
      Qstar_dat <- mutate(Qstar_dat, !!paste0('Qstar', gval, '_', ttt) := !!sym(all_Q[ttt]))
    }
  }

  tstart <- min(c(t0+1, tt-1))

  for (t in tstart:1) {
    if (t == 1) {
      Qstar_dat <- estimate_Qstar_j(Qstar_dat, folds, id, x, NULL, all_a[1], NULL, all_y[1], NULL, all_mu[t+1], paste0('Qstar', gval, '_', t+1), gval, t, lrnr)
    } else {
      Qstar_dat <- estimate_Qstar_j(Qstar_dat, folds, id, x, all_a[t-1], all_a[t], all_y[t-1], all_y[t], all_s[1:min(tt_s, t-1)], all_mu[t+1], paste0('Qstar', gval, '_', t+1), gval, t, lrnr)
    }
  }
  if (slim) {
    return(select(Qstar_dat, !!id, paste0('Qstar', gval, '_', 1:tt)))
  } else return(Qstar_dat)

}
