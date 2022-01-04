estimate_mu_j <- function(data, folds, id, x, g, a_jm1, a_j, y_j, y_jm1, sbar_jm1, gval, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_jm1) == 1)
  }
  if (!is.null(y_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_jm1) == 1)
  }
  if (!is.null(a_j)) {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(a_j) == 1 & !!sym(g) == gval)
  } else {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval)
  }
  estimate_binary(at_risk_data, folds, id, c(x, sbar_jm1), y_j, 'include_in_training', lrnr, paste0('mu', gval, '_', j))
}

estimate_mu_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, gval, lrnr, slim = FALSE) {
  tt <- length(all_y)
  t0 <- length(all_s)
  mu_js <- map(1:tt, function(t) {
    if (t == 1) {
      estimate_mu_j(data, folds, id, x, g, NULL, all_a[t], all_y[t], NULL, NULL, gval, t, lrnr)
    } else {
      estimate_mu_j(data, folds, id, x, g, all_a[t-1], all_a[t], all_y[t], all_y[t-1],all_s[1:min(t0, (t-1))], gval, t, lrnr)
    }

  })
  # mu_js_slim <- map(mu_js, ~select(., !!id, contains(paste0('mu', gval))))
  all_ids <- select(data, !!id)
  out_mu <- left_join(all_ids, mu_js[[1]], by = id)
  for (j in 2:tt) {
    out_mu <- left_join(out_mu, mu_js[[j]], by = id)
  }
  if (slim) {
    return(out_mu)
  } else {
    return(left_join(data, out_mu, by = id))
  }
}
