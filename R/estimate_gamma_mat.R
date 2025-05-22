estimate_gamma_j <- function(data, folds, id, x, g, a_j, a_m1, y_m1, sbar_m1, gval, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_m1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_m1) == 1)
  }
  if (!is.null(y_m1)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_m1) == 1)
  }
  at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval)
  estimate_binary(at_risk_data, folds, id, c(x, sbar_m1), a_j, "include_in_training", lrnr, paste0("gamma", gval, "_", j))
}

estimate_gamma_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, gval, lrnr, slim = FALSE) {
  tt <- length(all_a)
  gamma_js <- map(1:tt, function(t) {
    if (t == 1) {
      estimate_gamma_j(data, folds, id, x, g, all_a[1], NULL, NULL, NULL, gval, t, lrnr)
    } else {
      estimate_gamma_j(data, folds, id, x, g, all_a[t], all_a[t - 1], all_y[t - 1], all_s[1:(t - 1)], gval, t, lrnr)
    }
  })

  all_ids <- select(data, !!id)
  out_gamma <- left_join(all_ids, gamma_js[[1]])
  for (j in 2:tt) {
    out_gamma <- left_join(out_gamma, gamma_js[[j]], by = id)
  }
  if (slim) {
    return(out_gamma)
  } else {
    return(left_join(data, out_gamma, by = id))
  }
}
