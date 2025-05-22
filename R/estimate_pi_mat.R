estimate_pi_j <- function(data, folds, id, x, g, a_j, y_j, sbar_j, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_j)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_j) == 1)
  }
  if (!is.null(y_j)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_j) == 1)
  }
  at_risk_data <- mutate(at_risk_data, include_in_training = TRUE)
  at_risk_data <- estimate_binary(at_risk_data, folds, id, c(x, sbar_j), g, "include_in_training", lrnr, paste0("pi1_", j))
  at_risk_data <- mutate(at_risk_data, !!paste0("pi0_", j) := 1 - !!sym(paste0("pi1_", j)))
  at_risk_data
}

estimate_pi_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, lrnr, t0 = length(all_s), slim = FALSE) {
  tt <- length(all_y)
  pi_js <- map(1:t0, function(t) {
    estimate_pi_j(data, folds, id, x, g, all_a[t], all_y[t], all_s[1:min(t, tt - 1)], t + 1, lrnr)
    ### this is min(t, tt-1) because s does not go all the way through to the final timepoint, so the largest its index can be is tt-1
  })
  # pi_js_slim <- map(pi_js, ~select(., !!id, contains('pi1'), contains('pi0')))
  all_ids <- select(data, !!id)
  out_pi <- mutate(all_ids, pi1_1 = 1)
  for (j in 1:t0) {
    out_pi <- left_join(out_pi, pi_js[[j]], by = id)
  }
  if (tt > t0 + 1) {
    for (ttt in (t0 + 2):tt) {
      out_pi <- mutate(out_pi, !!paste0("pi1_", ttt) := 1)
    }
  }
  if (slim) {
    return(out_pi)
  } else {
    return(left_join(data, out_pi, by = id))
  }
}
