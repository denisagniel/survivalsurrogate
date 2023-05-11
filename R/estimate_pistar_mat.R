
estimate_pistar_mat <- function(data, folds, id, x, g, all_a, all_y, all_s, lrnr, t0 = length(all_s), slim = FALSE) {
  tt <- length(all_y)
  pistar_js <- map(1:t0, function(t) {
    if (t == 1) {
      estimate_pistar_j(data, folds, id, x, g, all_a[t],  all_y[t], NULL, t+1, lrnr)
    } else {
      estimate_pistar_j(data, folds, id, x, g, all_a[t],  all_y[t], all_s[1:(t-1)], t+1, lrnr)
    }

  })
  # pistar_js_slim <- map(pistar_js, ~select(., !!id, contains('pistar')))
  all_ids <- select(data, !!id)
  out_pistar <- mutate(all_ids, pistar_1 = 1)
  for (j in 1:t0) {
    out_pistar <- left_join(out_pistar, pistar_js[[j]], by = id)
  }
  if (tt > t0+1) {
    for (ttt in (t0+2):tt) {
      out_pistar <- mutate(out_pistar, !!paste0('pistar_', ttt) := 1)
    }
  }
  if (slim) {
    return(out_pistar)
  } else {
    return(left_join(data, out_pistar, by = id))
  }

}

estimate_pistar_j <- function(data, folds, id, x, g, a_j, y_j, sbar_j, j, lrnr) {
  at_risk_data <- data
  if (!is.null(a_j)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_j) == 1)
  }
  if (!is.null(y_j)) {
    at_risk_data <- filter(at_risk_data, !!sym(y_j) == 1)
  }
  at_risk_data <- mutate(at_risk_data, include_in_training = TRUE)
  estimate_binary(at_risk_data, folds, id, c(x, sbar_j), g, 'include_in_training', lrnr, paste0('pistar_', j))
}
