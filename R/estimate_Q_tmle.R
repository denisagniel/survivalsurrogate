estimate_Qj_tmle <- function(data, folds, id, x, g, a_jm1, a_j, y_jm1, y_j, sbar_jm1, Q_jp1, gammabar_j, e, gval, j, lrnr, epsilon = 1e-7) {
  at_risk_data <- data
  if (!is.null(a_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_jm1) == 1)
  }
  if (length(sbar_jm1) > 0) {
    s_jm1 <- sbar_jm1[length(sbar_jm1)]
    at_risk_data <- filter(at_risk_data, !is.na(!!sym(s_jm1)))
  }

    if (!is.null(a_j)) {
      at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval & !!sym(a_j) == 1)
    } else {
      at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(g) == gval,
                             a_j = 1)
      a_j <- 'a_j'
    }


  at_risk_data <- mutate(at_risk_data, Q_y = ifelse(!!sym(y_j) == 0, 0, !!sym(Q_jp1)))
  at_risk_data <- mutate(at_risk_data, Q_y = ifelse(!include_in_training, 0, Q_y))
  Q_nm <- paste0('Q', gval, '_', j)
  if (any(is.na(at_risk_data$Q_y))) browser()
  Q_j <- estimate_cont(at_risk_data, folds, id, c(x, sbar_jm1), 'Q_y', 'include_in_training', lrnr, Q_nm)
  updated_data <- inner_join(at_risk_data, Q_j, by = id)

  updated_data <- mutate(rowwise(updated_data), wt_j = ifelse(gval == 1, !!sym(g)/!!sym(e)*!!sym(a_j)/prod(c_across(gammabar_j)), (1-!!sym(g))/(1-!!sym(e))*!!sym(a_j)/prod(c_across(gammabar_j))),
                         !!Q_nm := pmax(pmin(!!sym(Q_nm), 1-epsilon), epsilon)) ### this is a kluge
  updated_data <- ungroup(updated_data)
  tmle_fm <- glue::glue('Q_y ~ offset(qlogis({Q_nm}))')
  tmle_fit <- glm(tmle_fm, weights = wt_j, data = updated_data %>%
                    filter(include_in_training), family = binomial)
  updated_data <- mutate(updated_data, !!Q_nm := plogis(qlogis(!!sym(Q_nm)) + coef(tmle_fit)))

  out_Q <- select(updated_data, !!id, !!Q_nm)
  out_Q <- left_join(data, out_Q, by = id)
  out_Q
}

estimate_Q_tmle <- function(data, folds, id, x, g, all_a, all_y, all_s, all_gamma, e, gval, lrnr, slim = FALSE) {
  tt <- length(all_y)
  tt_s <- length(all_s)
  Q_dat <- mutate(data, !!paste0('Q', gval, '_', tt+1) := 1)

  for (t in tt:1) {
    if (t == 0) {
      Q_dat <- estimate_Qj_tmle(Q_dat, folds, id, x, g, NULL, NULL, NULL, NULL, NULL, paste0('Q', gval, '_1'), NULL, e, gval, t, lrnr)
    } else if (t == 1) {
      Q_dat <- estimate_Qj_tmle(Q_dat, folds, id, x, g, NULL, all_a[1], NULL, all_y[1], NULL, paste0('Q', gval, '_', t+1), NULL, e, gval, t, lrnr)
    } else {
      Q_dat <- estimate_Qj_tmle(data = Q_dat,
                                folds = folds,
                                id = id,
                                x = x,
                                g = g,
                                a_jm1 = all_a[t-1],
                                a_j = all_a[t],
                                y_jm1 = all_y[t-1],
                                y_j = all_y[t],
                                sbar_jm1 = all_s[1:min(tt_s, t-1)],
                                Q_jp1 = paste0('Q', gval, '_', t+1),
                                gammabar_j = all_gamma[1:t],
                                e = e,
                                gval = gval,
                                j = t,
                                lrnr = lrnr)
    }
  }
  if (slim) {
    return(select(Q_dat, !!id, paste0('Q', gval, '_', 1:tt)))
  } else
    Q_dat
}
