estimate_Qstar_yj_tmle <- function(data, folds, id, x, g, a_jm1, a_j, y_jm1, y_j, sbar_jm1, Q_jp1, gammabar_j, pistar_jm1, pi_jm1, e, gval, j, lrnr_b, lrnr_c, epsilon = 1e-7) {
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
  if (all(abs(at_risk_data$Q_y - 1) < 1e-3 | abs(at_risk_data$Q_y) < 1e-3)) {
    Q_j <- estimate_binary(at_risk_data, folds, id, c(x, sbar_jm1), 'Q_y', 'include_in_training', lrnr_b, Q_nm)
  } else {
    Q_j <- estimate_cont(at_risk_data, folds, id, c(x, sbar_jm1), 'Q_y', 'include_in_training', lrnr_c, Q_nm)
  }

  updated_data <- inner_join(at_risk_data, Q_j, by = id)
  updated_data <- mutate(updated_data,
    # rowwise(updated_data), wt_j = ifelse(gval == 1, !!sym(g)/!!sym(e)*!!sym(a_j)/prod(c_across(gammabar_j))*prod(c_across(pistar_jm1))/prod(c_across(pi_jm1)), (1-!!sym(g))/(1-!!sym(e))*!!sym(a_j)/prod(c_across(gammabar_j))*prod(c_across(pistar_jm1))/prod(c_across(pi_jm1))),
                         !!Q_nm := pmax(pmin(!!sym(Q_nm), 1-epsilon), epsilon)) ### this is a kluge

  updated_data <- mutate(updated_data,
         static_part = gval*!!sym(g)/!!sym(e)*!!sym(a_j) + (1-gval)*(1-!!sym(g))/(1-!!sym(e))*!!sym(a_j),
         gamma_part = 1,
         pistar_part = 1,
         pi_part = 1)
  for(gam in gammabar_j) {
    # print(gam)
    updated_data <- mutate(updated_data,
           gamma_part = gamma_part/!!sym(gam))
  }
  for(pist in pistar_jm1) {
    # print(pist)
    updated_data <- mutate(updated_data,
           pistar_part = pistar_part*!!sym(pist))
  }
  for(pi in pi_jm1) {
    # print(pi)
    updated_data <- mutate(updated_data,
           pi_part = pi_part/!!sym(pi))
  }
  updated_data <- mutate(updated_data, wt_j = static_part*gamma_part*pistar_part*pi_part)

  # browser()

  # updated_data <- ungroup(updated_data)
  tmle_fm <- glue::glue('Q_y ~ offset(qlogis({Q_nm}))')
  tmle_fit <- glm(tmle_fm, weights = wt_j, data = updated_data %>%
                    filter(include_in_training), family = binomial)
  updated_data <- mutate(updated_data, !!Q_nm := plogis(qlogis(!!sym(Q_nm)) + coef(tmle_fit)))

  out_Q <- select(updated_data, !!id, !!Q_nm)
  out_Q <- left_join(data, out_Q, by = id)
  out_Q
}

estimate_Qstar_sj_tmle <- function(data, folds, id, x, g, a_jm1, a_j, y_jm1, y_j, sbar_jm1, Q_jp1, gammabar_j, pistar_jm1, pi_jm2, e, gval, j, lrnr_b, lrnr_c, epsilon = 1e-7) {
  at_risk_data <- data
  if (!is.null(a_jm1)) {
    at_risk_data <- filter(at_risk_data, !!sym(a_jm1) == 1)
  }
  if (length(sbar_jm1) > 0) {
    s_jm1 <- sbar_jm1[length(sbar_jm1)]
    at_risk_data <- filter(at_risk_data, !is.na(!!sym(s_jm1)))
  }

  if (!is.null(a_j)) {
    at_risk_data <- mutate(at_risk_data, include_in_training = !!sym(a_j) == 1)
  } else {
    at_risk_data <- mutate(at_risk_data, include_in_training = !is.na(!!sym(Q_jp1)),
                           a_j = 1)
    a_j <- 'a_j'
  }


  at_risk_data <- mutate(at_risk_data, Q_s = !!sym(Q_jp1))
  at_risk_data <- mutate(at_risk_data, Q_s = ifelse(!include_in_training, 0, Q_s))
  Q_nm <- paste0('Qstar', gval, '_', j)
  if (any(is.na(at_risk_data$Q_s))) browser()
  if (all(abs(at_risk_data$Q_s - 1) < 1e-3 | abs(at_risk_data$Q_s) < 1e-3)) {
    Q_j <- estimate_binary(at_risk_data, folds, id, c(x, sbar_jm1), 'Q_s', 'include_in_training', lrnr_b, Q_nm)
  } else {
    Q_j <- estimate_cont(at_risk_data, folds, id, c(x, sbar_jm1), 'Q_s', 'include_in_training', lrnr_c, Q_nm)
  }

  updated_data <- inner_join(at_risk_data, Q_j, by = id)

  updated_data <- mutate(updated_data,
                         # rowwise(updated_data), wt_j = ifelse(gval == 1, 1/!!sym(e)*!!sym(a_j)/prod(c_across(gammabar_j))*prod(c_across(pistar_jm1))/prod(c_across(pi_jm2)), 1/(1-!!sym(e))*!!sym(a_j)/prod(c_across(gammabar_j))*prod(c_across(pistar_jm1))/prod(c_across(pi_jm2))),
                         !!Q_nm := pmax(pmin(!!sym(Q_nm), 1-epsilon), epsilon)) ### this is a kluge
  updated_data <- mutate(updated_data,
         static_part = gval/!!sym(e)*!!sym(a_j) + (1-gval)/(1-!!sym(e))*!!sym(a_j),
         gamma_part = 1,
         pistar_part = 1,
         pi_part = 1)
  for(gam in gammabar_j) {
    # print(gam)
    updated_data <- mutate(updated_data,
           gamma_part = gamma_part/!!sym(gam))
  }
  for(pist in pistar_jm1) {
    # print(pist)
    updated_data <- mutate(updated_data,
           pistar_part = pistar_part*!!sym(pist))
  }
  for(pi in pi_jm2) {
    # print(pi)
    updated_data <- mutate(updated_data,
           pi_part = pi_part/!!sym(pi))
  }
  updated_data <- mutate(updated_data, wt_j = static_part*gamma_part*pistar_part*pi_part)

  # browser()

  # updated_data <- ungroup(updated_data)
  tmle_fm <- glue::glue('Q_s ~ offset(qlogis({Q_nm}))')
  tmle_fit <- glm(tmle_fm, weights = wt_j, data = updated_data %>%
                    filter(include_in_training), family = binomial)
  updated_data <- mutate(updated_data, !!Q_nm := plogis(qlogis(!!sym(Q_nm)) + coef(tmle_fit)))

  out_Q <- select(updated_data, !!id, !!Q_nm)
  out_Q <- left_join(data, out_Q, by = id)
  out_Q
}

estimate_Qstarj_tmle <- function(data, folds, id, x, g, a_jm1, a_j, y_jm1, y_j, sbar_jm1, Q_jp1, gammabar_j, e, pistar_jm1, pi_jm1, pi_jm2, gval, j, lrnr_b, lrnr_c, epsilon = 1e-7) {
  data_j <- estimate_Qstar_sj_tmle(data = data,
                                   folds = folds,
                                   id = id,
                                   x = x,
                                   g = g,
                                   a_jm1 = a_jm1,
                                   a_j = a_j,
                                   y_jm1 = y_jm1,
                                   y_j = y_j,
                                   sbar_jm1 = sbar_jm1,
                                   Q_jp1 = Q_jp1,
                                   gammabar_j = gammabar_j,
                                   e = e,
                                   pistar_jm1= pistar_jm1,
                                   pi_jm2 = pi_jm2,
                                   gval = gval,
                                   j = j,
                                   lrnr_b = lrnr_b,
                                   lrnr_c = lrnr_c,
                                   epsilon = epsilon)
  out_j <- estimate_Qstar_yj_tmle(data = data_j,
                                  folds = folds,
                                  id = id,
                                  x = x,
                                  g = g,
                                  a_jm1 = a_jm1,
                                  a_j = a_j,
                                  y_jm1 = y_jm1,
                                  y_j = y_j,
                                  sbar_jm1 = sbar_jm1,
                                  Q_jp1 = paste0('Qstar', gval, '_', j),
                                  gammabar_j = gammabar_j,
                                  e = e,
                                  pistar_jm1= pistar_jm1,
                                  pi_jm1 = pi_jm1,
                                  gval = gval,
                                  j = j,
                                  lrnr_b = lrnr_b,
                                  lrnr_c = lrnr_c,
                                  epsilon = epsilon)
  out_j
}

estimate_Qstar_tmle <- function(data, folds, id, x, g, all_a, all_y, all_s, all_gamma, pistar, pi, e, gval, lrnr_b, lrnr_c, slim = FALSE) {
  tt <- length(all_y)
  tt_s <- length(all_s)
  Q_dat <- mutate(data, !!paste0('Q', gval, '_', tt+1) := 1)

  for (t in tt:1) {
    if (t == 0) {
      Q_dat <- estimate_Qstarj_tmle(data = Q_dat, folds = folds, id = id,
                                    x = x, g = g, a_jm1 = NULL, a_j = NULL,
                                    y_jm1 = NULL, y_j = NULL,
                                    sbar_jm1 = NULL, Q_jp1 = paste0('Q', gval, '_1'),
                                    gammabar_j = NULL, pistar_jm1 = NULL, pi_jm1 = NULL, pi_jm2 = NULL, e = e, gval = gval, j = t, lrnr_b = lrnr_b,
                                    lrnr_c = lrnr_c,)
    } else if (t == 1) {
      Q_dat <- estimate_Qstarj_tmle(data = Q_dat, folds = folds, id = id,
                                    x = x, g = g, a_jm1 = NULL, a_j = all_a[1], y_jm1 = NULL,
                                    y_j = all_y[1], sbar_jm1 = NULL,
                                    Q_jp1 = paste0('Q', gval, '_', t+1),
                                    gammabar_j = NULL, pistar_jm1 = NULL, pi_jm1 = NULL, pi_jm2 = NULL, e = e, gval = gval, j = t, lrnr_b = lrnr_b,
                                    lrnr_c = lrnr_c,)
    } else if (t == 2) {
      Q_dat <- estimate_Qstarj_tmle(data = Q_dat,
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
                                    pistar_jm1 = pistar[1:(t-1)],
                                    pi_jm1 = pi[1:(t-1)],
                                    pi_jm2 = NULL,
                                    e = e,
                                    gval = gval,
                                    j = t,
                                    lrnr_b = lrnr_b,
                                    lrnr_c = lrnr_c,)
      } else {
      Q_dat <- estimate_Qstarj_tmle(data = Q_dat,
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
                                pistar_jm1 = pistar[1:(t-1)],
                                pi_jm1 = pi[1:(t-1)],
                                pi_jm2 = pi[1:(t-2)],
                                e = e,
                                gval = gval,
                                j = t,
                                lrnr_b = lrnr_b,
                                lrnr_c = lrnr_c,)
    }
  }
  if (slim) {
    return(select(Q_dat, !!id, paste0('Q', gval, '_', 1:tt)))
  } else
    Q_dat
}
