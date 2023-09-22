tmle_delta_s <- function(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, pi = NULL, pistar = NULL, truncate_e = 1e-12, verbose = FALSE) {
  tt <- length(y)
  analysis_data <- data
  if (all(is.null(e))) {
    if (verbose) {
      print('Propensity scores not provided in `e`. Estimating them.')
    }
    analysis_data <- estimate_propensity(data = analysis_data,
                                         folds = folds,
                                         id = id,
                                         x = x,
                                         g = g,
                                         lrnr = binary_lrnr)
    e <- 'e'
  }
  if (all(is.null(pi))) {
    if (verbose) {
      print('Pis not provided in `pi`. Estimating them.')
    }
    analysis_data <- estimate_pi_mat(data = analysis_data,
                                     folds = folds,
                                     id = id,
                                     x = x,
                                     g = g,
                                     all_a = a,
                                     all_y = y,
                                     all_s = s,
                                     lrnr = binary_lrnr)
    pi <- paste0('pi1_', 1:tt)
  }
  if (all(is.null(pistar))) {
    if (verbose) {
      print('Pistars not provided in `pistar`. Estimating them.')
    }
    analysis_data <- estimate_pistar_mat(data = analysis_data,
                                         folds = folds,
                                         id = id,
                                         x = x,
                                         g = g,
                                         all_a = a,
                                         all_y = y,
                                         all_s = s,
                                         lrnr = binary_lrnr)
    pistar <- paste0('pistar_', 1:tt)
  }
  if (all(is.null(gamma1))) {
    if (!all(is.null(a))) {
      if (verbose) {
        print('Censoring probabilities under treatment not provided in `gamma1`. Estimating them.')
      }
      analysis_data <- estimate_gamma_mat(data = analysis_data,
                                          folds = folds,
                                          id = id,
                                          x = x,
                                          g = g,
                                          all_a = a,
                                          all_y = y,
                                          all_s = s,
                                          gval = 1,
                                          lrnr = binary_lrnr)
    } else {
      warning('No censoring information provided. Assuming no censoring.')
      for (t in 1:tt) {
        analysis_data <- mutate(analysis_data, !!glue('gamma1_{t}') := 1)
      }
    }
    gamma1 <- paste0('gamma1_', 1:tt)
  }
  if (all(is.null(gamma0))) {
    if (!all(is.null(a))) {
      if (verbose) {
        print('Censoring probabilities under control not provided in `gamma0`. Estimating them.')
      }
      analysis_data <- estimate_gamma_mat(data = analysis_data,
                                          folds = folds,
                                          id = id,
                                          x = x,
                                          g = g,
                                          all_a = a,
                                          all_y = y,
                                          all_s = s,
                                          gval = 0,
                                          lrnr = binary_lrnr)
    } else {
      for (t in 1:tt) {
        analysis_data <- mutate(analysis_data, !!glue('gamma0_{t}') := 1)
      }
    }
    gamma0 <- paste0('gamma0_', 1:tt)
  }
  # browser()
  analysis_data <- estimate_Qstar_tmle(data = analysis_data,
                                   folds = folds,
                                   id = id,
                                   x = x,
                                   g = g,
                                   all_a = a,
                                   all_y = y,
                                   all_s = s,
                                   all_gamma = gamma1,
                                   pistar = pistar,
                                   pi = pi,
                                   e = e,
                                   gval = 1,
                                   lrnr_c = cont_lrnr,
                                   lrnr_b = binary_lrnr)
  Q1 <- paste0('Q1_', 1:tt)
  Qstar1 <- paste0('Qstar1_', 1:tt)
  analysis_data <- estimate_Qstar_tmle(data = analysis_data,
                                   folds = folds,
                                   id = id,
                                   x = x,
                                   g = g,
                                   all_a = a,
                                   all_y = y,
                                   all_s = s,
                                   all_gamma = gamma0,
                                   pistar = pistar,
                                   pi = pi,
                                   e = e,
                                   gval = 0,
                                   lrnr_c = cont_lrnr,
                                   lrnr_b = binary_lrnr)
  Q0 <- paste0('Q0_', 1:tt)
  Qstar0 <- paste0('Qstar0_', 1:tt)
  analysis_data <- clean_up_ds(analysis_data, a, y,
                               e = e, pistar = pistar, pi1 = pi,
                               truncate_e = truncate_e)


  gamma1_m <- ds_to_matrix(analysis_data, gamma1)
  gamma0_m <- ds_to_matrix(analysis_data, gamma0)
  y_m <- ds_to_matrix(analysis_data, y)
  if (all(is.null(a))) {
    a_m <- matrix(1, nrow(y_m), ncol(y_m))
  } else {
    a_m <- ds_to_matrix(analysis_data, a)
  }
  Q0_m <- ds_to_matrix(analysis_data, Q0)
  Q1_m <- ds_to_matrix(analysis_data, Q1)
  Qstar0_m <- ds_to_matrix(analysis_data, Qstar0)
  Qstar1_m <- ds_to_matrix(analysis_data, Qstar1)


  pi_m <- ds_to_matrix(analysis_data, pi)
  pistar_m <- ds_to_matrix(analysis_data, pistar)

  if_ds <- transmute(analysis_data,
                     !!id := id,
                     Q1_1,
                     Q0_1,
                     eif = eif_delta_s_tml(y = y_m,
                                         a = a_m,
                                         g = !!sym(g),
                                         e = !!sym(e),
                                         gamma0 = gamma0_m,
                                         Q0 = Q0_m,
                                         Qstar0 = Qstar0_m,
                                         gamma1 = gamma1_m,
                                         Q1 = Q1_m,
                                         Qstar1 = Qstar1_m,
                                         pi = pi_m,
                                         pistar = pistar_m
                                         ))

  summarise(if_ds,
            tmle_est = mean(Q1_1 - Q0_1),
            tmle_se = sd(eif)/sqrt(n()),
            if_data = list(if_ds))
}
