plugin_delta <- function(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, mu1 = NULL, mu0 = NULL, Q1 = NULL, Q0 = NULL, truncate_e = 1e-12, verbose = FALSE) {
if(verbose) {
  tt <- length(y)
  if (all(is.null(mu1))) {
    if (verbose) {
      cat('Hazards under treatment not provided in `mu1`. Estimating them.')
    }
    analysis_data <- estimate_mu_mat(data = data,
                                folds = folds,
                                id = id,
                                x = x,
                                g = g,
                                all_a = a,
                                all_y = y,
                                all_s = s,
                                gval = 1,
                                lrnr = binary_lrnr)
    mu1 <- paste0('mu1_', 1:tt)
  } else analysis_data <- data
  if (all(is.null(mu0))) {
    if (verbose) {
      cat('Hazards under control not provided in `mu0`. Estimating them.')
    }
    analysis_data <- estimate_mu_mat(data = analysis_data,
                                folds = folds,
                                id = id,
                                x = x,
                                g = g,
                                all_a = a,
                                all_y = y,
                                all_s = s,
                                gval = 0,
                                lrnr = binary_lrnr)
    mu0 <- paste0('mu0_', 1:tt)
  }
  if (all(is.null(Q1))) {
    if (verbose) {
      cat('Q functions under treatment not provided in `Q1`. Estimating them.')
    }
    analysis_data <- estimate_Q_mat(data = analysis_data,
                             folds = folds,
                             id = id,
                             x = x,
                             g = g,
                             all_a = a,
                             all_y = y,
                             all_s = s,
                             all_mu = mu1,
                             gval = 1,
                             lrnr = cont_lrnr)
    Q1 <- paste0('Q1_', 1:tt)
  }
  if (all(is.null(Q0))) {
    if (verbose) {
      cat('Q functions under control not provided in `Q0`. Estimating them.')
    }
    analysis_data <- estimate_Q_mat(data = analysis_data,
                             folds = folds,
                             id = id,
                             x = x,
                             g = g,
                             all_a = a,
                             all_y = y,
                             all_s = s,
                             all_mu = mu0,
                             gval = 0,
                             lrnr = cont_lrnr)
    Q0 <- paste0('Q0_', 1:tt)
  }
  if (all(is.null(e))) {
    if (verbose) {
      cat('Propensity scores not provided in `e`. Estimating them.')
    }
    analysis_data <- estimate_propensity(data = analysis_data,
                                         folds = folds,
                                         id = id,
                                         x = x,
                                         g = g,
                                         lrnr = binary_lrnr)
    e <- 'e'
  }

  if (all(is.null(gamma1))) {
    if (!all(is.null(a))) {
      if (verbose) {
        cat('Censoring probabilities under treatment not provided in `gamma1`. Estimating them.')
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
        cat('Censoring probabilities under control not provided in `gamma0`. Estimating them.')
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

  analysis_data <- clean_up_ds(analysis_data, a, y,
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
  mu1_m <- ds_to_matrix(analysis_data, mu1)
  mu0_m <- ds_to_matrix(analysis_data, mu0)


  if_ds <- transmute(analysis_data,
                     !!id := id,
                          eif = eif_delta(y = y_m,
                                          a = a_m,
                                          g = !!sym(g),
                                          e = !!sym(e),
                                          gamma0 = gamma0_m,
                                          mu0 = mu0_m,
                                          Q0 = Q0_m,
                                          gamma1 = gamma1_m,
                                          mu1 = mu1_m,
                                          Q1 = Q1_m))
  return(summarise(if_ds,
            plugin_est = mean(eif),
            plugin_se = sd(eif)/sqrt(n()),
            if_data = list(if_ds)))

}
if(!verbose){
	invisible(capture.output(
  suppressMessages({

   tt <- length(y)
  if (all(is.null(mu1))) {
    if (verbose) {
      cat('Hazards under treatment not provided in `mu1`. Estimating them.')
    }
    analysis_data <- estimate_mu_mat(data = data,
                                folds = folds,
                                id = id,
                                x = x,
                                g = g,
                                all_a = a,
                                all_y = y,
                                all_s = s,
                                gval = 1,
                                lrnr = binary_lrnr)
    mu1 <- paste0('mu1_', 1:tt)
  } else analysis_data <- data
  if (all(is.null(mu0))) {
    if (verbose) {
      cat('Hazards under control not provided in `mu0`. Estimating them.')
    }
    analysis_data <- estimate_mu_mat(data = analysis_data,
                                folds = folds,
                                id = id,
                                x = x,
                                g = g,
                                all_a = a,
                                all_y = y,
                                all_s = s,
                                gval = 0,
                                lrnr = binary_lrnr)
    mu0 <- paste0('mu0_', 1:tt)
  }
  if (all(is.null(Q1))) {
    if (verbose) {
      cat('Q functions under treatment not provided in `Q1`. Estimating them.')
    }
    analysis_data <- estimate_Q_mat(data = analysis_data,
                             folds = folds,
                             id = id,
                             x = x,
                             g = g,
                             all_a = a,
                             all_y = y,
                             all_s = s,
                             all_mu = mu1,
                             gval = 1,
                             lrnr = cont_lrnr)
    Q1 <- paste0('Q1_', 1:tt)
  }
  if (all(is.null(Q0))) {
    if (verbose) {
      cat('Q functions under control not provided in `Q0`. Estimating them.')
    }
    analysis_data <- estimate_Q_mat(data = analysis_data,
                             folds = folds,
                             id = id,
                             x = x,
                             g = g,
                             all_a = a,
                             all_y = y,
                             all_s = s,
                             all_mu = mu0,
                             gval = 0,
                             lrnr = cont_lrnr)
    Q0 <- paste0('Q0_', 1:tt)
  }
  if (all(is.null(e))) {
    if (verbose) {
      cat('Propensity scores not provided in `e`. Estimating them.')
    }
    analysis_data <- estimate_propensity(data = analysis_data,
                                         folds = folds,
                                         id = id,
                                         x = x,
                                         g = g,
                                         lrnr = binary_lrnr)
    e <- 'e'
  }

  if (all(is.null(gamma1))) {
    if (!all(is.null(a))) {
      if (verbose) {
        cat('Censoring probabilities under treatment not provided in `gamma1`. Estimating them.')
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
        cat('Censoring probabilities under control not provided in `gamma0`. Estimating them.')
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

  analysis_data <- clean_up_ds(analysis_data, a, y,
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
  mu1_m <- ds_to_matrix(analysis_data, mu1)
  mu0_m <- ds_to_matrix(analysis_data, mu0)


  if_ds <- transmute(analysis_data,
                     !!id := id,
                          eif = eif_delta(y = y_m,
                                          a = a_m,
                                          g = !!sym(g),
                                          e = !!sym(e),
                                          gamma0 = gamma0_m,
                                          mu0 = mu0_m,
                                          Q0 = Q0_m,
                                          gamma1 = gamma1_m,
                                          mu1 = mu1_m,
                                          Q1 = Q1_m))
  return(summarise(if_ds,
            plugin_est = mean(eif),
            plugin_se = sd(eif)/sqrt(n()),
            if_data = list(if_ds)))
  })))
  }


}
