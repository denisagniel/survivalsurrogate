ipw_delta_s <- function(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, e = NULL, pi = NULL, pistar = NULL, gamma1 = NULL, gamma0 = NULL, truncate_pi = 1e-12, truncate_e = 1e-12, verbose = FALSE, d_mix = NULL, d_1 = NULL, d_0 = NULL) {
  analysis_data <- data
  tt <- length(y)
  tt_s <- length(s)
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
  if (all(is.null(pi)) & (any(is.null(d_mix)) | any(is.null(d_1)) | any(is.null(d_0)) )) {
    if (verbose) {
      print('Pis not provided in `pi` and surrogate densities not provided. Estimating pi.')
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
  if (all(is.null(pistar)) & (any(is.null(d_mix)) | any(is.null(d_1)) | any(is.null(d_0)) )) {
    if (verbose) {
      print('Pistars not provided in `pistar` and surrogate densities not provided. Estimating them.')
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

  # analysis_data <- mutate(analysis_data,
  #                         w1_1 = 1/gamma1_1,
  #                         w0_1 = 1/gamma0_1)
  # for (t in 2:(tt_s+1)) { ## do some error checking for this / compare it to the length of y
  #   w1_j <- glue('w1_{t}')
  #   w0_j <- glue('w0_{t}')
  #   pi_jm1 <- glue('pi_{t-1}')
  #   pistar_jm1 <- glue('pistar_{t-1}')
  #   gamma1_j <- glue('gamma1_{t}')
  #   gamma0_j <- glue('gamma0_{t}')
  #   analysis_data <- mutate(analysis_data,
  #                           !!w1_j := !!sym(pistar_jm1)/!!sym(pi_jm1)/!!sym(gamma1_j),
  #                           !!w0_j := (1-!!sym(pistar_jm1))/(1-!!sym(pi_jm1))/!!sym(gamma0_j))
  # }
  analysis_data <- clean_up_ds(analysis_data, a, y, e = e,
                               truncate_pi = truncate_pi,
                               truncate_e = truncate_e)
  # rowwise_data <- rowwise(analysis_data)
  #
  # rowwise_data <- mutate(rowwise_data,
  #                        ipw_num_1 = prod(c_across(any_of(c(g, a, y)) | contains('pistar'))),
  #                        ipw_denom_1 = prod(c_across(contains('gamma1') | starts_with('pi_')))

  if (!all(is.null(a))) {
    a_m <- ds_to_matrix(analysis_data, a)
  } else {
    a_m <- matrix(1, 1, 1)
  }
  gamma1_m <- ds_to_matrix(analysis_data, gamma1)
  gamma0_m <- ds_to_matrix(analysis_data, gamma0)

  if (!is.null(pistar) & !is.null(pi)) {
    pistar1_m <- ds_to_matrix(analysis_data, pistar)
    pistar0_m <- 1 - pistar1_m
    pi1_m <- ds_to_matrix(analysis_data, pi)
    pi0_m <- 1 - pi1_m
    pistar1_m[,1] <- pistar0_m[,1] <- pi1_m[,1] <- pi0_m[,1] <- 1


    analysis_data <- mutate(analysis_data,
                            ipw_if1 = !!sym(g)/!!sym(e)*!!sym(y[tt])*matprod(a_m)*matprod(pistar1_m) / matprod(pi1_m) / matprod(gamma1_m),
                            ipw_if0 = (1-!!sym(g))/(1-!!sym(e))*!!sym(y[tt])*matprod(a_m)*matprod(pistar0_m) / matprod(pi0_m) / matprod(gamma0_m),
                            ipw_if = ifelse(!!sym(g) == 1, ipw_if1, -ipw_if0)
    )
    if_ds <- select(analysis_data, !!id, ipw_if)
    summarise(analysis_data,
              ipw_delta_s = mean(ipw_if),
              ipw_se = sd(ipw_if)/sqrt(n()),
              if_data = list(if_ds))
  } else {
    if (any(is.null(d_mix)) | any(is.null(d_1)) | any(is.null(d_0)) ) {
      stop("Surrogate densities expected but not provided.")
    }
    d_mix_m <- ds_to_matrix(analysis_data, d_mix)
    d1_m <- ds_to_matrix(analysis_data, d_1)
    d0_m <- ds_to_matrix(analysis_data, d_0)
    analysis_data <- mutate(analysis_data,
                            ipw_if1 = !!sym(g)/!!sym(e)*!!sym(y[tt])*matprod(a_m)*matprod(d_mix_m)/matprod(d1_m) / matprod(gamma1_m),
                            ipw_if0 = (1-!!sym(g))/(1-!!sym(e))*!!sym(y[tt])*matprod(a_m)*matprod(d_mix_m)/matprod(d0_m) / matprod(gamma0_m),
                            ipw_if = ifelse(!!sym(g) == 1, ipw_if1, -ipw_if0)
    )
    if_ds <- select(analysis_data, !!id, ipw_if)
    summarise(analysis_data,
              ipw_delta_s = mean(ipw_if),
              ipw_se = sd(ipw_if)/sqrt(n()),
              if_data = list(if_ds))
  }





}
