gcomp_delta_s <- function(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, mu1 = NULL, mu0 = NULL, Q1 = NULL, Q0 = NULL, Qstar1 = NULL, Qstar0 = NULL, verbose = FALSE) {
  tt <- length(y)
  t0 <- length(s)
  if (all(is.null(mu1))) {
    if (verbose) {
      print('Hazards under treatent not provided in `mu1`. Estimating them.')
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
      print('Hazards under control not provided in `mu0`. Estimating them.')
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

  #############
  ## do something smarter with this - only need to estimate some Q's sometimes
  #############
  if (all(is.null(Q1))) {
    if (verbose) {
      print('Q functions under treatent not provided in `Q1`. Estimating them.')
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
      print('Q functions under treatent not provided in `Q1`. Estimating them.')
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

  if (all(is.null(Qstar1))) {
    if (verbose) {
      print('Qstar functions under treatment not provided in `Qstar1`. Estimating them.')
    }
    analysis_data <- estimate_Qstar_mat(data = analysis_data,
                                     folds = folds,
                                     id = id,
                                     x = x,
                                     g = g,
                                     all_a = a,
                                     all_y = y,
                                     all_s = s,
                                     all_Q = Q1,
                                     all_mu = mu1,
                                     gval = 1,
                                     lrnr = cont_lrnr)
    Qstar1 <- paste0('Qstar1_', 1:tt)
  }
  if (all(is.null(Qstar0))) {
    if (verbose) {
      print('Qstar functions under control not provided in `Qstar0`. Estimating them.')
    }
    analysis_data <- estimate_Qstar_mat(data = analysis_data,
                                     folds = folds,
                                     id = id,
                                     x = x,
                                     g = g,
                                     all_a = a,
                                     all_y = y,
                                     all_s = s,
                                     all_Q = Q0,
                                     all_mu = mu0,
                                     gval = 0,
                                     lrnr = cont_lrnr)
    Qstar0 <- paste0('Qstar0_', 1:tt)
  }
  if_ds <- transmute(analysis_data,
                     !!id,
                     gcomp_if = !!sym(Qstar1[1])*!!sym(mu1[1]) - !!sym(Qstar0[1])*!!sym(mu0[1]))
  summarise(if_ds,
            gcomp_est = mean(gcomp_if),
            gcomp_se = sd(gcomp_if)/sqrt(n()),
            if_data = list(if_ds))
}
