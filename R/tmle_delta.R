#' Estimate Overall Treatment Effect using Targeted Minimum Loss-Based Estimation
#'
#' This function estimates the overall treatment effect using targeted minimum loss-based estimation (TMLE). TMLE is a statistical method that combines propensity score estimation with outcome regression to estimate treatment effects in observational studies.
#'
#' @param data A data frame containing the variables used in the analysis.
#' @param folds A vector specifying the fold assignment for cross-fitting.
#' @param id A character vector specifying the column name in the data frame that represents the unique identifier for each observation.
#' @param x A character vector specifying the column name(s) in the data frame that represent the pre-treatment covariate(s) used in the propensity score estimation and outcome regression.
#' @param g A character vector specifying the column name in the data frame that represents the treatment indicator variable.
#' @param a An optional character vector specifying the column name(s) in the data frame that represents an indicator of remaining uncensored. If not provided, all non-censoring indicators are assumed to be 1 (treatment group).
#' @param y A character vector specifying the column names in the data frame that represents the indicators of not having experienced the outcome.
#' @param s A character vector specifying the column name(s) in the data frame that represents the longitudinal surrogate marker.
#' @param binary_lrnr An optional learner object used for binary regression. If not provided, a default learner will be used.
#' @param cont_lrnr An optional learner object used for continuous regression. If not provided, a default learner will be used.
#' @param e An optional character vector specifying the column name(s) in the data frame that represent the propensity score. If not provided, the propensity scores will be estimated using the `estimate_propensity` function.
#' @param gamma1 An optional character vector specifying the column name(s) in the data frame that represent the censoring probabilities under treatment. If not provided, the censoring probabilities under treatment will be estimated using the `estimate_gamma_mat` function.
#' @param gamma0 An optional character vector specifying the column name(s) in the data frame that represent the censoring probabilities under control. If not provided, the censoring probabilities under control will be estimated using the `estimate_gamma_mat` function.
#' @param truncate_e A numeric value specifying the truncation threshold for propensity scores. Propensity scores below this threshold will be truncated to avoid numerical instability. Default is 1e-12.
#' @param verbose A logical value indicating whether to print additional information during the estimation process. Default is FALSE.
#'
#' @return A tibble containing the estimated treatment effect (`tmle_est`), the standard error of the estimated treatment effect (`tmle_se`), and the influence function data used for estimation (`if_data`).
#'
#' @examples
#' data <- read.csv("data.csv")
#' result <- tmle_delta(data = data,
#' folds = folds,
#' id = "id",
#' x = c("covariate1", "covariate2"),
#' g = "treatment",
#' a = c("noncensored1", "noncensored2", "noncensored3"),
#' y = c("outcome1", "outcome2", "outcome3")
#' s = c("surrogate1", "surrogate2"))
#'
#' @references
#' 1. van der Laan MJ, Rose S. Targeted Learning: Causal Inference for Observational and Experimental Data. Springer Science & Business Media; 2011.
#' 2. Gruber S, van der Laan MJ. tmle: An R package for targeted maximum likelihood estimation. J Stat Softw. 2010;35(13):1-33.
#' 3. Schuler MS, Rose S. Targeted maximum likelihood estimation for causal inference in observational studies. Am J Epidemiol. 2017;185(1):65-73.
tmle_delta <- function(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, truncate_e = 1e-12, verbose = FALSE) {
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
  analysis_data <- estimate_Q_tmle(data = analysis_data,
                                   folds = folds,
                                   id = id,
                                   x = x,
                                   g = g,
                                   all_a = a,
                                   all_y = y,
                                   all_s = s,
                                   all_gamma = gamma1,
                                   e = e,
                                   gval = 1,
                                   lrnr_c = cont_lrnr,
                                   lrnr_b = binary_lrnr)
  Q1 <- paste0('Q1_', 1:tt)
  analysis_data <- estimate_Q_tmle(data = analysis_data,
                                   folds = folds,
                                   id = id,
                                   x = x,
                                   g = g,
                                   all_a = a,
                                   all_y = y,
                                   all_s = s,
                                   all_gamma = gamma0,
                                   e = e,
                                   gval = 0,
                                   lrnr_c = cont_lrnr,
                                   lrnr_b = binary_lrnr)
  Q0 <- paste0('Q0_', 1:tt)
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


  if_ds <- transmute(analysis_data,
                     !!id := id,
                     Q1_1,
                     Q0_1,
                     eif = eif_delta_tml(y = y_m,
                                     a = a_m,
                                     g = !!sym(g),
                                     e = !!sym(e),
                                     gamma0 = gamma0_m,
                                     Q0 = Q0_m,
                                     gamma1 = gamma1_m,
                                     Q1 = Q1_m))

  summarise(if_ds,
            tmle_est = mean(Q1_1 - Q0_1),
            tmle_se = sd(eif)/sqrt(n()),
            if_data = list(if_ds))
}
