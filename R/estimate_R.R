estimate_R <- function(delta_if, delta_s_if, delta = NULL, delta_s = NULL, se_type = 'asymptotic', n_boot = NULL, alpha = 0.05) {
  if (is.null(delta)) delta  <- mean(delta_if)
  if (is.null(delta_s)) delta_s <- mean(delta_s_if)
  n <- length(delta_if)

  if (se_type == 'asymptotic') {
    sigma_sq <- delta^(-2)*mean(delta_if^2) + delta_s^2*delta^(-4)*mean(delta_s_if^2) -
      2*delta_s*delta^(-3)*mean(delta_if*delta_s_if)
    tibble(estimand = c('Delta', 'Delta_S', 'R'),
           estimate = c(delta, delta_s, 1 - delta_s/delta),
           se = c(sd(delta_if)/sqrt(n),
                  sd(delta_s_if)/sqrt(n),
                  sqrt(sigma_sq)/sqrt(n)),
           ci_l = estimate - qnorm(1 - alpha/2)*se,
           ci_h = estimate + qnorm(1 - alpha/2)*se)
  } else if (se_type == 'bootstrap') {
    if (is.null(n_boot)) stop('If se_type = "bootstrap", must provide number of bootstraps in n_boot.')
    # browser()
    gmat <- rBeta2009::rdirichlet(n_boot, rep(1, n))*n
    boot_res_l <- list()
    for (b in 1:n_boot) {
      wt_b <- gmat[b,]
      delta_b <- sum(delta_if*wt_b)/sum(wt_b)
      delta_s_b <- sum(delta_s_if*wt_b)/sum(wt_b)
      boot_res_l[[b]] <- tibble(estimand = c('Delta', 'Delta_S', 'R'),
                              boot_estimate = c(delta_b,
                                                delta_s_b,
                                                1 - delta_s_b/delta_b))
    }
    # browser()
    boot_res <- bind_rows(boot_res_l) %>%
      group_by(estimand) %>%
      summarise(estimate = unique(case_when(estimand == 'Delta' ~ delta,
                                     estimand == 'Delta_S' ~ delta_s,
                                     estimand == 'R' ~ 1 - delta_s/delta)),
                se = sd(boot_estimate),
                ci_l = quantile(boot_estimate, alpha/2),
                ci_h = quantile(boot_estimate, 1 - alpha/2))


    boot_res
  }

}
