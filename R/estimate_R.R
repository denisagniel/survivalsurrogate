estimate_R <- function(delta_if, delta_s_if, delta = NULL, delta_s = NULL) {
  if (is.null(delta)) delta  <- mean(delta_if)
  if (is.null(delta_s)) delta_s <- mean(delta_s_if)
  n <- length(delta_if)

  sigma_sq <- delta^(-2)*mean(delta_if^2) + delta_s^2*delta^(-4)*mean(delta_s_if^2) -
    2*delta_s*delta^(-3)*mean(delta_if*delta_s_if)
  tibble(estimand = c('Delta', 'Delta_S', 'R'),
         estimate = c(delta, delta_s, 1 - delta_s/delta),
         se = c(sd(delta_if)/sqrt(n),
                sd(delta_s_if)/sqrt(n),
                sqrt(sigma_sq)/sqrt(n)))
}
