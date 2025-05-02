eif_delta_s <- function(y, a, g, e, pi, pistar, gamma0, mu0, Q0, gamma1, mu1, Q1, t0) {
  tt <- ncol(y)

  ## check that pi and pistar are the correct width
  ## check that pi and pistar are uniformly 1 after t0

  part_j1_mu <- map(1:tt, ~eif_delta_s_partj_mu(., y, a, gamma1, mu1, pi, pistar, Q1, t0)) %>%
    bind_cols
  part_j0_mu <- map(1:tt, ~eif_delta_s_partj_mu(., y, a, gamma0, mu0, 1-pi, 1-pistar, Q0, t0)) %>%
    bind_cols
  part_j1_q <- map(1:t0, ~eif_delta_s_partj_q(., y, a, gamma1, mu1, pi, pistar, Q1)) %>%
    bind_cols
  part_j0_q <- map(1:t0, ~eif_delta_s_partj_q(., y, a, gamma0, mu0, 1-pi, 1-pistar, Q0)) %>%
    bind_cols


  g/e*(rowSums(part_j1_mu)) -
    (1-g)/(1-e)*(rowSums(part_j0_mu)) +
    1/e*rowSums(part_j1_q) - 1/(1-e)*rowSums(part_j0_q) +
    mu1[,1]*Q1[,1] - mu0[,1]*Q0[,1]
}

eif_delta_s_partj_mu <- function(j, y, a, gamma, mu, pi, pistar, Q, t0) {
  pi <- ifelse(pistar == 0 & pi == 0, 1, pi)
  pi <- ifelse(pistar == 1 & pi == 0, 1, pi)
  pi <- ifelse(pistar == 0 & pi == 1, 1, pi)
  pistar <- ifelse(pistar == 0 & pi == 0, 1, pistar)
  pistar <- ifelse(pistar == 1 & pi == 0, 1, pistar)
  pistar <- ifelse(pistar == 0 & pi == 1, 1, pistar)
  aypistar_gammapi <- pistar/gamma/pi*a*y
  if (j == 1) {
    ayp_gp_jm1 <- 1
  } else {
    ayp_gp_jm1 <- apply(aypistar_gammapi[,1:(j-1),drop = FALSE], 1, prod)
  }

  apq_gp_j <- a[,j]*pistar[,j]/gamma[,j]/pi[,j]
  # Qtilde <- Q
  # tt <- ncol(y)
  # if (tt >= t0+2) {
  #   Qtilde[,(t0+1):(tt-1)] <- mu[,(t0+2):tt]*Q[,(t0+2):tt]
  # }
  # Qtilde[,tt] <- 1

  out <- tibble(!!glue('mu_part{j}') := ayp_gp_jm1*apq_gp_j*(y[,j]*Q[,j] - mu[,j]*Q[,j]))
  if (any(is.na(out))) browser()
  out
}

eif_delta_s_partj_q <- function(j, y, a, gamma, mu, pi, pistar, Q) {
  pi <- ifelse(pistar == 0 & pi == 0, 1, pi)
  pi <- ifelse(pistar == 1 & pi == 0, 1, pi)
  pistar <- ifelse(pistar == 0 & pi == 0, 1, pistar)
  pistar <- ifelse(pistar == 1 & pi == 0, 1, pistar)
  aypistar_gammapi <- a*y*pistar/gamma/pi
  ayp_gp <- apply(aypistar_gammapi[,1:j,drop=FALSE], 1, prod)

  out <- tibble(!!glue('Q_part{j}') := ayp_gp*pistar[,j+1]*(mu[,j+1]*Q[,j+1] - Q[,j]))
  if (any(is.na(out))) browser()
  out
}

