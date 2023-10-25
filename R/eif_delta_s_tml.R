eif_delta_s_tml <- function(y, a, g, e, gamma0, Q0, Qstar0, gamma1, Q1, Qstar1, pi, pistar) {
  tt <- ncol(y)
  pi1 <- ifelse(pi == 0, 1, pi)
  pi0 <- ifelse(pi == 1, 1, 1-pi)
  part_sj1 <- map(1:tt, ~eif_delta_tml_part_sj(., y, a, gamma1, Q1, Qstar1, pi1, pistar, tt)) %>%
    bind_cols
  part_sj0 <- map(1:tt, ~eif_delta_tml_part_sj(., y, a, gamma0, Q0, Qstar0, pi0, 1-pistar, tt)) %>%
    bind_cols
  part_yj1 <- map(1:tt, ~eif_delta_tml_part_yj(., y, a, gamma1, Q1, Qstar1, pi1, pistar, tt)) %>%
    bind_cols
  part_yj0 <- map(1:tt, ~eif_delta_tml_part_yj(., y, a, gamma0, Q0, Qstar0, pi0, 1-pistar, tt)) %>%
    bind_cols
  g/e*rowSums(part_yj1) - (1-g)/(1-e)*rowSums(part_yj0) + 1/e*rowSums(part_sj1) - 1/(1-e)*rowSums(part_sj0) + Q1[,1] - Q0[,1]
}

eif_delta_tml_part_yj <- function(j, y, a, gamma, Q, Qstar, pi, pistar, tt) {
  if (j == 1) {
    H_yj <- 1
  } else {
      ay_gamma_pi <- a*y/gamma*pistar/pi
      H_yj <- apply(ay_gamma_pi[,1:(j-1),drop=FALSE], 1, prod)
  }
  apq_gp_j <- a[,j]*pistar[,j]/gamma[,j]/pi[,j]
  if (j == tt) {
    out <- tibble(!!glue('part{j}') := H_yj*apq_gp_j*(y[,j] - Q[,j]))
  } else {
    out <- tibble(!!glue('part{j}') := H_yj*apq_gp_j*(y[,j]*Qstar[,j] - Q[,j]))
  }
  if (any(is.na(out))) browser()
  out
}


eif_delta_tml_part_sj <- function(j, y, a, gamma, Q, Qstar, pi, pistar, tt) {
  pi <- ifelse(pistar == 0 & pi == 0, 1, pi)
  pi <- ifelse(pistar == 1 & pi == 0, 1, pi)
  pistar <- ifelse(pistar == 0 & pi == 0, 1, pistar)
  pistar <- ifelse(pistar == 1 & pi == 0, 1, pistar)
  aypistar_gammapi <- a*y*pistar/gamma/pi
  ayp_gp <- apply(aypistar_gammapi[,1:j,drop=FALSE], 1, prod)

  if (j == tt) {
    out <- tibble(!!glue('part{j}') := 0)
  } else {
    out <- tibble(!!glue('part{j}') := ayp_gp*pistar[,j+1]*(Q[,j+1] - Qstar[,j]))
  }
  if (any(is.na(out))) browser()
  out
}
