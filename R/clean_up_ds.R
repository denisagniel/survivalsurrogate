clean_up_ds <- function(data, a = NULL, y, truncate_e = 0, truncate_pi = 0,
                        e = "e", pistar = "pistar", pi1 = "pi1", pi0 = "pi0",
                        mu1 = "mu1", mu0 = "mu0") {
  clean_data <- dplyr::mutate_at(data, dplyr::vars(
    dplyr::any_of(c(a, y, e)),
    dplyr::contains("Q1"),
    dplyr::contains("Q0"),
    dplyr::contains("Qstar1"),
    dplyr::contains("Qstar0")
  ), ~ ifelse(is.na(.), 0, .))


  if (length(pistar) == 1) {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::contains(pistar)), ~ ifelse(is.na(.), 0, .))
  } else {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::any_of(pistar)), ~ ifelse(is.na(.), 0, .))
  }
  if (length(mu1) == 1) {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::contains(mu1)), ~ ifelse(is.na(.), 0, .))
  } else {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::any_of(mu1)), ~ ifelse(is.na(.), 0, .))
  }
  if (length(mu0) == 1) {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::contains(mu0)), ~ ifelse(is.na(.), 0, .))
  } else {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::any_of(mu0)), ~ ifelse(is.na(.), 0, .))
  }

  if (length(pi1) == 1) {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::contains(pi1)), ~ ifelse(is.na(.), 1, .))
  } else {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::any_of(pi1)), ~ ifelse(is.na(.), 1, .))
  }
  if (length(pi0) == 1) {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::contains(pi0)), ~ ifelse(is.na(.), 1, .))
  } else {
    clean_data <- dplyr::mutate_at(clean_data, dplyr::vars(dplyr::any_of(pi0)), ~ ifelse(is.na(.), 1, .))
  }

  clean_data <- dplyr::mutate_at(
    clean_data, dplyr::vars(
      dplyr::contains("gamma1"),
      dplyr::contains("gamma0")
    ),
    ~ ifelse(is.na(.), 1, .)
  )

  if (truncate_pi > 0) {
    if (truncate_pi > 1) stop("Truncation point `truncate_pi` must be less than 1.")
    if (length(pi1) == 1) {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::contains(pi1)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    } else {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::any_of(pi1)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    }
    if (length(pi0) == 1) {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::contains(pi0)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    } else {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::any_of(pi0)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    }
    if (length(pistar) == 1) {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::contains(pistar)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    } else {
      clean_data <- dplyr::mutate_at(
        clean_data, dplyr::vars(dplyr::any_of(pi1)),
        ~ case_when(
          . < truncate_pi ~ truncate_pi,
          . > 1 - truncate_pi ~ 1 - truncate_pi,
          TRUE ~ .
        )
      )
    }
  }
  if (truncate_e > 0) {
    if (truncate_e > 1) stop("Truncation point `truncate_e` must be less than 1.")
    clean_data <- dplyr::mutate_at(
      clean_data,
      dplyr::vars(dplyr::contains("gamma"), e),
      ~ case_when(
        . < truncate_e ~ truncate_e,
        . > 1 - truncate_e ~ 1 - truncate_e,
        TRUE ~ .
      )
    )
  }
  clean_data
}
