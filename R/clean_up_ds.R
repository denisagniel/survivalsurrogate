clean_up_ds <- function(data, a = NULL, y, truncate_e = 0, truncate_pi = 0,
                        e = 'e', pistar = 'pistar', pi1 = 'pi1', pi0 = 'pi0',
                        mu1 = 'mu1', mu0 = 'mu0') {
  clean_data <- mutate_at(data, vars(any_of(c(a, y)), contains(pistar),
                                     contains(mu1),
                                     contains(mu0),
                                     contains('Q1'),
                                     contains('Q0'),
                                     contains('Qstar1'),
                                     contains('Qstar0'),
                                     contains('dmix'),
                                     contains('d_1'), contains('d_0')
                                              ), ~ifelse(is.na(.), 0, .))
  clean_data <- mutate_at(clean_data, vars(contains(pi1),
                                           contains(pi0),
                                              contains('gamma1'),
                                              contains('gamma0')),
                          ~ifelse(is.na(.), 1, .))

  if (truncate_pi > 0) {
    if (truncate_pi > 1) stop('Truncation point `truncate_pi` must be less than 1.')
    clean_data <- mutate_at(clean_data, vars(contains(pi1),
                                             contains(pi0),
                                             contains(pistar),
                                             contains('d_1'),
                                             contains('d_0'),
                                             contains('dmix')),
                            ~case_when(. < truncate_pi ~ truncate_pi,
                                       . > 1 - truncate_pi ~ 1 - truncate_pi,
                                       TRUE ~ .))
  }
  if (truncate_e > 0) {
    if (truncate_e > 1) stop('Truncation point `truncate_e` must be less than 1.')
    clean_data <- mutate_at(clean_data,
                            vars(contains('gamma'), e),
                         ~case_when(. < truncate_e ~ truncate_e,
                                       . > 1 - truncate_e ~ 1 - truncate_e,
                                       TRUE ~ .))
  }
  clean_data
}
