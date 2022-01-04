estimate_propensity <- function(data, folds, id, x, g, lrnr, slim = FALSE) {
  data <- mutate(data, include_in_training = TRUE)
  ps_data <- estimate_binary(data, folds, x, g, 'include_in_training', lrnr, 'e')
  if (slim) {
    return(select(ps_data, !!id, e))
  } else return(ps_data)
}

