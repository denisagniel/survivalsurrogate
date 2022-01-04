estimate_propensity <- function(data, folds, id, x, g, lrnr, slim = FALSE) {
  data <- mutate(data, include_in_training = TRUE)
  ps_data <- estimate_binary(data, folds, id, x, g, 'include_in_training', lrnr, 'e')
  if (slim) {
    return(ps_data)
  } else return(left_join(data, ps_data, by = id))
}

