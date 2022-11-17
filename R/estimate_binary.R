estimate_binary <- function(data, folds, id, x, y, include_train, lrnr, task_name) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))

  xy_dat <- select(data, any_of(c(x, y)))
  if (!inherits(pull(xy_dat, !!y), 'factor')) {
    xy_dat <- mutate_at(xy_dat, vars(y), as.factor)
  }
  this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name)

  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)
  incl <- pull(data, include_train)
  predictions <- map_df(unique_folds,
                        ~learn_fold_prob(task = this_task,
                                         train_ids = which(all_folds != . & incl),
                                         test_ids = which(all_folds == .),
                                         lrnr = lrnr))
  predictions <- rename(predictions, !!task_name := pred)
  out_ds <- inner_join(data, predictions, by = 'row_id')
  select(out_ds, !!id, !!task_name)
}
