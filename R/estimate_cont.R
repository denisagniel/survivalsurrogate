estimate_cont <- function(data, folds, x, y, include_train, lrnr, task_name) {
  data <- mutate(data, row_id = 1:nrow(data))

  xy_dat <- select(data, any_of(c(x, y)))
  this_task <- mlr3::as_task_regr(xy_dat, target = y, id = task_name)

  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)
  incl <- pull(data, include_train)
  predictions <- map_df(unique_folds,
                        ~learn_fold_cont(task = this_task,
                                         train_ids = which(all_folds != . & incl),
                                         test_ids = which(all_folds == .),
                                         lrnr = lrnr))
  predictions <- rename(predictions, !!task_name := pred)
  inner_join(data, predictions)
}
