learn_fold_prob <- function(task, train_ids, test_ids, lrnr) {
  lrnr$train(task, row_ids = train_ids)
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  tibble(
    row_id = predicted_vals$row_ids,
    pred = predicted_vals$prob[, "1"]
  )
}

learn_fold_cont <- function(task, train_ids, test_ids, lrnr) {
  lrnr$train(task, row_ids = train_ids)
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  tibble(
    row_id = predicted_vals$row_ids,
    pred = predicted_vals$response
  )
}
