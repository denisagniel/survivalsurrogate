matprod <- function(m) apply(m, 1, prod)

ds_to_matrix <- function(ds, prefix, tt = NULL) {
  if (is.null(tt)) {
    as.matrix(select(ds, all_of(paste0(prefix))))
  } else {
    as.matrix(select(ds, all_of(paste0(prefix, 1:tt))))
  }
}
