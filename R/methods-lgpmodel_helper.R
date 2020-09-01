#' @rdname model_getters
get_covariate_info_cont <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- vn$x_cont
  num_nan <- rowSums(model@stan_input$x_cont_mask == 1)
  df <- data.frame(nam, num_nan)
  colnames(df) <- c("Variable", "#Missing")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}

#' @rdname model_getters
get_covariate_info_cat <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- vn$x_cat
  num_levels <- model@stan_input$x_cat_num_levels
  levels <- model@var_info$x_cat_levels
  level_names <- c()
  J <- length(nam)
  for (j in seq_len(J)) {
    a <- levels[[j]]
    if (length(a) > 4) {
      a <- "..."
    } else {
      a <- paste(a, collapse = ", ")
    }
    level_names[j] <- a
  }
  df <- data.frame(nam, num_levels, level_names)
  colnames(df) <- c("Factor", "#Levels", "Values")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}
