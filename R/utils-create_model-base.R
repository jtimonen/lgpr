# Create common Stan input needed for all models
create_model.base <- function(formula, data, prior, options, prior_only,
                              verbose) {

  # Data, formula and common Stan inputs
  data <- convert_to_data_frame(data)
  lgp_formula <- create_lgpformula(formula, data, verbose)
  parsed_input <- standata_common(
    data, lgp_formula, options, prior, prior_only, verbose
  )

  # Variable names
  var_names <- list(
    y = lgp_formula@y_name,
    x = rownames(dollar(parsed_input, "X")),
    z = rownames(dollar(parsed_input, "Z"))
  )

  # Create the 'lgpmodel' object
  new("lgpmodel",
    model_formula = lgp_formula,
    data = data,
    parsed_input = parsed_input,
    var_names = var_names,
    info = creation_info()
  )
}

# Checks if formula is in advanced format and translates if not
create_lgpformula <- function(formula, data, verbose = FALSE) {
  advanced <- is_advanced_formula(formula)
  if (!advanced) formula <- formula_to_advanced(formula, data, verbose)
  fp <- as.character(formula)
  formula_str <- paste(fp[2], fp[1], fp[3])
  log_info(paste0("Formula interpreted as: ", formula_str), verbose)
  parse_formula_advanced(formula)
}

# Create common Stan input needed for all models
standata_common <- function(data, lgp_formula, opts, prior, prior_only, vrb) {
  opts <- standata_common_options(opts, prior_only, vrb)
  covs <- standata_covariates(data, lgp_formula)
  comps <- standata_components(lgp_formula, covs)
  expanding <- standata_expanding(covs, comps)
  si <- c(opts, covs, comps, expanding)
  pri <- standata_common_prior(prior, si, vrb)
  lst <- c(si, pri)
  return(lst)
}

# Misc info for created objects
creation_info <- function() {
  list(
    created = date(),
    lgpr_version = utils::packageVersion("lgpr")
  )
}
