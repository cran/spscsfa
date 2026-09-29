#' Print SPSC-SFA results
#'
#' @param x An object of class `spsc_sfa`.
#' @param ... Additional arguments (currently ignored).
#' @return The original object of class \code{"spsc_sfa"}, returned invisibly to avoid duplicate printing
#'  after printing its summary. The returned object is a list containing the 
#'  components \code{smooth_coefficients}, \code{mle_results}, \code{technical_efficiency},
#'  \code{inefficiency}, \code{model_summary}, \code{formula_parsed}, and \code{data_info}.
#'  These components contain the estimated smooth coefficients and recovered
#'  intercept, maximum-likelihood results, technical-efficiency estimates and
#'  confidence intervals, JLMS inefficiency estimates, model-specification
#'  information, parsed formula components, and data-cleaning information.
#' @export
print.spsc_sfa <- function(x, ...) {
  cat("\n")
  cat("==================== SPSC-SFA RESULTS ====================\n")
  cat("\nFormula:", deparse(x$model_summary$formula), "\n")
  cat("Number of observations:", x$model_summary$n_obs, "\n")

  if (x$model_summary$data_supplied) {
    cat("\n--- Data Information ---\n")
    cat("Data frame supplied: YES\n")
    if (!is.null(x$data_info)) {
      cat("Rows used:", x$data_info$n_used, "out of", x$data_info$n_total, "\n")
      cat("Variables from data:", paste(x$data_info$variables_in_formula, collapse = ", "), "\n")
    }
  }

  cat("\nData source: from formula and data\n")

  cat("\n--- Parametric variables (x) ---\n")
  cat("Count:", x$model_summary$n_params_x, "\n")
  cat("  ", paste(x$model_summary$x_names, collapse = ", "), "\n")

  cat("\n--- Nonparametric variables (z) ---\n")
  cat("Count:", x$model_summary$n_params_z, "\n")
  for (i in seq_along(x$model_summary$z_names)) {
    z_name <- x$model_summary$z_names[i]
    is_factor <- x$model_summary$z_is_factor[i]
    is_ordered <- x$model_summary$z_is_ordered[i]

    if (is_factor) {
      type <- if (is_ordered) "ordered factor" else "unordered factor"
    } else {
      type <- "continuous"
    }
    cat("  ", z_name, " (", type, ")\n", sep = "")
  }

  if (!is.null(x$model_summary$categorical_conversion_warning)) {
    cat("\n", x$model_summary$categorical_conversion_warning)
  }

  cat("\n--- STEP 1: SMOOTH COEFFICIENTS ---\n")
  if (isTRUE(x$smooth_coefficients$bandwidth_user_supplied)) {
    cat("Bandwidth source: user-supplied\n")
  } else {
    cat("Bandwidth source: LSCV\n")
  }
  cat("Training observations:", x$smooth_coefficients$n_training, "\n")
  cat("Evaluation observations:", x$smooth_coefficients$n_evaluation, "\n")
  if (x$smooth_coefficients$used_evaluation_data) {
    cat("Coefficients estimated on: evaluation data\n")
  } else {
    cat("Coefficients estimated on: training data\n")
  }

  cat("\n--- STEP 2: MLE RESULTS ---\n")
  if (x$model_summary$inefficiency_vars_null) {
    cat("Inefficiency variables: ALL z variables from formula's z part (no third part)\n")
  } else {
    cat("Inefficiency variables: specified in the third part of the formula\n")
  }
  cat("Number of inefficiency variable slope coefficients, delta_1:", x$model_summary$n_inefficiency_vars, "\n\n")

  print(x$mle_results$delta_table)

  cat("\n--- Variance Parameters ---\n")
  print(x$mle_results$summary_table)

  cat("\n--- TECHNICAL EFFICIENCY (BC Method) ---\n")
  cat("Mean TE:", round(x$technical_efficiency$summary["Mean"], 6), "\n")
  cat("Std Dev TE:", round(x$technical_efficiency$summary["SD"], 6), "\n")
  cat("Min TE:", round(x$technical_efficiency$summary["Min"], 6), "\n")
  cat("Max TE:", round(x$technical_efficiency$summary["Max"], 6), "\n")
  cat("Median TE:", round(x$technical_efficiency$summary["Median"], 6), "\n")

  cat("\n--- INEFFICIENCY (JLMS Method) ---\n")
  cat("Mean Inefficiency:", round(x$inefficiency$summary["Mean"], 6), "\n")
  cat("Std Dev Inefficiency:", round(x$inefficiency$summary["SD"], 6), "\n")
  cat("Min Inefficiency:", round(x$inefficiency$summary["Min"], 6), "\n")
  cat("Max Inefficiency:", round(x$inefficiency$summary["Max"], 6), "\n")
  cat("Median Inefficiency:", round(x$inefficiency$summary["Median"], 6), "\n")
  cat("\n")
  invisible(x)
}

#' Summary method for SPSC-SFA
#'
#' @param object An object of class `spsc_sfa`.
#' @param ... Additional arguments (currently ignored).
#' @return See \code{\link{print.spsc_sfa}}.
#' @export
summary.spsc_sfa <- function(object, ...) {
  print.spsc_sfa(object)
  invisible(object)
}
