#' Summarize SPSC-SFA results
#'
#' @description
#' \code{summary} method for objects of class \code{"spscsfa"}. Following the
#' convention of \code{\link[base]{summary}}, it \emph{computes and returns} a
#' summary object rather than printing one, so that the result can be assigned,
#' subsetted and reused programmatically. The summary is rendered by
#' \code{\link{print.summary.spscsfa}}, which is dispatched automatically when
#' the object is auto-printed at the console.
#'
#' @param object An object of class \code{"spscsfa"}, as returned by
#'   \code{\link{spscsfa}}.
#' @param digits Integer; the minimum number of significant digits to use when
#'   the summary is printed. Stored in the returned object and used as the
#'   default by \code{\link{print.summary.spscsfa}}, where it can be
#'   overridden.
#' @param ... Further arguments passed to or from other methods. Accepted for
#'   compatibility with the generic \code{\link[base]{summary}} and otherwise
#'   unused.
#'
#' @return An object of class \code{"summary.spscsfa"}: a list with components
#'   \describe{
#'     \item{formula, n_obs, data_supplied, data_info}{Model formula, number of
#'       observations used, whether a data frame was supplied, and the
#'       data cleaning information.}
#'     \item{n_params_x, x_names}{Number and names of the parametric variables.}
#'     \item{n_params_z, z_names, z_is_factor, z_is_ordered}{Number, names and
#'       types of the nonparametric variables.}
#'     \item{categorical_conversion_warning}{Note about categorical inefficiency
#'       determinants converted to binary variables, or \code{NULL}.}
#'     \item{spscreg}{The estimator of the smooth coefficients:
#'       \code{"lc"} (local-constant) or \code{"ll"} (local-linear).}
#'     \item{bandwidth, bandwidth_user_supplied, n_training, n_evaluation,
#'       used_evaluation_data}{Step 1 smooth-coefficient information.}
#'     \item{smooth_coefficients}{Numeric matrix with one row for each smooth
#'       coefficient (\eqn{\theta(z)}, the intercept, and then each \eqn{\beta(z)})
#'       and the columns \code{Min.}, \code{1st Qu.}, \code{Median},
#'       \code{Mean}, \code{3rd Qu.} and \code{Max.}: the summary statistics
#'       of the Step 1 estimates over the evaluation points (the rows of
#'       \code{ezdat} if it was supplied to \code{\link{spscsfa}}, the training
#'       observations otherwise), or \code{NULL} if the estimates are not
#'       available.}
#'     \item{inefficiency_vars_null, n_inefficiency_vars, sigma_u_formula}{
#'       Specification of the inefficiency equation.}
#'     \item{alpha}{Numeric matrix with one row, \code{Intercept (alpha)}, and
#'       the columns \code{Min.}, \code{1st Qu.}, \code{Median}, \code{Mean},
#'       \code{3rd Qu.} and \code{Max.}: the summary statistics of the
#'       intercept \eqn{\alpha(z) = \theta(z) + E(u|z)} over the training
#'       observations, or \code{NULL} if it is not available.}
#'     \item{coefficients}{Numeric matrix of the \code{delta} coefficients with
#'       columns \code{Estimate}, \code{Std. Error}, \code{z value} and
#'       \code{Pr(>|z|)}.}
#'     \item{variance}{Numeric matrix, in the same layout, for
#'       \code{sigma2_v}.}
#'     \item{logLik, optim_method}{Maximized log-likelihood and optimizer used.}
#'     \item{te, inefficiency}{Named numeric vectors of summary statistics
#'       (\code{Mean}, \code{SD}, \code{Min}, \code{Max}, \code{Median}) for the
#'       BC technical efficiency and JLMS inefficiency estimates.}
#'     \item{sig_level}{Significance level used for the technical efficiency
#'       confidence intervals.}
#'     \item{gradients}{\code{NULL} unless the model was estimated with
#'       \code{gradients = TRUE}. Otherwise a list with the elements
#'       \code{theta_beta_grad} (a list with one matrix for each continuous
#'       \code{z} variable), \code{alpha_grad} and \code{inefficiency_grad}
#'       (matrices), each with one row for each gradient and the columns
#'       \code{Min.}, \code{1st Qu.}, \code{Median}, \code{Mean},
#'       \code{3rd Qu.} and \code{Max.}: the summary statistics of that
#'       gradient. For \code{theta_beta_grad} they
#'       are computed over the evaluation points: the rows of \code{ezdat} if it
#'       was supplied to \code{\link{spscsfa}}, the training observations
#'       otherwise. For \code{alpha_grad} and \code{inefficiency_grad} they are
#'       always computed over the training observations.}
#'     \item{digits}{Default number of significant digits for printing.}
#'   }
#'
#' @seealso \code{\link{print.summary.spscsfa}}, \code{\link{print.spscsfa}},
#'   \code{\link{spscsfa}}
#' @export
summary.spscsfa <- function(object, digits = max(3L, getOption("digits") - 3L),
                             ...) {
  if (!.is_spscsfa(object)) {
    stop("'object' must be of class 'spscsfa'")
  }

  ms <- object$model_summary
  sc <- object$smooth_coefficients
  ml <- object$mle_results

  n_z <- ml$n_z
  if (is.null(n_z)) n_z <- ms$n_inefficiency_vars

  ans <- list(
    formula        = ms$formula,
    n_obs          = ms$n_obs,
    data_supplied  = isTRUE(ms$data_supplied),
    data_info      = object$data_info,
    n_params_x     = ms$n_params_x,
    x_names        = ms$x_names,
    n_params_z     = ms$n_params_z,
    z_names        = ms$z_names,
    z_is_factor    = ms$z_is_factor,
    z_is_ordered   = ms$z_is_ordered,
    categorical_conversion_warning = ms$categorical_conversion_warning,

    ## --- Step 1: smooth coefficients ------------------------------------
    spscreg                 = if (is.null(sc$spscreg)) "lc" else sc$spscreg,
    bandwidth               = sc$bandwidth,
    bandwidth_user_supplied = isTRUE(sc$bandwidth_user_supplied),
    n_training              = sc$n_training,
    n_evaluation            = sc$n_evaluation,
    used_evaluation_data    = isTRUE(sc$used_evaluation_data),
    smooth_coefficients     = .spsc_col_stats(sc$theta_beta),

    ## --- Step 2: maximum likelihood --------------------------------------
    inefficiency_vars_null = isTRUE(ms$inefficiency_vars_null),
    n_inefficiency_vars    = ms$n_inefficiency_vars,
    sigma_u_formula        = if (isTRUE(n_z == 0)) {
      "exp(delta_0)"
    } else {
      "exp(delta_0 + delta_1'z)"
    },
    alpha        = .spsc_col_stats(sc$alpha),
    coefficients = .spsc_coefmat(object),
    variance     = .spsc_varmat(object),
    logLik       = ml$log_likelihood,
    optim_method = ml$optim_method,

    ## --- efficiency -------------------------------------------------------
    te           = object$technical_efficiency$summary,
    inefficiency = object$inefficiency$summary,
    sig_level    = object$technical_efficiency$sig_level,

    ## --- gradients (only if estimated with gradients = TRUE) ---------------
    gradients    = .spsc_grad_summary(object),

    ## --- printing default --------------------------------------------------
    digits = digits
  )

  class(ans) <- "summary.spscsfa"
  ans
}

#' @method summary spsc_sfa
#' @export
#' @noRd
summary.spsc_sfa <- summary.spscsfa

#' Print an SPSC-SFA summary
#'
#' @description
#' \code{print} method for objects of class \code{"summary.spscsfa"} created by
#' \code{\link{summary.spscsfa}}. Numeric output is controlled by the standard
#' arguments of \code{\link[base]{print}}: \code{digits} and
#' \code{signif.stars} are handled explicitly, and any further arguments in
#' \code{...} (for example \code{na.print}, \code{print.gap}, \code{quote} or
#' \code{right}) are passed down to \code{\link[stats]{printCoefmat}} and
#' \code{\link[base]{print.default}}.
#'
#' @param x An object of class \code{"summary.spscsfa"}.
#' @param digits Integer; the minimum number of significant digits. Defaults to
#'   the value stored by \code{\link{summary.spscsfa}}.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#'   next to the p-values of the inefficiency determinants. Defaults to
#'   \code{getOption("show.signif.stars")}.
#' @param signif.legend Logical; if \code{TRUE}, the legend for the significance
#'   stars is printed. Defaults to \code{signif.stars}.
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}} and
#'   \code{\link[base]{print.default}}.
#'
#' @return The argument \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.spscsfa}}, \code{\link{print.spscsfa}}
#' @method print summary.spscsfa
#' @export
print.summary.spscsfa <- function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   signif.legend = signif.stars,
                                   ...) {
  ## Use the digits setting from summary() unless the caller provides another digits value.
  if (missing(digits) && !is.null(x$digits)) digits <- x$digits

  ## Points over which the Step 1 estimates (and their gradients) are
  ## summarized: the evaluation points if 'ezdat' was used, else the training data
  eval_pts <- if (isTRUE(x$used_evaluation_data)) {
    "the evaluation points (the rows of ezdat)"
  } else {
    "the training observations"
  }

  cat("\n")
  cat("==================== SPSC-SFA RESULTS ====================\n")
  cat("\nFormula:", .spsc_deparse(x$formula), "\n")
  cat("Number of observations:", x$n_obs, "\n")

  if (isTRUE(x$data_supplied)) {
    cat("\n--- Data Information ---\n")
    cat("Data frame supplied: YES\n")
    if (!is.null(x$data_info)) {
      cat("Rows used:", x$data_info$n_used, "out of", x$data_info$n_total, "\n")
      cat("Variables from data:",
          paste(x$data_info$variables_in_formula, collapse = ", "), "\n")
    }
  }

  cat("\nData source: from formula and data\n")

  cat("\n--- Parametric variables (x) ---\n")
  cat("Count:", x$n_params_x, "\n")
  cat("  ", paste(x$x_names, collapse = ", "), "\n")

  cat("\n--- Nonparametric variables (z) ---\n")
  cat("Count:", x$n_params_z, "\n")
  for (i in seq_along(x$z_names)) {
    type <- if (isTRUE(x$z_is_factor[i])) {
      if (isTRUE(x$z_is_ordered[i])) "ordered factor" else "unordered factor"
    } else {
      "continuous"
    }
    cat("  ", x$z_names[i], " (", type, ")\n", sep = "")
  }

  if (!is.null(x$categorical_conversion_warning)) {
    cat("\n", x$categorical_conversion_warning, sep = "")
  }

  cat("\n--- STEP 1: SMOOTH COEFFICIENTS ---\n")
  cat("Smooth-coefficient estimator:",
      if (identical(x$spscreg, "ll")) "local-linear (LL)" else "local-constant (LC)",
      "\n")
  cat("Bandwidth source:",
      if (isTRUE(x$bandwidth_user_supplied)) "user-supplied" else "LSCV", "\n")
  cat("Training observations:", x$n_training, "\n")
  cat("Evaluation observations:", x$n_evaluation, "\n")
  cat("Coefficients estimated on:",
      if (isTRUE(x$used_evaluation_data)) "evaluation data" else "training data",
      "\n")

  if (!is.null(x$smooth_coefficients)) {
    cat("\nSmooth Coefficients Summary (over ", eval_pts, "):\n", sep = "")
    print(x$smooth_coefficients, digits = digits, ...)
  }

  cat("\n--- STEP 2: MLE RESULTS ---\n")
  if (isTRUE(x$inefficiency_vars_null)) {
    cat("Inefficiency variables: ALL z variables from formula's z part",
        "(no third part)\n")
  } else {
    cat("Inefficiency variables: specified in the third part of the formula\n")
  }
  cat("Number of inefficiency variable slope coefficients, delta_1:",
      x$n_inefficiency_vars, "\n")

  if (!is.null(x$alpha)) {
    cat("\nIntercept (alpha) Summary (over the training observations):\n")
    print(x$alpha, digits = digits, ...)
  }

  if (!is.null(x$coefficients)) {
    cat("\nInefficiency determinants (sigma_u = ", x$sigma_u_formula, "):\n\n",
        sep = "")
    stats::printCoefmat(x$coefficients, digits = digits,
                        signif.stars = signif.stars,
                        signif.legend = signif.legend, ...)
  }

  if (!is.null(x$variance)) {
    cat("\n--- Variance Parameters ---\n")
    stats::printCoefmat(x$variance, digits = digits,
                        signif.stars = FALSE, signif.legend = FALSE, ...)
  }

  if (!is.null(x$logLik)) {
    cat("\nLog-likelihood:", format(x$logLik, digits = digits))
    if (!is.null(x$optim_method)) cat("   (optimizer:", x$optim_method, ")")
    cat("\n")
  }

  cat("\n--- TECHNICAL EFFICIENCY (BC Method) ---\n")
  if (!is.null(x$sig_level)) {
    cat("Confidence level for TE intervals: ",
        format(100 * (1 - x$sig_level), digits = digits), "%\n", sep = "")
  }
  if (!is.null(x$te)) print(x$te, digits = digits, ...)

  cat("\n--- INEFFICIENCY (JLMS Method) ---\n")
  if (!is.null(x$inefficiency)) print(x$inefficiency, digits = digits, ...)

  ## Gradients: reported whenever they were estimated (gradients = TRUE)
  if (!is.null(x$gradients)) {
    g <- x$gradients
    est <- if (identical(x$spscreg, "ll")) "local-linear" else "local-constant"

    cat("\n--- GRADIENTS ---\n")

    if (length(g$theta_beta_grad) > 0L) {
      cat("\nSmooth coefficients, theta (intercept) and beta (", est,
          " estimates): gradients with respect to each continuous z variable\n",
          "Summary statistics over ", eval_pts, ".\n", sep = "")
      for (zn in names(g$theta_beta_grad)) {
        cat("\nWith respect to ", zn, ":\n", sep = "")
        print(g$theta_beta_grad[[zn]], digits = digits, ...)
      }
    } else {
      cat("\nNo continuous z variable: no gradients of the smooth coefficients.\n")
    }

    ## alpha_grad and inefficiency_grad are both summarized over the training
    ## observations
    train_note <- "Summary statistics over the training observations.\n"

    if (!is.null(g$alpha_grad)) {
      cat("\n", train_note, sep = "")
      cat("Intercept alpha: gradients with respect to each continuous z variable\n")
      print(g$alpha_grad, digits = digits, ...)
    }

    if (!is.null(g$inefficiency_grad)) {
      if (is.null(g$alpha_grad)) cat("\n", train_note, sep = "")
      cat("\nJLMS inefficiency: gradients with respect to the inefficiency",
          "determinants\n")
      print(g$inefficiency_grad, digits = digits, ...)
    }
  }

  cat("\n")
  invisible(x)
}

#' @method print summary.spsc_sfa
#' @export
#' @noRd
print.summary.spsc_sfa <- print.summary.spscsfa

#' Print SPSC-SFA results
#'
#' @description
#' \code{print} method for objects of class \code{"spscsfa"}. It formats the
#' fit through \code{\link{summary.spscsfa}}, so the display is identical to
#' that of \code{summary(x)} and responds to the same standard arguments of
#' \code{\link[base]{print}}. As required for \code{print} methods, the object
#' itself is returned invisibly.
#'
#' @param x An object of class \code{"spscsfa"}, as returned by
#'   \code{\link{spscsfa}}.
#' @param digits Integer; the minimum number of significant digits to use.
#' @param ... Further arguments passed to
#'   \code{\link{print.summary.spscsfa}}, and from there to
#'   \code{\link[stats]{printCoefmat}} and \code{\link[base]{print.default}}
#'   (for example \code{signif.stars}, \code{na.print} or \code{right}).
#'
#' @return The original object of class \code{"spscsfa"}, returned invisibly.
#'   It is a list containing the components \code{smooth_coefficients},
#'   \code{mle_results}, \code{technical_efficiency}, \code{inefficiency},
#'   \code{model_summary}, \code{formula_parsed} and \code{data_info}. These
#'   components contain the estimated smooth coefficients and recovered
#'   intercept, maximum likelihood results, technical efficiency estimates and
#'   confidence intervals, JLMS inefficiency estimates, model specification
#'   information, parsed formula components, and data cleaning information. If
#'   the model was estimated with \code{gradients = TRUE},
#'   \code{smooth_coefficients_grad} and \code{inefficiency_grad} are also returned.
#'
#' @seealso \code{\link{summary.spscsfa}},
#'   \code{\link{print.summary.spscsfa}}, \code{\link{spscsfa}}
#' @export
print.spscsfa <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x, digits = digits), digits = digits, ...)
  invisible(x)
}

#' @method print spsc_sfa
#' @export
#' @noRd
print.spsc_sfa <- print.spscsfa

## ---------------------------------------------------------------------------
## Internal helpers (not exported, not documented)
## ---------------------------------------------------------------------------

## Collapse a language object (e.g. a formula) into a single-line string.
.spsc_deparse <- function(x) {
  if (is.null(x)) return("NULL")
  paste(deparse(x), collapse = " ")
}

## Standard four-column coefficient matrix for the inefficiency (delta)
## equation: Estimate / Std. Error / z value / Pr(>|z|).  The layout and the
## column names are the ones expected by stats::printCoefmat().
.spsc_coefmat <- function(object) {
  tab <- object$mle_results$delta_table
  cn <- c("Estimate", "Std.Error", "t.value", "p.value")
  
  if (!is.null(tab) && all(cn %in% names(tab))) {
    cm <- as.matrix(tab[, cn, drop = FALSE])
    rn <- if (!is.null(tab$Coefficient)) as.character(tab$Coefficient) else rownames(tab)
  } else {
    est <- object$mle_results$delta_coef
    if (is.null(est)) return(NULL)
    se <- object$mle_results$delta_se
    if (is.null(se)) se <- NA_real_
    zval <- est / se
    cm <- cbind(est, se, zval, 2 * (1 - stats::pnorm(abs(zval))))
    rn <- names(est)
    if (is.null(rn)) {
      rn <- c("Intercept", paste0("delta_", seq_len(max(0L, length(est) - 1L))))
    }
  }
  
  storage.mode(cm) <- "double"
  dimnames(cm) <- list(rn, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  cm
}

## One-row matrix with the variance parameter sigma2_v, same layout as above.
.spsc_varmat <- function(object) {
  est <- object$mle_results$sigma2_v
  if (is.null(est)) return(NULL)
  se <- object$mle_results$sigma2_v_se
  if (is.null(se)) se <- NA_real_
  zval <- est / se
  vm <- matrix(as.double(c(est, se, zval, 2 * (1 - stats::pnorm(abs(zval))))),
               nrow = 1L)
  dimnames(vm) <- list("sigma2_v",
                       c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  vm
}

## Summary statistics (the ones of summary() of a numeric vector: Min., 1st Qu.,
## Median, Mean, 3rd Qu., Max.) of each column of a matrix: one row for each
## column.  Used for the smooth coefficients, theta (the intercept) and beta,
## estimated in Step 1, and for each gradient.
.spsc_col_stats <- function(m) {
  if (is.null(m)) return(NULL)
  m <- as.matrix(m)
  if (ncol(m) == 0L || nrow(m) == 0L) return(NULL)
  
  st <- t(apply(m, 2L, function(v) {
    q <- stats::quantile(v, probs = c(0, 0.25, 0.5, 0.75, 1),
                         na.rm = TRUE, names = FALSE)
    c(Min. = q[1L], `1st Qu.` = q[2L], Median = q[3L],
      Mean = mean(v, na.rm = TRUE), `3rd Qu.` = q[4L], Max. = q[5L])
  }))
  rownames(st) <- colnames(m)
  st
}

## Summary statistics of all the gradients of a fitted object, or NULL if the
## model was estimated without gradients.
.spsc_grad_summary <- function(object) {
  if (!isTRUE(object$model_summary$gradients)) return(NULL)
  
  scg <- object$smooth_coefficients_grad
  
  list(
    theta_beta_grad = if (length(scg$theta_beta_grad) > 0L) {
      lapply(scg$theta_beta_grad, .spsc_col_stats)
    } else {
      NULL
    },
    alpha_grad = if (!is.null(scg$alpha_grad)) .spsc_col_stats(scg$alpha_grad),
    inefficiency_grad = if (!is.null(object$inefficiency_grad)) {
      .spsc_col_stats(object$inefficiency_grad)
    }
  )
}
