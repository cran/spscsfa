.is_spscsfa <- function(x) inherits(x, c("spscsfa", "spsc_sfa"))

#' Extract technical efficiency scores
#'
#' @param x An object of class \code{"spscsfa"}.
#' @return A numeric vector of technical efficiency estimates.
#' @export
te <- function(x) {
  if (!.is_spscsfa(x)) stop("x must be of class 'spscsfa'")
  x$technical_efficiency$TE
}

#' Extract JLMS inefficiency estimates
#'
#' @param x An object of class `spscsfa`.
#' @return A numeric vector of JLMS conditional inefficiency estimates.
#' @export
inefficiency <- function(x) {
  if (!.is_spscsfa(x)) stop("x must be of class 'spscsfa'")
  x$inefficiency$u_JLMS
}

#' Extract smooth coefficients
#'
#' @description
#' Extracts the smooth coefficients of a fitted model. \code{smooth_coef} is an alias of
#' \code{smoothcoef}.
#'
#' @param x An object of class `spscsfa`.
#' @return A list with two components:
#'   \item{theta_beta}{The matrix of the smooth coefficients of the first step,
#'     \eqn{\theta(z)} and \eqn{\beta(z)}: one row for each evaluation point
#'     (those of \code{ezdat} if it was supplied to \code{\link{spscsfa}}, the
#'     training data otherwise), and the columns \code{Intercept (theta)} and the
#'     names of the \code{x} variables. They are local-constant or local-linear
#'     estimates, as chosen by \code{spscreg}.}
#'   \item{alpha}{A matrix with one column: the estimates of the intercept
#'     \eqn{\alpha(z) = \theta(z) + E(u|z)} at the training data.}
#' @export
smoothcoef <- function(x) {
  if (!.is_spscsfa(x)) stop("x must be of class 'spscsfa'")
  list(
    theta_beta = x$smooth_coefficients$theta_beta,
    alpha = x$smooth_coefficients$alpha
  )
}

#' @rdname smoothcoef
#' @usage NULL
#' @export
smooth_coef <- smoothcoef

#' Extract gradient estimates
#'
#' @description
#' Extracts the gradients estimated by \code{\link{spscsfa}} with
#' \code{gradients = TRUE}: the gradients of the smooth coefficients
#' \eqn{\theta} and \eqn{\beta} and of the intercept \eqn{\alpha} with
#' respect to the continuous nonparametric variables \code{z}, and the gradients
#' of the JLMS inefficiency estimates with respect to the inefficiency
#' determinants. The smooth coefficients are local-constant if the model was
#' estimated with \code{spscreg = "lc"} and local-linear if it was estimated
#' with \code{spscreg = "ll"}.
#'
#' @param x An object of class \code{"spscsfa"} estimated with
#'   \code{gradients = TRUE}.
#' @param zvar Optional character vector with the names of the variables whose
#'   gradients are extracted, as they are written in \code{formula} (they are
#'   the names of \code{x$model_summary$z_names} for the \code{z} variables
#'   and \code{x$mle_results$inefficiency_columns} for the inefficiency
#'   determinants), or their positions among the \code{z} variables. If
#'   \code{NULL} (default), all the gradients are returned. For example,
#'   \code{grad(x, zvar = "z1")} extracts all the gradients with respect to
#'   \code{z1}.
#'
#'   If \code{z1} is not a continuous \code{z} variable (it is a factor), there
#'   is no gradient of the smooth coefficients with respect to it: only
#'   \code{inefficiency_grad}, for the set of binary variables converted from
#'   \code{z1} in the inefficiency equation, can be extracted. Note that 
#'   \code{inefficiency_grad} is not available for a continuous \code{z} variable 
#'   that is not an inefficiency determinant.
#'
#' @return A list with three elements:
#'   \item{theta_beta_grad}{The gradients of the smooth coefficients: a list
#'     with one matrix for each continuous \code{z} variable (named after it).
#'     The matrix for \code{z_l} has the columns
#'     \code{d(Intercept (theta))/d(z_l)}, \code{d(x1)/d(z_l)}, \ldots, \code{d(xp)/d(z_l)}, 
#'     that is, the gradients of \eqn{\theta, \beta_1, \ldots, \beta_p} with respect
#'     to \code{z_l}, with one row for each evaluation point (\code{ezdat} if
#'     it was supplied, the training data otherwise).}
#'   \item{alpha_grad}{The gradients of \eqn{\alpha} at the training data: a
#'     matrix with one column for each continuous \code{z} variable, named
#'     \code{d(Intercept (alpha))/d(z_l)}.}
#'   \item{inefficiency_grad}{The gradients of the JLMS inefficiency estimates
#'     with respect to the inefficiency determinants: a matrix with one column
#'     for each determinant (the binary variables for a categorical one), or
#'     \code{NULL} if the model has no inefficiency determinant.}
#'   With \code{zvar}, the elements are restricted to the variables chosen. An
#'   element that is not available for them is \code{NULL}.
#'
#' @seealso \code{\link{smoothcoef}}, \code{\link{spscsfa}}
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' z1 <- runif(n)
#' z2 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' y <- cos(z1) + sin(z1) * x + rnorm(n, 0, 0.1) -
#'      exp(0.1 + 0.1 * z1) * abs(rnorm(n, 0, 1))
#' dat <- data.frame(y, x, z1, z2)
#' result <- spscsfa(y ~ x | z1 + z2, data = dat, gradients = TRUE,
#'                    verbose = FALSE)
#'
#' # all the gradients
#' g <- grad(result)
#' names(g)
#'
#' # the gradients with respect to z1
#' g1 <- grad(result, zvar = "z1")
#' head(g1$theta_beta_grad$z1)
#'
#' # z2 is not continuous: only inefficiency_grad can be extracted, and a note
#' # is printed
#' g2 <- grad(result, zvar = "z2")
#' head(g2$inefficiency_grad)
#' 
#' @export
grad <- function(x, zvar = NULL) {

  if (!.is_spscsfa(x)) stop("x must be of class 'spscsfa'")

  scg <- x$smooth_coefficients_grad

  if (!isTRUE(x$model_summary$gradients) ||
      (is.null(scg) && is.null(x$inefficiency_grad))) {
    stop("'x' has no gradients: estimate the model with ",
         "spscsfa(..., gradients = TRUE)", call. = FALSE)
  }

  out <- list(
    theta_beta_grad = scg$theta_beta_grad,
    alpha_grad = scg$alpha_grad,
    inefficiency_grad = x$inefficiency_grad
  )

  if (is.null(zvar)) {
    return(out)
  }

  # ---- variables with gradients -------------------------------------------

  z_names <- x$model_summary$z_names

  # continuous z variables: those of theta_beta_grad (the columns of alpha_grad
  # are in the same order)
  z_cont <- names(scg$theta_beta_grad)

  # inefficiency determinants and the columns of inefficiency_grad
  ie_cols <- x$mle_results$inefficiency_columns
  if (is.null(ie_cols)) {
    ie_cols <- stats::setNames(as.list(colnames(x$inefficiency_grad)),
                               colnames(x$inefficiency_grad))
  }
  ie_names <- names(ie_cols)

  if (is.numeric(zvar)) {
    if (anyNA(zvar) || any(zvar < 1 | zvar > length(z_names))) {
      stop("'zvar' must select among the z variables: ",
           paste(z_names, collapse = ", "), call. = FALSE)
    }
    zvar <- z_names[zvar]
  }

  if (!is.character(zvar) || length(zvar) == 0L || anyNA(zvar)) {
    stop("'zvar' must be the name of a variable", call. = FALSE)
  }

  zvar <- unique(zvar)
  unknown <- setdiff(zvar, union(z_names, ie_names))

  if (length(unknown) > 0L) {
    stop("unknown variable(s) in 'zvar': ", paste(unknown, collapse = ", "),
         ". Available: ", paste(union(z_names, ie_names), collapse = ", "),
         call. = FALSE)
  }

  # ---- restrict the gradients to zvar ---------------------------------------

  zc <- zvar[zvar %in% z_cont]
  ie <- zvar[zvar %in% ie_names]

  res <- list(theta_beta_grad = NULL, alpha_grad = NULL,
              inefficiency_grad = NULL)

  if (length(zc) > 0L) {
    res$theta_beta_grad <- scg$theta_beta_grad[zc]
    res$alpha_grad <- scg$alpha_grad[, match(zc, z_cont), drop = FALSE]
  }

  if (length(ie) > 0L && !is.null(x$inefficiency_grad)) {
    cols <- unlist(ie_cols[ie], use.names = FALSE)
    res$inefficiency_grad <- x$inefficiency_grad[, cols, drop = FALSE]
  }

  # ---- notes on what is not available ---------------------------------------

  for (v in zvar) {
    is_z <- v %in% z_names
    is_cont <- v %in% z_cont
    is_ie <- v %in% ie_names

    if (is_cont && is_ie) next

    if (is_z && !is_cont && is_ie) {
      msg <- paste0("'", v, "' is not a continuous z variable, so there are no ",
                    "gradients of the smooth coefficients (theta_beta_grad) or ",
                    "of alpha (alpha_grad) with respect to it. Only ",
                    "inefficiency_grad, for the binary variable(s) converted ",
                    "from it (", paste(ie_cols[[v]], collapse = ", "),
                    "), can be extracted")
    } else if (is_z && !is_cont) {
      msg <- paste0("'", v, "' is not a continuous z variable and not an ",
                    "inefficiency determinant: there is no gradient with ",
                    "respect to it")
    } else if (is_cont) {
      msg <- paste0("'", v, "' is not an inefficiency determinant, so ",
                    "inefficiency_grad is not available for it")
    } else {
      msg <- paste0("'", v, "' is not a z variable, so there are no gradients ",
                    "of the smooth coefficients (theta_beta_grad) or of alpha ",
                    "(alpha_grad) with respect to it. Only inefficiency_grad ",
                    "can be extracted")
    }

    message("Note: ", msg, ".")
  }

  res
}
