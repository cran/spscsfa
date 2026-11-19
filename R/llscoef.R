#' Local-linear semiparametric smooth coefficient estimation
#'
#' @description
#' Estimates the semiparametric smooth coefficient (SPSC) model by local-linear
#' kernel regression, at the training data or at new evaluation points. The bandwidth 
#' object of \code{\link{llscoefbw}} can be passed as \code{bws}. It is used 
#' by \code{\link{spscsfa}} for the estimation of the smooth coefficients.
#'
#' @details
#' The model is \eqn{y_i = x_i'\beta(z_i) + u_i}, where \eqn{x_i} consists of a
#' constant followed by the variables in the first part of \code{formula} and
#' \eqn{z_i} holds the variables in the second part. A third part, as accepted
#' by \code{\link{spscsfa}} for the inefficiency determinants, is ignored.
#' Observations with missing values are removed.
#'
#' For a given \eqn{z}, the local-linear estimator following Geng and Sun (2019) is
#' \deqn{\hat\gamma(z) = \Big(\sum_i S_i S_i' K_i(z)\Big)^{-1}
#'       \sum_i S_i y_i K_i(z),}
#' where \eqn{S_i = (x_i', (x_i \otimes (z_i - z))')'} involves the continuous
#' variables in \eqn{z} only and \eqn{K_i(z)} is a generalized product kernel
#' (Second-order Gaussian for continuous, Aitchison-Aitken for unordered and Li-Racine for
#' ordered variables, the defaults of \code{np::npksum}). The first elements of
#' \eqn{\hat\gamma(z)} estimate the intercept and slope smooth coefficients and the
#' remaining ones their gradients. Variables of the second part of
#' \code{formula} are classified by their class: numeric variables are
#' continuous, factors are unordered categorical and ordered factors are
#' ordered categorical.
#'
#' \describe{
#'   \item{Formula and bandwidths}{\code{bws} is the first argument, and 
#'     what it is decides which of the other arguments are
#'     needed. The formula is given once, to \code{llscoefbw}:
#'     \itemize{
#'       \item \code{bws} is the object returned by \code{llscoefbw} with
#'         \code{LSCV = TRUE}: \code{formula} must not be given,
#'         and \code{data} is not needed. The object keeps the formula and the
#'         call to \code{llscoefbw}, whose \code{data} are evaluated again in
#'         the environment the call was made from. \code{data} replaces the
#'         training data, \code{newdata} gives the evaluation points.
#'       \item \code{bws} is a formula: the bandwidths are selected first,
#'         by \code{llscoefbw(bws, data = data)}, and used in the same call.
#'       \item \code{bws} is a numeric vector of local-linear bandwidths, such
#'         as the result of \code{llscoefbw} with
#'         \code{LSCV = FALSE}: \code{formula} is required, and
#'         \code{data} unless the variables are in the environment of
#'         \code{formula}. A bandwidth object of \code{np::npscoefbw} is accepted
#'         too, and its bandwidths are used, if it was made with
#'         \code{regtype = "ll"}, fixed bandwidths (\code{bwtype = "fixed"},
#'         \code{bwscaling = FALSE}) and the default kernels (\code{ckertype},
#'         \code{ckerorder}, \code{ckerbound}, \code{ukertype} and
#'         \code{okertype} arguments). However, one made with \code{regtype = "lc"}, the default 
#'         \code{regtype} argument of \code{np::npscoefbw}, has bandwidths for the 
#'         local-constant estimator, and gives an error like the other objects that 
#'         do not qualify. A user-supplied vector that does not qualify cannot be told 
#'         apart from other vectors, but such a vector is not suitable, either.
#'     }
#'     A bandwidth object returned by \code{llscoefbw} cannot be used together
#'     with a \code{formula}. This is because an object of class \code{"llsmoothcoefficientbw"} 
#'     stores the formula and data used to select the bandwidths. Therefore, 
#'     formula is not allowed to appear again to guarantee internal consistency.}
#'   \item{Evaluation points}{The estimation always uses
#'     the training data: those of the call to \code{llscoefbw} for a bandwidth
#'     object, otherwise those given by \code{formula} and \code{data}.
#'     \code{newdata} only decides where the smooth coefficients are evaluated. It
#'     needs the variables of the \code{x} and \code{z} parts of \code{formula}
#'     (the response is not needed and is ignored). The \code{x} variables
#'     matter for the fitted values only, and any values, such as the training
#'     means, can be used when only the smooth coefficients are of interest.
#'     Factors in \code{newdata} are matched to the levels of the training
#'     data, and rows of \code{newdata} with missing values are dropped.}
#'   \item{Residuals}{The residuals are those of the training data, that is, 
#'   \eqn{y} minus the fitted values at the training points, also when \code{newdata} is supplied.}
#'   \item{Gradients}{The gradients of the smooth coefficients with respect to
#'     the continuous variables in \code{z}, one column for each coefficient
#'     and variable, in X-major order: all the gradients of the first
#'     coefficient (the intercept) come first, then those of the second one,
#'     and so on, each in the order of the continuous \code{z} variables.}
#' }
#'
#' @param bws A bandwidth specification: an object of class
#'   \code{"llsmoothcoefficientbw"} returned by \code{\link{llscoefbw}}, a
#'   formula \code{y ~ x | z} (the bandwidths are then selected by
#'   \code{llscoefbw}), or a numeric vector of local-linear bandwidths, one for
#'   each variable in the \code{z} part of \code{formula}. If the bandwidths
#'   have names, as in the result of \code{llscoefbw}, they are matched to the
#'   names of the \code{z} variables, in any order: a warning is given if their
#'   order differs from that of the \code{z} variables in the formula, and an
#'   error if the names do not match. Bandwidths without names are taken in the
#'   order of the \code{z} variables. A bandwidth object of \code{np::npscoefbw}
#'   made with \code{regtype = "ll"}, fixed bandwidths and the default kernels
#'   is also accepted, and its bandwidths are used; one made with another
#'   \code{regtype}, \code{bwtype}, \code{bwscaling} or kernel is not, and a
#'   vector taken from it is not suitable, either. See \sQuote{Details} for the
#'   arguments each kind of \code{bws} needs.
#' @param formula A formula \code{y ~ x | z} (or \code{y ~ x | z | q}), where
#'   \code{y} is the dependent variable, \code{x} are the parametric regressors
#'   (the constant is added automatically) and \code{z} are the variables that
#'   the smooth coefficients depend on. Required if \code{bws} is a numeric
#'   vector. Must not be given if \code{bws} is a bandwidth object from
#'   \code{llscoefbw}, which already keeps the formula.
#' @param data Optional data frame containing the variables in \code{formula}.
#'   If omitted, the data of the call to \code{llscoefbw} are used for a
#'   bandwidth object, and otherwise the variables are taken from the
#'   environment of \code{formula}.
#' @param newdata Optional data frame with the evaluation points. If omitted,
#'   the training data are used. See \sQuote{Details}.
#' @param betas Logical; if \code{TRUE}, the smooth coefficients are computed
#'   and returned. Defaults to \code{FALSE}.
#' @param residuals Logical; if \code{TRUE}, the residuals are computed and
#'   returned. Defaults to \code{FALSE}.
#' @param gradients Logical; if \code{TRUE}, the gradients of the smooth
#'   coefficients are computed and returned. Defaults to \code{FALSE}.
#' @param na.action As in \code{\link{llscoefbw}}: what to do with rows that
#'   have missing values. It is used by the methods for a formula and for a
#'   vector of bandwidths. A bandwidth object from \code{llscoefbw} carries its
#'   own, which is used, so \code{na.action} is not an argument then. It applies to
#'   \code{newdata} too.
#' @param object,x An object of class \code{"llsmoothcoefficient"}.
#' @param ... Further arguments, not used.
#'
#' @return \code{llscoef} returns an object of class
#'   \code{"llsmoothcoefficient"}, a list with components
#'   \item{mean}{The fitted values at the evaluation points.}
#'   \item{beta}{If \code{betas = TRUE}, the matrix of smooth coefficients, with
#'     one row for each evaluation point and the columns \code{"Intercept"} and
#'     the \code{x} variables.}
#'   \item{resid}{If \code{residuals = TRUE}, the residuals, \eqn{y} minus the
#'     fitted values at the training points.}
#'   \item{grad}{If \code{gradients = TRUE}, the matrix of gradients of the
#'     smooth coefficients with respect to the continuous variables in
#'     \code{z}, in X-major order, with one row for each evaluation point.}
#'   Together with \code{bw}, \code{x_names}, \code{z_names}, \code{n_train},
#'   \code{n_eval}, \code{eval_train} (evaluation and training data are the same),
#'   \code{omit} (the rows with missing values that were dropped from the
#'   evaluation data: from \code{newdata} if it is given, and otherwise from the training data) and
#'   \code{call}. With \code{na.action = na.exclude}, \code{mean},
#'   \code{beta} and \code{grad} are padded with \code{NA} at the rows of the
#'   evaluation data that were dropped (those in \code{omit}), and
#'   \code{resid} at the dropped rows of the training data, so that they line
#'   up with \code{newdata} or with the data. \code{n_train} and \code{n_eval}
#'   still count the rows that were used.
#'   
#' @note For an object of class \code{"llsmoothcoefficient"}, the smooth coefficients and residuals can be extracted with
#'   \code{\link{coef}} and \code{\link{residuals}}, respectively, and the fitted values with
#'   \code{\link{fitted}}.
#'
#' @references
#' Xin Geng and Kai Sun (2019). Gradient estimation of the local-constant
#' semiparametric smooth coefficient model. \emph{Economics Letters}, 185,
#' 108684.
#'
#' Tristen Hayfield and Jeffrey S. Racine (2008). Nonparametric econometrics: The np
#' package. \emph{Journal of Statistical Software}, 27(5), 1--32.
#'
#' Qi Li, Cliff J. Huang, Dong Li and Tsu-Tan Fu (2002). Semiparametric smooth
#' coefficient models. \emph{Journal of Business & Economic Statistics}, 20(3),
#' 412--422.
#'
#' @seealso \code{\link{llscoefbw}}, \code{\link{lcscoefgrad}}, \code{\link{spscsfa}}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n)
#' z1 <- runif(n)
#' z2 <- runif(n)
#' z3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' y <- cos(z1 + z2) + sin(z1 - z2) * x + rnorm(n, 0, 0.1)
#' dat <- data.frame(y, x, z1, z2, z3)
#'
#' bw <- llscoefbw(y ~ x | z1 + z2 + z3, data = dat)
#'
#' # estimation at the training data: the formula and the data
#' # are those of llscoefbw()
#' fit <- llscoef(bw, betas = TRUE, residuals = TRUE, gradients = TRUE)
#' fit
#' head(coef(fit))
#' head(residuals(fit))
#' head(fit$grad)
#' \donttest{
#' # smooth coefficients on a grid of z1, with z2 and z3 fixed
#' grid <- data.frame(x = mean(x), z1 = seq(0, 1, length.out = 50), z2 = 0.5,
#'                    z3 = factor("A", levels = levels(dat$z3)))
#' fit.grid <- llscoef(bw, newdata = grid, betas = TRUE)
#' plot(grid$z1, coef(fit.grid)[, "x"], xlab=expression(z[1]), ylab=expression(beta(z[1])), 
#'       type = "l")
#'
#' # bandwidths selected inside llscoef(), from a formula
#' fit.f <- llscoef(y ~ x | z1 + z2 + z3, data = dat, betas = TRUE)
#'
#' # given a bandwidth vector that is not an object of class 
#' # "llsmoothcoefficientbw": the formula is needed
#' fit.b <- llscoef(bw$bw, formula = y ~ x | z1 + z2 + z3, data = dat, betas = TRUE)
#' }
#' @export
llscoef <- function(bws, ...) {
  UseMethod("llscoef")
}

#' @rdname llscoef
#' @export
llscoef.llsmoothcoefficientbw <- function(bws, formula = NULL, data = NULL,
                                          newdata = NULL, betas = FALSE,
                                          residuals = FALSE, gradients = FALSE,
                                          ...) {
  .llscoef_unused(...)

  if (!is.null(formula)) {
    stop("'formula' must not be given when 'bws' is an object returned by ",
         "llscoefbw(): the object keeps the formula that its bandwidths were ",
         "selected for", call. = FALSE)
  }

  formula <- bws$formula

  if (!inherits(formula, "formula")) {
    stop("'bws' does not hold the formula of llscoefbw(): give a numeric vector ",
         "of local-linear bandwidths as 'bws' and a 'formula'", call. = FALSE)
  }

  if (is.null(data)) {
    data <- .llscoef_bw_data(bws)
  }

  # The na.action function used by llscoefbw() is kept in the bandwidth object
  .llscoef_fit(formula, bws, data, newdata, betas, residuals, gradients,
               call = match.call(), n_obs_bw = bws$n_obs,
               na.action = .llscoef_check_na_action(bws$na_action))
}

#' @rdname llscoef
#' @export
llscoef.formula <- function(bws, data = NULL, newdata = NULL, betas = FALSE,
                            residuals = FALSE, gradients = FALSE,
                            na.action = stats::na.omit, ...) {
  .llscoef_unused(...)

  na.action <- .llscoef_check_na_action(na.action)

  # given a formula, the bandwidths are selected first
  bw <- llscoefbw(bws, data = data, na.action = na.action)

  .llscoef_fit(bws, bw, data, newdata, betas, residuals, gradients,
               call = match.call(), n_obs_bw = bw$n_obs, na.action = na.action)
}

#' @rdname llscoef
#' @export
llscoef.default <- function(bws, formula, data = NULL, newdata = NULL,
                            betas = FALSE, residuals = FALSE,
                            gradients = FALSE, na.action = stats::na.omit,
                            ...) {
  .llscoef_unused(...)

  if (is.null(bws)) {
    stop("'bws' is required: an object returned by llscoefbw(), a formula or ",
         "a numeric vector of bandwidths", call. = FALSE)
  }

  # An object of class "scbandwidth" is used through its bandwidths if they
  # are those of the local-linear estimator, fixed bandwidths and the default kernels
  if (inherits(bws, "scbandwidth")) {
    bws <- .llscoef_np_bw(bws, "bws")
  }

  if (!is.numeric(bws)) {
    stop("'bws' must be an object returned by llscoefbw(), a formula or a ",
         "numeric vector of bandwidths", call. = FALSE)
  }

  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' such as y ~ x | z or y ~ x1 + x2 | z1 + z2 is required ",
         "when 'bws' is a vector of bandwidths (an object returned by ",
         "llscoefbw() with LSCV = TRUE keeps its formula and ",
         "accepts no other)",
         call. = FALSE)
  }

  .llscoef_fit(formula, bws, data, newdata, betas, residuals, gradients,
               call = match.call(),
               na.action = .llscoef_check_na_action(na.action))
}

# Estimation, given the formula, the training data and the bandwidths. 'bws' is
# a numeric vector or an object with a component bw; 'n_obs_bw' is the number of
# observations the bandwidths were computed with, if known.
.llscoef_fit <- function(formula, bws, data, newdata, betas, residuals,
                         gradients, call, n_obs_bw = NULL,
                         na.action = stats::na.omit) {

  betas <- .llscoef_flag(betas, "betas")
  residuals <- .llscoef_flag(residuals, "residuals")
  gradients <- .llscoef_flag(gradients, "gradients")

  # ----------------------------------------------------------
  # Training data: y, X = [1, x] and Z, extracted as in spscsfa()
  # ----------------------------------------------------------

  dat <- .llscoef_model_data(formula, data, na.action)
  X <- dat$X
  Y <- dat$Y
  Z <- dat$Z

  if (ncol(Z) < 1L) {
    stop("the nonparametric (z) part of 'formula' has no variables",
         call. = FALSE)
  }

  Z.is.cont <- vapply(Z, is.numeric, logical(1))
  Z.is.fac <- vapply(Z, is.factor, logical(1))

  if (any(!Z.is.cont & !Z.is.fac)) {
    stop("z variables must be numeric or factors: ",
         paste(names(Z)[!Z.is.cont & !Z.is.fac], collapse = ", "),
         call. = FALSE)
  }

  # ----------------------------------------------------------
  # Bandwidths
  # ----------------------------------------------------------

  if (is.list(bws) && !is.null(bws$bw)) {
    bws <- bws$bw
  }

  if (!is.null(n_obs_bw) && !identical(as.integer(n_obs_bw), nrow(Z))) {
    warning("the training data have ", nrow(Z), " observations but the ",
            "bandwidths were computed with ", n_obs_bw, ": check 'data'",
            call. = FALSE)
  }

  # Named bandwidths are matched to the z variables by name; a warning is
  # given if their order is not that of the z variables in 'formula'
  bws <- .llscoef_match_bw_names(bws, names(Z))

  bws <- .llscoef_check_bw(bws, Z, .llscoef_bw_max(Z))

  # ------------------
  # Evaluation points
  # ------------------

  has.newdata <- !is.null(newdata)

  # Rows dropped because of missing values of training data or
  # evaluation data if given
  omit <- dat$omit
  eval.na <- dat$na.action

  if (has.newdata) {
    ev <- .llscoef_newdata(formula, newdata, na.action)
    X.eval <- ev$X
    Z.eval <- ev$Z
    omit <- ev$omit
    eval.na <- ev$na.action

    if (nrow(Z.eval) < 1L) {
      stop("'newdata' has no rows without missing values", call. = FALSE)
    }
    if (ncol(X.eval) != ncol(X)) {
      stop("'newdata' gives a different number of x columns than the training ",
           "data (check the levels of factor x variables)", call. = FALSE)
    }
  } else {
    X.eval <- X
    Z.eval <- Z
  }

  # ----------------------------------------------------------
  # Estimation: [beta(z), gradients(z)] at the evaluation points
  # ----------------------------------------------------------

  k <- ncol(X)
  p <- sum(Z.is.cont)

  est <- sscm.ll.mixed(X, Y, Z, bws, Z.eval)

  beta.eval <- est[, seq_len(k), drop = FALSE]
  mean <- as.vector(rowSums(X.eval * beta.eval))

  # Residuals: training data, Y minus the fitted values at the training points
  resid <- NULL

  if (residuals) {
    est.train <- if (has.newdata) sscm.ll.mixed(X, Y, Z, bws) else est
    resid <- as.vector(Y - rowSums(X * est.train[, seq_len(k), drop = FALSE]))
  }

  coef.names <- c("Intercept", colnames(X)[-1L])

  beta <- NULL

  if (betas) {
    beta <- beta.eval
    dimnames(beta) <- list(NULL, coef.names)
  }

  # Gradients with respect to the continuous z variables, X-major order
  grad <- NULL

  if (gradients) {
    grad <- est[, k + seq_len(k * p), drop = FALSE]
    grad.names <- character(0)
    if (p > 0L) {
      grad.names <- as.vector(outer(
        names(Z)[Z.is.cont],
        coef.names,
        function(z, b) paste0("d(", b, ")/d(", z, ")")
      ))
    }
    dimnames(grad) <- list(NULL, grad.names)
  }

  # With na.exclude the results are padded with NA at the rows that were
  # dropped, so that they line up with either evaluation or training data. 
  # The residuals are those of the training data.
  mean <- stats::napredict(eval.na, mean)

  if (!is.null(beta)) {
    beta <- stats::napredict(eval.na, beta)
  }

  if (!is.null(grad)) {
    grad <- stats::napredict(eval.na, grad)
  }

  if (!is.null(resid)) {
    resid <- stats::naresid(dat$na.action, resid)
  }

  structure(
    list(
      mean = mean,
      beta = beta,
      resid = resid,
      grad = grad,
      bw = stats::setNames(bws, names(Z)),
      x_names = coef.names[-1L],
      z_names = names(Z),
      n_train = nrow(Z),
      n_eval = nrow(Z.eval),
      eval_train = !has.newdata,
      betas = betas,
      residuals = residuals,
      gradients = gradients,
      omit = omit,
      call = call
    ),
    class = "llsmoothcoefficient"
  )
}

#' @rdname llscoef
#' @export
print.llsmoothcoefficient <- function(x, ...) {
  cat("\nLocal-Linear Smooth Coefficient Model",
      "\nRegression data: ", x$n_train, " training points",
      if (x$eval_train) "" else paste0(" and ", x$n_eval, " evaluation points"),
      ", in ", length(x$z_names), " smooth variable(s)\n", sep = "")
  print(stats::setNames(x$bw, x$z_names))
  cat("\n")
  invisible(x)
}

#' @rdname llscoef
#' @export
coef.llsmoothcoefficient <- function(object, ...) {
  if (is.null(object$beta)) {
    stop("the smooth coefficients were not computed: use betas = TRUE",
         call. = FALSE)
  }
  object$beta
}

#' @rdname llscoef
#' @export
residuals.llsmoothcoefficient <- function(object, ...) {
  if (is.null(object$resid)) {
    stop("the residuals were not computed: use residuals = TRUE",
         call. = FALSE)
  }
  object$resid
}

#' @rdname llscoef
#' @export
fitted.llsmoothcoefficient <- function(object, ...) {
  object$mean
}


.llscoef_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
  }
  x
}

# Training data of the call to llscoefbw() that produced 'bws', evaluated again
# in the environment the call was made from. NULL if that call had no 'data'; 
# the variables are then taken from the environment of the formula.
.llscoef_bw_data <- function(bws) {

  data.expr <- bws$call[["data"]]

  if (is.null(data.expr)) {
    return(NULL)
  }

  env <- environment(bws$call)

  if (is.null(env)) {
    env <- environment(bws$formula)
  }

  if (is.null(env)) {
    env <- globalenv()
  }

  tryCatch(
    eval(data.expr, envir = env),
    error = function(e) {
      stop("could not find the training data '",
           paste(deparse(data.expr), collapse = ""),
           "' of the call to llscoefbw() (", conditionMessage(e),
           "): supply them with 'data'", call. = FALSE)
    }
  )
}

# An object returned by np::npscoefbw() is accepted by llscoef() only if it 
# was made with regtype = "ll", fixed bandwidths (not scale factors) and the default kernels.
.llscoef_np_bw <- function(bw_obj, arg) {
  regtype <- .spsc_np_regtype(bw_obj)

  if (!identical(regtype, "ll")) {
    stop("'", arg, "' is an object returned by npscoefbw() with regtype = \"",
         regtype, "\"",
         if (identical(regtype, "lc")) {
           ", whose bandwidths are for the local-constant estimator"
         },
         ": the local-linear estimator needs an object returned by ",
         "npscoefbw() with regtype = \"ll\", an object returned by ",
         "llscoefbw() or a numeric vector of bandwidths", call. = FALSE)
  }

  .spsc_validate_np_bw_settings(bw_obj, arg, "local-linear")

  bw <- bw_obj[["bw"]]

  if (!is.numeric(bw)) {
    stop("'", arg, "' does not hold numeric bandwidths", call. = FALSE)
  }

  bw
}

# 'llscoef()' has no further arguments: supplying any is an error, as for a
# function without '...'
.llscoef_unused <- function(...) {
  args <- list(...)

  if (length(args) > 0L) {
    nms <- names(args)
    nms <- nms[!is.na(nms) & nzchar(nms)]
    stop("unused argument(s)",
         if (length(nms) > 0L) paste0(": ", paste(nms, collapse = ", ")),
         call. = FALSE)
  }

  invisible(NULL)
}

# Bandwidths with names are matched to the z variables by name, in any order,
# as in the result of llscoefbw(), whose names are those of the z variables in
# the order of the formula given to llscoefbw(). A warning is given if that
# order differs from the order of the z variables in the formula given to
# llscoef(). Bandwidths without names are taken in the order of the z
# variables.
.llscoef_match_bw_names <- function(bws, znames) {

  bw.names <- names(bws)

  if (is.null(bw.names)) {
    return(bws)
  }

  if (anyNA(bw.names) || any(!nzchar(bw.names)) || anyDuplicated(bw.names)) {
    stop("'bws' has missing, empty or duplicated names: name every bandwidth ",
         "after a z variable (", paste(znames, collapse = ", "),
         ") or give no names", call. = FALSE)
  }

  not.z <- setdiff(bw.names, znames)
  no.bw <- setdiff(znames, bw.names)

  if (length(not.z) > 0L || length(no.bw) > 0L) {
    detail <- c(
      if (length(not.z) > 0L) {
        paste0("not a z variable: ", paste(not.z, collapse = ", "))
      },
      if (length(no.bw) > 0L) {
        paste0("no bandwidth for: ", paste(no.bw, collapse = ", "))
      }
    )
    stop("names of 'bws' do not match the z variables (",
         paste(znames, collapse = ", "), "); ",
         paste(detail, collapse = "; "), call. = FALSE)
  }

  # The names are those of the z variables, so only their order can differ
  if (!identical(as.character(bw.names), as.character(znames))) {
    warning("the z variables in 'formula' (", paste(znames, collapse = ", "),
            ") are in a different order from the one the bandwidths were ",
            "computed for (", paste(bw.names, collapse = ", "),
            "); the bandwidths were matched to the z variables by name",
            call. = FALSE)
  }

  bws[znames]
}

# X = [1, x] and Z at the rows of 'newdata'. The variables are looked up in
# 'newdata' by the x and z parts of 'formula'; the response is not needed.
# Rows with missing values are dropped from 'na.action' ('omit' gives them,
# 'na.action' is the object that napredict() uses to pad results).
.llscoef_newdata <- function(formula, newdata, na.action = stats::na.omit) {

  parsed_formula <- Formula::Formula(formula)
  x_formula <- stats::formula(parsed_formula, lhs = 0, rhs = 1)
  z_formula <- stats::formula(parsed_formula, lhs = 0, rhs = 2)

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  x_mf <- stats::model.frame(x_formula, data = newdata,
                             na.action = stats::na.pass)
  z_mf <- stats::model.frame(z_formula, data = newdata,
                             na.action = stats::na.pass)

  x <- stats::model.matrix(attr(x_mf, "terms"), data = x_mf)

  if (ncol(x) > 0L && colnames(x)[1L] == "(Intercept)") {
    x <- x[, -1L, drop = FALSE]
  }

  # ensure that the x and z variables use exactly the same observations 
  # after applying the specified na.action rule
  ev.df <- if (ncol(x_mf) > 0L) cbind(x_mf, z_mf) else z_mf
  na.obj <- attr(match.fun(na.action)(ev.df), "na.action")

  keep <- rep_len(TRUE, nrow(z_mf))
  if (length(na.obj) > 0L) {
    keep[as.integer(na.obj)] <- FALSE
  }

  Z <- z_mf
  attr(Z, "terms") <- NULL

  X <- cbind(rep(1, nrow(z_mf)), x)[keep, , drop = FALSE]
  Z <- Z[keep, , drop = FALSE]

  # na.pass leaves missing values in the data, which the estimation cannot use
  if (anyNA(X) || any(vapply(Z, anyNA, logical(1)))) {
    stop("missing values remain in 'newdata' after 'na.action': use na.omit, ",
         "na.exclude or na.fail", call. = FALSE)
  }

  list(X = X, Z = Z, omit = as.integer(na.obj), na.action = na.obj)
}

# ==============================================================================
# Local-linear estimation at the evaluation points
# ==============================================================================

# Arguments
#   X, Y     training data (n.train rows); 
#            X includes the column of ones
#   Z.train  training data of the variables the coefficients depend on
#   bw       bandwidths, one for each column of Z.train
#   Z.eval   evaluation data; it may differ from the training data. It needs the
#            same columns as Z.train (factors may have fewer levels)
# Value
#   matrix with one row for each row of Z.eval: the estimates of the
#   coefficients (first ncol(X) columns) followed by their gradients with
#   respect to the continuous variables in Z (X-major order)

sscm.ll.mixed <- function(X, Y, Z.train, bw, Z.eval = Z.train) {

    Z.train <- as.data.frame(Z.train)
    Z.eval <- as.data.frame(Z.eval)

    n.train <- nrow(Z.train)
    n.eval <- nrow(Z.eval)

    if (NROW(X) != n.train || NROW(Y) != n.train) {
      stop("X and Y must have the same number of rows as Z.train (the training data)")
    }
    if (ncol(Z.eval) != ncol(Z.train)) {
      stop("Z.eval must have the same columns as Z.train")
    }

    Z.col.numeric <- sapply(1:ncol(Z.train),function(i){is.numeric(Z.train[,i])})

    # factors of Z.eval come from the levels of the training factors
    for (i in 1:ncol(Z.train)) {
      if (Z.col.numeric[i]) {
        if (!is.numeric(Z.eval[,i])) {
          stop(paste("column", i, "of Z.eval must be numeric, as in Z.train"))
        }
      } else if (is.factor(Z.train[,i])) {
        Z.eval.i <- factor(as.character(Z.eval[,i]), levels=levels(Z.train[,i]),
                           ordered=is.ordered(Z.train[,i]))
        if (any(is.na(Z.eval.i) & !is.na(Z.eval[,i]))) {
          stop(paste("column", i, "of Z.eval has levels that are not in Z.train"))
        }
        Z.eval[,i] <- Z.eval.i
      }
    }

    Z_c <- as.data.frame(Z.train[,Z.col.numeric])
    Z_c.eval <- as.data.frame(Z.eval[,Z.col.numeric])
    p <- ncol(Z_c)

    beta <- rbind()
    for (j in 1:n.eval) {
      Kvec <- np::npksum(bws=bw, txdat=Z.eval[j,,drop=FALSE], exdat=Z.train)$ksum

      S <- cbind(X)
      if (p > 0) {
        Zmat <- as.matrix(Z_c-sapply(Z_c.eval[j,], rep, n.train))
        for (k in 1:ncol(X)) {S <- cbind(S, X[,k]*Zmat)}
      }

      SK <- S*as.vector(Kvec)
      A <- crossprod(SK, S)
      B <- crossprod(SK, Y)

      # Invert A; if that fails, use ordinary diagonal ridging and report
      # the required ridge with a warning.
      A.inv <- tryCatch(solve(A), error = function(e) NULL)
      
      if (is.null(A.inv)) {
        ridging <- .ordinary_diagonal_ridging(
          A = A,
          n.train = n.train,
          n.diag = ncol(A),
          eval = j
        )
        A.inv <- ridging$A.inv
        ridge <- ridging$ridge

        if (is.null(A.inv)) {
          # using LSCV bandwidth, this should not happen
          stop("could not invert the weighted cross-product matrix at ",
               "evaluation point ", j, call. = FALSE)
        }
      }
      
      beta.temp <- as.vector(A.inv %*% B)

      beta <- rbind(beta, t(beta.temp))

    }

    beta

}
