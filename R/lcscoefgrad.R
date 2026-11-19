#' Gradients of the local-constant smooth coefficient estimator
#'
#' @description
#' Computes the analytic gradients, with respect to the continuous \code{z}
#' variables, of the smooth coefficients of the local-constant (LC)
#' semiparametric smooth coefficient (SPSC) model
#' \eqn{y_i = x_i'\beta(z_i) + u_i} following Geng and Sun (2019). It is called by
#' \code{\link{spscsfa}} or can be used independently. An \code{np::npscoefbw} 
#' bandwidth object with its default LC regression type and kernels can be passed 
#' into the \code{bws} argument of this function.
#'
#' @details
#' The LC estimator is
#' \deqn{\hat\beta(z) = \Big(\sum_i x_i x_i' K_i(z)\Big)^{-1}
#'   \sum_i x_i y_i K_i(z).}
#' Differentiating the locally weighted moment condition
#' \eqn{\sum_i x_i \hat u_i(z) K_i(z) = 0}, where
#' \eqn{\hat u_i(z) = y_i - x_i'\hat\beta(z)}, with respect to the \eqn{l}th
#' continuous \code{z} variable gives
#' \deqn{D_l\hat\beta(z) = \Big(\sum_i x_i x_i' K_i(z)\Big)^{-1}
#'   \sum_i x_i \hat u^*_{li}(z) K_i(z),}
#' that is, another LC regression, of the transformed residuals
#' \eqn{\hat u^*_{li}(z) = \hat u_i(z)(z_{li} - z_l)/h_l^2}, on the original
#' regressors. This is the result for a Gaussian kernel for the continuous
#' variables: \eqn{k\left(\frac{z_{li}-z_l}{h_l}\right) = \frac{1}{\sqrt{2\pi}}
#' \exp\left(-\frac{1}{2} \left(\frac{z_{li}-z_l}{h_l}\right)^2 \right)}. 
#' Factors enter the product kernel \eqn{K_i(z)} through kernels that
#' do not depend on the continuous variables, so they do not change the formula.
#'
#' The kernel weights are those of \code{np::npksum}. A
#' bandwidth object from \code{np::npscoefbw} must therefore be for
#' \code{regtype = "lc"}, with fixed bandwidths (not scale factors) and a
#' second-order Gaussian kernel without bounds for the continuous variables;
#' its kernels for factors are used as they are. If the weighted cross-product
#' matrix at an evaluation point is singular, a ridge is added to it and a warning 
#' is given.
#'
#' @param bws The bandwidths: the bandwidth object returned by
#'   \code{np::npscoefbw} (as \code{spscsfa} holds it), an object with a component
#'   \code{bw}, or a numeric vector with one bandwidth for each column of
#'   \code{tzdat}. A numeric vector is used with the default kernels of
#'   \code{np::npscoefbw}: second-order Gaussian for continuous variables, Aitchison-Aitken for
#'   unordered and Li-Racine for ordered factors. If it has names that are those
#'   of the columns of \code{tzdat}, it is matched to them by name; otherwise it
#'   is taken in their order.
#' @param txdat A data frame or matrix with the regressors \code{x}. The
#'   constant is added.
#' @param tydat A numeric vector with the dependent variable.
#' @param tzdat A data frame with the variables the coefficients depend on:
#'   numeric (continuous), factors (unordered) or ordered factors.
#' @param ezdat Optional data frame with the evaluation points, with the same
#'   variables as \code{tzdat}. If \code{NULL} (default), the gradients are
#'   evaluated at the rows of \code{tzdat}.
#'
#' @return A named list with one element for each continuous variable in
#'   \code{tzdat}, in the order of its columns and named after them. The
#'   element for the \eqn{l}th continuous variable is a matrix with one row for
#'   each evaluation point and one column for each smooth coefficient, in the
#'   order of the coefficients (the constant first, then the columns of
#'   \code{txdat}). Its columns are the gradients of all smooth coefficients
#'   with respect to \code{z_l} and are named as \code{llscoef} names its gradients:
#'   \code{d(Intercept)/d(z_l)}, \code{d(x1)/d(z_l)}, \ldots. The list is empty if
#'   there is no continuous variable in \code{tzdat}.
#'
#' @references
#' Xin Geng and Kai Sun (2019). Gradient estimation of the local-constant
#' semiparametric smooth coefficient model. \emph{Economics Letters}, 185,
#' 108684.
#'
#' @seealso \code{\link{llscoef}}, whose gradients are those of the
#'   local-linear estimator.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' z <- runif(n)
#' x <- rnorm(n)
#' y <- sin(2 * pi * z) - cos(2 * pi * z) * x + rnorm(n, sd = 0.5)
#' dat <- data.frame(y, x, z)
#'
#' bw <- np::npscoefbw(y ~ x | z, data = dat)
#' g <- lcscoefgrad(bw, txdat = dat["x"], tydat = dat$y, tzdat = dat["z"])
#'
#' # g$z has the columns d(Intercept)/d(z) and d(x)/d(z); the true gradients
#' # are 2 * pi * cos(2 * pi * z) and 2 * pi * sin(2 * pi * z)
#' head(g$z)
#'
#' \dontrun{
#' 
#' # gradients on a grid of z, through 'ezdat': plotted against z, with the
#' # true gradients (dashed)
#' ez <- data.frame(z = seq(0, 1, length.out = 100))
#' g.ez <- lcscoefgrad(bw, txdat = dat["x"], tydat = dat$y, tzdat = dat["z"],
#'                     ezdat = ez)
#'
#' op <- par(mfrow = c(1, 2))
#' plot(ez$z, g.ez$z[, "d(Intercept)/d(z)"], type = "l",
#'      xlab = expression(z),
#'      ylab = expression(partialdiff * beta[0](z) / partialdiff * z))
#' lines(ez$z, 2 * pi * cos(2 * pi * ez$z), lty = 2)
#' plot(ez$z, g.ez$z[, "d(x)/d(z)"], type = "l",
#'      xlab = expression(z),
#'      ylab = expression(partialdiff * beta[1](z) / partialdiff * z))
#' lines(ez$z, 2 * pi * sin(2 * pi * ez$z), lty = 2)
#' par(op)
#' }
#'
#' @export
lcscoefgrad <- function(bws, txdat, tydat, tzdat, ezdat = NULL) {

  # ----------------------------------------------------------
  # Training data: Y, X = [1, x] and Z
  # ----------------------------------------------------------

  X0 <- as.matrix(txdat)

  if (!is.numeric(X0)) {
    stop("'txdat' must be numeric: dummy-code factor regressors first",
         call. = FALSE)
  }

  Y <- as.matrix(tydat)

  if (!is.numeric(Y) || ncol(Y) != 1L) {
    stop("'tydat' must be a numeric vector", call. = FALSE)
  }

  Y <- as.vector(Y)
  Z.train <- as.data.frame(tzdat)
  n.train <- length(Y)

  if (nrow(X0) != n.train || nrow(Z.train) != n.train) {
    stop("'txdat', 'tydat' and 'tzdat' must have the same number of rows",
         call. = FALSE)
  }

  if (ncol(Z.train) < 1L) {
    stop("'tzdat' has no variables", call. = FALSE)
  }

  Z.is.cont <- vapply(Z.train, is.numeric, logical(1))
  Z.is.fac <- vapply(Z.train, is.factor, logical(1))

  if (any(!Z.is.cont & !Z.is.fac)) {
    stop("z variables must be numeric or factors: ",
         paste(names(Z.train)[!Z.is.cont & !Z.is.fac], collapse = ", "),
         call. = FALSE)
  }

  if (anyNA(X0) || anyNA(Y) || any(vapply(Z.train, anyNA, logical(1)))) {
    stop("'txdat', 'tydat' and 'tzdat' must not have missing values",
         call. = FALSE)
  }

  x.names <- colnames(X0)

  if (is.null(x.names)) {
    x.names <- paste0("x", seq_len(ncol(X0)))
  }

  coef.names <- c("Intercept", x.names)
  X <- cbind(1, X0)
  p <- ncol(X)

  # ----------------------------------------------------------
  # Evaluation points: the training data unless 'ezdat' is given
  # ----------------------------------------------------------

  if (is.null(ezdat)) {
    Z.eval <- Z.train
  } else {
    Z.eval <- as.data.frame(ezdat)

    if (ncol(Z.eval) != ncol(Z.train)) {
      stop("'ezdat' must have the same variables as 'tzdat'", call. = FALSE)
    }

    if (!identical(names(Z.eval), names(Z.train)) &&
        setequal(names(Z.eval), names(Z.train))) {
      Z.eval <- Z.eval[names(Z.train)]
    }

    if (any(vapply(Z.eval, anyNA, logical(1)))) {
      stop("'ezdat' must not have missing values", call. = FALSE)
    }

    # Factors take the levels of the training data
    for (i in seq_len(ncol(Z.train))) {
      if (Z.is.fac[i]) {
        ev <- factor(as.character(Z.eval[[i]]),
                     levels = levels(Z.train[[i]]),
                     ordered = is.ordered(Z.train[[i]]))

        if (anyNA(ev)) {
          stop("column ", i, " of 'ezdat' has levels that are not in 'tzdat'",
               call. = FALSE)
        }

        Z.eval[[i]] <- ev
      } else if (is.numeric(Z.eval[[i]])) {
        Z.eval[[i]] <- as.numeric(Z.eval[[i]])
      } else {
        stop("column ", i, " of 'ezdat' must be numeric, as in 'tzdat'",
             call. = FALSE)
      }
    }
  }

  n.eval <- nrow(Z.eval)

  # ----------------------------------------------------------
  # Bandwidths and kernels
  # ----------------------------------------------------------

  ukertype <- "aitchisonaitken"
  okertype <- "liracine"
  bw <- bws

  if (inherits(bws, "scbandwidth")) {
    # Bandwidth object of class "scbandwidth": the gradient below is that of the
    # local-constant estimator with fixed bandwidths and a second-order
    # Gaussian kernel for the continuous variables
    if (!identical(bws$regtype, "lc")) {
      stop("'bws' is for regtype = \"", bws$regtype, "\", but the gradients ",
           "are those of the local-constant estimator (regtype = \"lc\")",
           call. = FALSE)
    }

    if (!identical(bws$type, "fixed") || isTRUE(bws$scaling) ||
        !identical(bws$ckertype, "gaussian") ||
        !identical(as.numeric(bws$ckerorder), 2) ||
        !(is.null(bws$ckerbound) || identical(bws$ckerbound, "none"))) {
      stop("'bws' must have fixed bandwidths (not scale factors) and a ",
           "second-order Gaussian kernel without bounds for the continuous ",
           "z variables", call. = FALSE)
    }

    ukertype <- bws$ukertype
    okertype <- bws$okertype
    bw <- bws$bw
  } else if (is.list(bws) && !is.null(bws$bw)) {
    bw <- bws$bw
  }

  if (!is.numeric(bw)) {
    stop("'bws' must be a bandwidth object or a numeric vector of bandwidths",
         call. = FALSE)
  }

  # Named bandwidths are matched to the z variables by name
  if (!is.null(names(bw)) && !anyDuplicated(names(bw)) &&
      setequal(names(bw), names(Z.train))) {
    bw <- bw[names(Z.train)]
  }

  bw <- as.numeric(bw)

  if (length(bw) != ncol(Z.train)) {
    stop("'bws' has ", length(bw), " bandwidth(s) but 'tzdat' has ",
         ncol(Z.train), " variable(s)", call. = FALSE)
  }

  if (any(!is.finite(bw)) || any(bw < 0) || any(bw[Z.is.cont] == 0)) {
    stop("bandwidths must be positive for continuous variables and ",
         "non-negative for factors", call. = FALSE)
  }

  # ----------------------------------------------------------
  # Gradients, one evaluation point at a time
  # ----------------------------------------------------------

  cont <- which(Z.is.cont)
  n.cont <- length(cont)
  z.names <- names(Z.train)[cont]

  if (n.cont == 0L) {
    return(stats::setNames(list(), character(0)))
  }

  Z.train.cont <- as.matrix(Z.train[cont])
  Z.eval.cont <- as.matrix(Z.eval[cont])
  h2 <- bw[cont]^2

  # One matrix for each continuous z variable
  grad <- lapply(seq_len(n.cont), function(l) matrix(NA_real_, n.eval, p))

  for (j in seq_len(n.eval)) {

    # K_i(z_j) for all training points i: the training point of npksum() is
    # the evaluation point and its evaluation points are the training data
    Kvec <- as.vector(np::npksum(bws = bw,
                                 txdat = Z.eval[j, , drop = FALSE],
                                 exdat = Z.train,
                                 ukertype = ukertype,
                                 okertype = okertype)$ksum)

    if (!all(is.finite(Kvec))) {
      stop("non-finite kernel weights at evaluation point ", j, call. = FALSE)
    }

    XK <- X * Kvec
    A <- crossprod(XK, X)
    B <- crossprod(XK, Y)

    # Invert A; if that fails, use ordinary diagonal ridging and report
    # the required ridge with a warning.
    A.inv <- tryCatch(solve(A), error = function(e) NULL)

    if (is.null(A.inv)) {
      ridging <- .ordinary_diagonal_ridging(
        A = A,
        n.train = n.train,
        n.diag = p,
        eval = j
      )
      A.inv <- ridging$A.inv
      ridge <- ridging$ridge

      if (is.null(A.inv)) {
        stop("could not invert the weighted cross-product matrix at ",
             "evaluation point ", j, call. = FALSE)
      }
    }

    beta <- as.vector(A.inv %*% B)
    uhat <- Y - as.vector(X %*% beta)

    # u*_li = uhat_i (z_li - z_l) / h_l^2: one column for each continuous z
    Ustar <- uhat * sweep(sweep(Z.train.cont, 2L, Z.eval.cont[j, ], "-"),
                          2L, h2, "/")

    G <- A.inv %*% crossprod(XK, Ustar)

    for (l in seq_len(n.cont)) {
      grad[[l]][j, ] <- G[, l]
    }
  }

  names(grad) <- z.names

  for (l in seq_len(n.cont)) {
    dimnames(grad[[l]]) <- list(
      NULL,
      paste0("d(", coef.names, ")/d(", z.names[l], ")")
    )
  }

  grad
}
