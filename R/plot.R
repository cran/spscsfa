#' Plot smooth coefficients of an SPSC-SFA model
#'
#' Plot method for objects of class \code{"spscsfa"}. Each estimated smooth
#' coefficient, i.e. the intercept \eqn{\theta(Z)} and the slope coefficients
#' \eqn{\beta(Z)}, is plotted against one continuous nonparametric variable
#' \eqn{Z} at a time, over a grid running from \code{min(Z)} to \code{max(Z)}.
#' The smooth coefficients are those of the estimator used to fit the model:
#' local-constant if \code{spscreg = "lc"} and local-linear if
#' \code{spscreg = "ll"} in \code{\link{spscsfa}}.
#'
#' @details
#' For every continuous variable in the nonparametric part \eqn{Z} of the model
#' (the second part of the formula), the smooth coefficients are evaluated on
#' \code{n_grid} equally spaced points between the sample minimum and maximum
#' of that variable. The grid is supplied to \code{\link{spscsfa}} through its
#' \code{ezdat} argument, re-using the bandwidth of the fitted model, so no new
#' cross-validation is carried out. The same estimator (\code{spscreg}) as in
#' the fitted model is used for the grid: the local-constant estimator for \code{spscreg = "lc"} 
#' and the local-linear estimator for \code{spscreg = "ll"}.
#'
#' If the model has more than one \eqn{Z} variable, the other continuous
#' \eqn{Z} variables are held fixed at their sample averages while the variable
#' on the horizontal axis runs over its grid. One page is drawn for each
#' continuous \eqn{Z} variable, with one panel for each smooth coefficient.
#'
#' Discrete \eqn{Z} variables (factors and ordered factors) are ignored by the
#' plot method: no plots are drawn against them. The smooth coefficients still
#' depend on them, however, so they have to take some value when the
#' coefficients are evaluated; each discrete \eqn{Z} variable is held at its
#' most frequent category.
#'
#' A term in the second part of the formula is a variable as it enters the
#' model: for \code{log(z)} the grid runs over the values of \code{log(z)},
#' which also labels the axis, and \code{factor(z2)} is a discrete variable.
#'
#' Only the smooth coefficients of the first estimation step, i.e. the columns
#' of \code{smoothcoef(x)$theta_beta}, are plotted. In particular, the
#' intercept is \eqn{\theta(Z)} (\code{Intercept (theta)}) and not the
#' recovered intercept \eqn{\alpha(Z)}.
#'
#' @param x An object of class \code{"spscsfa"} returned by
#'   \code{\link{spscsfa}}.
#'
#' @param which Smooth coefficients to plot, given as positions or names of the
#'   columns of \code{smoothcoef(x)$theta_beta}. The default (\code{NULL})
#'   plots all of them: the intercept followed by one slope coefficient for
#'   each variable in the parametric part of the formula.
#'
#' @param zvar Continuous \eqn{Z} variables to plot against, given as names or
#'   as positions in the nonparametric part of the formula. The default
#'   (\code{NULL}) uses all continuous \eqn{Z} variables.
#'
#' @param n_grid Number of equally spaced grid points between \code{min(Z)} and
#'   \code{max(Z)} (default 100).
#'
#' @param rug Logical; if \code{TRUE} (default), the observed values of \eqn{Z}
#'   are marked by a rug along the horizontal axis.
#'
#' @param mfrow Integer vector of length 2 giving the number of rows and
#'   columns of panels per page. The default (\code{NULL}) arranges all smooth
#'   coefficients against one \eqn{Z} variable on a single, roughly square page. The
#'   graphical parameters are restored on exit.
#'
#' @param ask Logical; if \code{TRUE}, the user is asked before each new page is
#'   drawn. The default asks only if more than one page is drawn on an
#'   interactive device.
#'
#' @param ... Further graphical parameters passed to
#'   \code{\link[graphics]{plot}}, for example \code{col}, \code{lty},
#'   \code{lwd}, \code{ylim} or \code{main}. By default lines are drawn
#'   (\code{type = "l"}) and the axes are labelled with the name of the
#'   \eqn{Z} variable and of the smooth coefficient.
#'
#' @return Invisibly, a named list with one data frame for each \eqn{Z}
#'   variable plotted against. The first column of each data frame is the grid
#'   of that \eqn{Z} variable and the remaining columns are the smooth
#'   coefficients evaluated on the grid, so that the curves can be redrawn or
#'   customized.
#'
#' @seealso \code{\link{spscsfa}}, \code{\link{smoothcoef}}
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' z <- runif(n)
#' y <- cos(z) + sin(z) * x + rnorm(n, 0, 0.1) -
#'      exp(0.1 + 0.1 * z) * abs(rnorm(n, 0, 1))
#' dat <- data.frame(output = y, input = x, covariate = z)
#' result <- spscsfa(output ~ input | covariate, data = dat, verbose = FALSE)
#'
#' # the intercept and the slope coefficient plotted against 'covariate'
#' plot(result)
#'
#' # the plotted values are returned invisibly
#' curves <- plot(result, n_grid = 50, col = "blue")
#' head(curves$covariate)
#' \donttest{
#' # Two continuous z variables and one discrete z variable: one page per
#' # continuous z variable; the other continuous z variable is held at its
#' # average and the factor at its most frequent category
#' z1 <- runif(n)
#' z2 <- runif(n)
#' z3 <- factor(rep(c("A", "B"), length.out = n))
#' y <- cos(z1 + z2) + sin(z1 - z2) * x + rnorm(n, 0, 0.1) -
#'      exp(0.1 + 0.1 * z1) * abs(rnorm(n, 0, 1))
#' dat <- data.frame(y, x, z1, z2, z3)
#' result <- spscsfa(y ~ x | z1 + z2 + z3, data = dat, verbose = FALSE)
#' plot(result)
#'
#' # the slope coefficient only, against z2 only
#' plot(result, which = 2, zvar = "z2")
#'
#' # local-linear smooth coefficients (spscreg = "ll"), plotted the same way
#' result_ll <- spscsfa(y ~ x | z1 + z2 + z3, data = dat, spscreg = "ll", verbose = FALSE)
#' plot(result_ll)
#' }
#'
#' @importFrom graphics par plot rug
#' @importFrom grDevices dev.interactive devAskNewPage
#' @export
plot.spscsfa <- function(x,
                          which = NULL,
                          zvar = NULL,
                          n_grid = 100L,
                          rug = TRUE,
                          mfrow = NULL,
                          ask = NULL,
                          ...) {

  if (!.is_spscsfa(x)) stop("x must be of class 'spscsfa'")

  if (!is.null(mfrow) &&
      (!is.numeric(mfrow) || length(mfrow) != 2L || anyNA(mfrow) ||
       any(mfrow < 1))) {
    stop("'mfrow' must be a vector of two positive integers")
  }

  # smooth coefficients on the grids, obtained through 'ezdat' of spscsfa()
  pd <- spsc_plot_data(x, which = which, zvar = zvar, n_grid = n_grid)
  zs <- names(pd$curves)
  cf <- pd$coef_names

  # one page per z variable, one panel per smooth coefficient
  if (is.null(mfrow)) {
    n_col <- ceiling(sqrt(length(cf)))
    mfrow <- c(ceiling(length(cf) / n_col), n_col)
  }

  op <- par(mfrow = mfrow)
  on.exit(par(op), add = TRUE)

  n_pages <- length(zs) * ceiling(length(cf) / prod(mfrow))
  if (is.null(ask)) ask <- n_pages > 1L && dev.interactive()
  if (isTRUE(ask)) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask), add = TRUE)
  }

  dots <- list(...)
  for (zn in zs) {
    d <- pd$curves[[zn]]
    z_obs <- x$model_frame[[zn]]
    for (j in seq_along(cf)) {
      defaults <- list(type = "l", lwd = 2, xlab = zn, ylab = cf[j])
      args <- c(list(d[[1L]], d[[j + 1L]]), dots,
                defaults[setdiff(names(defaults), names(dots))])
      do.call(plot, args)
      if (isTRUE(rug)) rug(z_obs)
    }
  }

  invisible(pd$curves)
}

#' @method plot spsc_sfa
#' @export
#' @noRd
plot.spsc_sfa <- plot.spscsfa

## Internal helper: evaluate the smooth coefficients of a fitted model on
## one-dimensional grids of its continuous z variables.
##
## For each continuous z variable plotted, the evaluation points are an equally
## spaced grid from min(z) to max(z). All other continuous z variables are held
## at their means and discrete z variables at their most frequent category. The
## evaluation points of all z variables are stacked into a single 'ezdat' data
## frame, so that spscsfa() is called only once.
spsc_plot_data <- function(x, which = NULL, zvar = NULL, n_grid = 100L) {

  ms <- x$model_summary
  mf <- x$model_frame
  z_names <- ms$z_names

  if (is.null(mf) || !all(z_names %in% names(mf))) {
    stop("'x' does not contain the estimation data ('model_frame') needed ",
         "for plotting; re-estimate the model with the current version of ",
         "spscsfa()", call. = FALSE)
  }

  zdat <- mf[, z_names, drop = FALSE]
  is_cont <- vapply(zdat, is.numeric, logical(1L))
  cont <- z_names[is_cont]

  if (length(cont) == 0L) {
    stop("the model has no continuous nonparametric (z) variable to plot ",
         "against; discrete z variables are ignored by the plot method",
         call. = FALSE)
  }

  # z variables on the horizontal axis
  if (is.null(zvar)) {
    zsel <- cont
  } else {
    zsel <- if (is.character(zvar)) zvar else z_names[zvar]
    if (anyNA(zsel) || !all(zsel %in% z_names)) {
      stop("'zvar' must be names or positions of z variables of the model: ",
           paste(z_names, collapse = ", "), call. = FALSE)
    }
    if (!all(zsel %in% cont)) {
      stop("discrete z variables are ignored by the plot method: ",
           paste(setdiff(zsel, cont), collapse = ", "), call. = FALSE)
    }
    zsel <- unique(zsel)
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1L || is.na(n_grid) ||
      n_grid < 2) {
    stop("'n_grid' must be a single number greater than or equal to 2",
         call. = FALSE)
  }
  n_grid <- as.integer(n_grid)

  # smooth coefficients to plot; labels as in smoothcoef(x)$theta_beta
  tb0 <- x$smooth_coefficients$theta_beta
  n_coef <- NCOL(tb0)
  lab <- colnames(tb0)
  if (is.null(lab) || anyNA(lab) || any(!nzchar(lab))) {
    lab <- c("Intercept (theta)", ms$x_names)
    if (length(lab) != n_coef) lab <- paste0("coef", seq_len(n_coef))
  }

  if (is.null(which)) {
    csel <- seq_len(n_coef)
  } else {
    if (is.character(which)) {
      csel <- match(which, lab)
    } else {
      csel <- seq_len(n_coef)[which]
    }
    if (length(csel) == 0L || anyNA(csel)) {
      stop("'which' must select among the smooth coefficients: ",
           paste(lab, collapse = ", "), call. = FALSE)
    }
    csel <- unique(csel)
  }

  # values at which the z variables not on the horizontal axis are held
  fixed <- lapply(zdat, function(v) {
    if (is.numeric(v)) {
      mean(v)
    } else {
      tab <- table(v)
      v[match(names(tab)[which.max(tab)], as.character(v))]
    }
  })

  # stacked evaluation data: one block of n_grid rows per z variable plotted
  k <- length(zsel)
  n_eval <- k * n_grid
  ez <- lapply(fixed, rep, times = n_eval)

  zgrid <- vector("list", k)
  rows <- vector("list", k)
  for (j in seq_len(k)) {
    v <- zdat[[zsel[j]]]
    zgrid[[j]] <- seq(min(v), max(v), length.out = n_grid)
    rows[[j]] <- (j - 1L) * n_grid + seq_len(n_grid)
    cvals <- ez[[zsel[j]]]
    cvals[rows[[j]]] <- zgrid[[j]]
    ez[[zsel[j]]] <- cvals
  }
  ezdat <- data.frame(ez, check.names = FALSE)

  # Evaluate the smooth coefficients at the grid points via 'ezdat', with the
  # estimator of the fitted model (local-constant, or local-linear if the model
  # was estimated with spscreg = "ll"; objects made before 'spscreg' existed are
  # local-constant). The bandwidth of the fitted model is re-used, and the
  # maximum-likelihood step, whose results are not needed here, starts at the
  # estimated values. No gradients are needed.
  spscreg <- if (is.null(ms$spscreg)) "lc" else ms$spscreg

  # The model frame has the terms of the formula as columns, for instance
  # "factor(z2)": spscsfa() is called on a formula that refers to the columns
  refit <- spsc_plot_refit_args(ms$formula, mf, ezdat)

  fit_ez <- suppressWarnings(
    spscsfa(formula = refit$formula,
             data = refit$data,
             ezdat = refit$ezdat,
             sig_level = ms$sig_level,
             optim_method = x$mle_results$optim_method,
             stval = c(x$mle_results$delta_coef, x$mle_results$sigma2_v),
             spscreg = spscreg,
             bandwidth = x$smooth_coefficients$bandwidth,
             verbose = FALSE)
  )

  tb <- fit_ez$smooth_coefficients$theta_beta
  if (is.null(dim(tb)) || nrow(tb) != n_eval || ncol(tb) != n_coef) {
    stop("could not evaluate the smooth coefficients on the grid",
         call. = FALSE)
  }
  colnames(tb) <- lab

  curves <- lapply(seq_len(k), function(j) {
    d <- data.frame(zgrid[[j]], tb[rows[[j]], csel, drop = FALSE],
                    check.names = FALSE)
    names(d) <- c(zsel[j], lab[csel])
    rownames(d) <- NULL
    d
  })
  names(curves) <- zsel

  list(curves = curves, coef_names = lab[csel])
}

## Internal helper: arguments for calling spscsfa() again from the stored
## model frame.
##
## The columns of the model frame are named as the terms of the formula are
## written, for instance "factor(z2)" or "log(x1)". spscsfa() evaluates the
## terms of its formula, so given this model frame and the original formula it
## would look for z2 and x1, which are not columns of the model frame. The
## columns are therefore given syntactic names and each variable of the formula
## is replaced by the symbol of its column. The model is the same, and no term
## is evaluated again. 'ezdat', whose columns are named as the z variables of
## the model frame, is renamed in the same way.
spsc_plot_refit_args <- function(formula, mf, ezdat) {

  syntactic <- make.names(names(mf), unique = TRUE)
  map <- stats::setNames(syntactic, names(mf))

  parsed <- Formula::Formula(formula)
  n_rhs <- length(parsed)[2L]

  as_text <- function(expr) {
    paste(deparse(spsc_formula_symbols(expr, map), width.cutoff = 500L),
          collapse = " ")
  }

  lhs <- as_text(stats::formula(parsed, lhs = 1, rhs = 0)[[2L]])
  rhs <- vapply(seq_len(n_rhs), function(i) {
    as_text(stats::formula(parsed, lhs = 0, rhs = i)[[2L]])
  }, character(1))

  env <- environment(formula)
  if (is.null(env)) {
    env <- globalenv()
  }

  names(ezdat) <- unname(map[names(ezdat)])
  names(mf) <- syntactic

  list(
    formula = stats::as.formula(
      paste(lhs, "~", paste(rhs, collapse = " | ")),
      env = env
    ),
    data = mf,
    ezdat = ezdat
  )
}

## Internal helper: replaces each variable of a formula expression by the
## symbol of its column in the model frame. 'map' gives the symbol for the label
## of each variable. The operators of the formula (+, -, *, :, /, ^, %in% and
## parentheses) and numbers, as in x - 1, are kept.
spsc_formula_symbols <- function(expr, map) {

  if (!is.language(expr)) {
    return(expr)
  }

  # The label is the name that model.frame() gives the column of the variable
  label <- paste(deparse(expr, width.cutoff = 500L,
                         backtick = !is.symbol(expr) && is.language(expr)),
                 collapse = " ")

  if (label %in% names(map)) {
    return(as.name(map[[label]]))
  }

  if (is.call(expr)) {
    for (i in seq_along(expr)[-1L]) {
      expr[[i]] <- spsc_formula_symbols(expr[[i]], map)
    }
  }

  expr
}
