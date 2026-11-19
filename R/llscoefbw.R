#' Least-squares cross-validation bandwidths for the local-linear SPSC estimator
#'
#' @description
#' Selects the bandwidths of the local-linear semiparametric smooth-coefficient
#' (SPSC) estimator by least-squares cross-validation (LSCV). The result can be
#' passed to \code{\link{llscoef}} through its \code{bws} argument
#' or to \code{\link{spscsfa}} through its \code{bandwidth} argument.
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
#' \eqn{\hat\gamma(z)} estimate the smooth coefficients. The LSCV criterion is the
#' mean of the squared leave-one-out residuals.
#'
#' Variables of the second part of \code{formula} are classified by their
#' class: numeric variables are continuous, factors are unordered categorical
#' and ordered factors are ordered categorical.
#'
#' \describe{
#'   \item{Starting values}{Rule-of-thumb bandwidths:
#'     \eqn{1.059224\, sd(z)\, n^{-1/(4+p)}} for a continuous variable and
#'     \eqn{1.059224\, n^{-2/(4+p)}} for a categorical variable, where \eqn{n}
#'     is the number of observations and \eqn{p} the number of continuous
#'     variables in \eqn{z}. A categorical starting value outside the admissible
#'     range (see below) is replaced by half of the upper limit.}
#'   \item{Admissible bandwidths}{The criterion takes a
#'     very large value for inadmissible bandwidths. Continuous bandwidths must
#'     be positive, unordered categorical bandwidths must lie in
#'     \eqn{[0, (c-1)/c]} with \eqn{c} the number of levels, and ordered
#'     categorical bandwidths in \eqn{[0, 1]}. Singular local systems are
#'     ridged.}
#'   \item{Search}{Nelder-Mead, that is \code{\link[stats]{optim}} with
#'     \code{method = "Nelder-Mead"} and its default settings, started from the
#'     starting values.}
#'   \item{Scan and further searches}{The LSCV criterion is not a smooth
#'     function of the bandwidths: with few observations it can have local
#'     minima, for example at small bandwidths, where a Nelder-Mead search that
#'     starts nearby stops. Unless \code{max.search = 1}, the criterion is
#'     therefore also evaluated on a coarse scan: the starting values with
#'     every continuous bandwidth multiplied by each value in \code{grid_scan}
#'     (categorical bandwidths are kept at their starting values). By default,
#'     \code{grid_scan = 2^c(-2:5)}, which takes 8 evaluations. The scan moves
#'     all the continuous bandwidths together, and is skipped if no variable in
#'     the \code{z} part is continuous.}
#'  \item{Progress}{Progress reporting during bandwidth selection is provided by the 
#'     \pkg{progress} package.}
#' }
#'
#' @param formula A formula \code{y ~ x | z} (or \code{y ~ x | z | q}), where
#'   \code{y} is the dependent variable, \code{x} are the parametric regressors
#'   (the constant is added automatically) and \code{z} are the variables that
#'   the smooth coefficients depend on.
#' @param bws Optional numeric vector of local-linear bandwidths, one for each
#'   variable in the \code{z} part of \code{formula}, in that order, or an
#'   object returned by \code{llscoefbw}. An object returned by
#'   \code{np::npscoefbw} made with \code{regtype = "ll"}, fixed bandwidths and the
#'   default kernels is also accepted, and its bandwidths are used; one made
#'   with another \code{regtype}, \code{bwtype}, \code{bwscaling} or kernel is
#'   not, and a vector taken from it is not suitable, either. With
#'   \code{LSCV = FALSE} these bandwidths are returned without a
#'   search. Otherwise they replace the rule-of-thumb starting values.
#' @param data Optional data frame containing the variables in \code{formula}.
#'   If omitted, the variables are taken from the environment of
#'   \code{formula}.
#' @param LSCV Logical; if \code{TRUE} (default) the bandwidths are
#'   selected by LSCV. If \code{FALSE}, \code{bws} is checked, named and
#'   returned.
#' @param na.action A function, or the name of a function, that decides what to do
#'   with rows that have missing values in any part of \code{formula}, as in
#'   \code{\link[stats]{model.frame}}. The default, \code{na.omit}, drops them.
#'   \code{na.exclude} also drops them, but \code{\link{llscoef}} then pads its
#'   results with \code{NA} at those rows, so that they line up with the data.
#'   \code{na.fail} gives an error. \code{na.pass} is not allowed: an error is
#'   given if missing values remain. It is kept in the returned object, and
#'   \code{llscoef} uses it.
#' @param max.search Maximum number of Nelder-Mead searches, a positive integer
#'   (default 2): the one from the starting values and up to
#'   \code{max.search - 1} more, from the most promising points of a coarse scan of
#'   the bandwidths (see \sQuote{Details}). With \code{max.search = 1} there is no
#'   scan and a single search is made from the starting values. It is ignored
#'   with \code{LSCV = FALSE}.
#' @param grid_scan A strictly increasing numeric vector of finite positive
#'   multipliers for the coarse scan (default \code{2^c(-2:5)}). Each multiplier
#'   is applied to all continuous starting bandwidths, while categorical
#'   bandwidths are kept at their starting values. It affects the search only when
#'   \code{LSCV = TRUE}, \code{max.search > 1}, and at least one continuous
#'   variable is present in the \code{z} part of \code{formula}.
#'
#' @return With \code{LSCV = FALSE}, a numeric vector of
#'   bandwidths, one for each variable in the \code{z} part of \code{formula},
#'   in that order and named after those variables.
#'
#'   With \code{LSCV = TRUE}, an object of class
#'   \code{"llsmoothcoefficientbw"}: a list with components
#'   \describe{
#'     \item{\code{bw}}{the bandwidths, a numeric vector named after the
#'       variables in the \code{z} part of \code{formula}.}
#'     \item{\code{formula}}{the \code{formula} that was supplied.}
#'     \item{\code{n_obs}}{the number of observations used.}
#'     \item{\code{n_dim}}{the number of variables in the \code{z} part.}
#'     \item{\code{lscv}}{the value of the LSCV criterion at \code{bw}.}
#'     \item{\code{n_feval}}{the number of evaluations of the criterion, those
#'       of the scan and of every search made.}
#'     \item{\code{scan_factor}}{the multiplier used to initialize the
#'       Nelder-Mead search that produced the selected bandwidths. It is
#'       \code{1} when the search from the original starting values gives the
#'       best result; otherwise it is the corresponding value in \code{grid_scan}.}
#'     \item{\code{n_conti}, \code{n_unord}, \code{n_ord}}{the numbers of
#'       continuous, unordered categorical and ordered categorical variables
#'       in the \code{z} part.}
#'     \item{\code{call}}{the call to \code{llscoefbw}, with the environment
#'       it was made from.}
#'     \item{\code{obs_omit}, \code{n_obs_omit}}{the rows of the training data
#'       that were dropped because of missing values in any part of
#'       \code{formula}, and their number. The row numbers refer to \code{data}, 
#'       or to the variables in the environment of \code{formula}.}
#'     \item{\code{na_action}}{the \code{na.action} function used.}
#'   }
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
#' @seealso \code{\link{llscoef}}, \code{\link{spscsfa}}
#'
#' @examples
#' \dontrun{
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
#' bw
#' bw$bw
#'
#' # the formula and the data are those of llscoefbw()
#' fit <- llscoef(bw, betas = TRUE)
#'
#' # a single search from the starting values, without the scan
#' llscoefbw(y ~ x | z1 + z2 + z3, data = dat, max.search = 1)
#'
#' # check and return given bandwidths without a search
#' llscoefbw(y ~ x | z1 + z2 + z3, bws = bw$bw, data = dat, LSCV = FALSE)
#' }
#' @export
llscoefbw <- function(formula, bws, data, LSCV, na.action,
                      max.search, grid_scan) {
  UseMethod("llscoefbw")
}

#' @rdname llscoefbw
#' @export
llscoefbw.default <- function(formula, bws = NULL, data = NULL,
                              LSCV = TRUE,
                              na.action = stats::na.omit,
                              max.search = 2L,
                              grid_scan = 2^c(-2:5)) {
  stop("'formula' must be a formula such as y ~ x | z or y ~ x1 + x2 | z1 + z2",
       call. = FALSE)
}

#' @rdname llscoefbw
#' @export
llscoefbw.formula <- function(formula, bws = NULL, data = NULL,
                              LSCV = TRUE,
                              na.action = stats::na.omit,
                              max.search = 2L,
                              grid_scan = 2^c(-2:5)) {

  if (!is.logical(LSCV) || length(LSCV) != 1L ||
      is.na(LSCV)) {
    stop("'LSCV' must be TRUE or FALSE", call. = FALSE)
  }

  na.action <- .llscoef_check_na_action(na.action)

  max.search <- .llscoef_check_max_search(max.search)
  grid_scan <- .llscoef_check_grid_scan(grid_scan)

  # An object of class "scbandwidth" is used through its bandwidths if they
  # are those of the local-linear estimator: regtype = "ll" (not "lc"), fixed
  # bandwidths and the default kernels
  if (inherits(bws, "scbandwidth")) {
    bws <- .llscoef_np_bw(bws, "bws")
  }

  # The call and the environment it was made from
  cl <- match.call()
  call.env <- parent.frame()
  environment(cl) <- call.env

  # ----------------------------------------------------------
  # Data: y, X = [1, x] and Z, extracted as in spscsfa()
  # ----------------------------------------------------------

  dat <- .llscoef_model_data(formula, data, na.action)
  X <- dat$X
  Y <- dat$Y
  Z <- dat$Z

  n <- nrow(Z)

  if (ncol(Z) < 1L) {
    stop("the nonparametric (z) part of 'formula' has no variables",
         call. = FALSE)
  }

  # Continuous (numeric) versus categorical (factor) variables, from Z
  Z.is.cont <- vapply(Z, is.numeric, logical(1))
  Z.is.fac <- vapply(Z, is.factor, logical(1))

  if (any(!Z.is.cont & !Z.is.fac)) {
    stop("z variables must be numeric or factors: ",
         paste(names(Z)[!Z.is.cont & !Z.is.fac], collapse = ", "),
         call. = FALSE)
  }

  Z.is.ord <- vapply(Z, is.ordered, logical(1))
  n_conti <- sum(Z.is.cont)
  n_unord <- sum(Z.is.fac & !Z.is.ord)
  n_ord <- sum(Z.is.ord)

  bw.max <- .llscoef_bw_max(Z)

  # ----------------------------------------------------------
  # Starting values: rule-of-thumb bandwidths, in the order of Z
  # ----------------------------------------------------------

  p <- sum(Z.is.cont)

  stval <- numeric(ncol(Z))

  for (j in seq_along(stval)) {
    if (Z.is.cont[j]) {
      if (!(stats::sd(Z[[j]]) > 0)) {
        stop("continuous variable '", names(Z)[j], "' is constant",
             call. = FALSE)
      }
      stval[j] <- rot_c(Z[[j]], n, p)
    } else {
      stval[j] <- rot_d(n, p)
    }
  }

  # Supplied bandwidths
  if (!is.null(bws)) {
    if (is.list(bws) && !is.null(bws$bw)) {
      bws <- bws$bw
    }
    bws <- .llscoef_check_bw(bws, Z, bw.max)
    if (!LSCV) {
      return(stats::setNames(bws, names(Z)))
    }
    stval <- bws
  } else if (!LSCV) {
    stop("'LSCV = FALSE' requires 'bws'", call. = FALSE)
  }

  # A categorical start must lie strictly inside its admissible range
  clamp <- is.finite(bw.max) & !(stval < bw.max)
  stval[clamp] <- 0.5 * bw.max[clamp]

  # ----------------------------------------------------------
  # Search
  # ----------------------------------------------------------

  muffle_optim_warning <- function(w) {
    cl <- conditionCall(w)
    from.optim <- (is.call(cl) && identical(cl[[1L]], quote(stats::optim))) ||
      grepl("Nelder-Mead", conditionMessage(w), fixed = TRUE)
    if (length(stval) == 1L && from.optim) {
      invokeRestart("muffleWarning")
    }
  }

  # The LSCV criterion as a function of the vector of bandwidths
  crit <- function(b) ss.ll.mixed.fast.progress(b, X = X, Y = Y, Z = Z)

  # One progress session covers the scan and every search
  fit <- .spscsfa_progress_select_bandwidth(
    .llscoef_search(crit, stval, Z.is.cont, max.search, grid_scan,
                    muffle_optim_warning)
  )

  if (!is.finite(fit$value) || fit$value >= sqrt(.Machine$double.xmax)) {
    stop("bandwidth search failed: no admissible bandwidths with a finite ",
         "cross-validation value were found", call. = FALSE)
  }

  if (fit$convergence != 0L) {
    warning("optim did not converge (code ", fit$convergence,
            "); the bandwidths may not be optimal", call. = FALSE)
  }

  structure(
    list(
      bw = stats::setNames(fit$par, names(Z)),
      formula = formula,
      n_obs = n,
      n_dim = ncol(Z),
      lscv = fit$value,
      n_feval = as.integer(fit$n_feval),
      scan_factor = fit$scan_factor,
      n_conti = n_conti,
      n_unord = n_unord,
      n_ord = n_ord,
      call = cl,
      obs_omit = dat$omit,
      n_obs_omit = length(dat$omit),
      na_action = na.action
    ),
    class = "llsmoothcoefficientbw"
  )
}

#' @rdname llscoefbw
#' @param x An object of class \code{"llsmoothcoefficientbw"}.
#' @param ... Further arguments, ignored.
#' @export
print.llsmoothcoefficientbw <- function(x, ...) {
  cat("\nLocal-Linear Smooth Coefficient Model",
      "\nRegression Data (", x$n_obs, " observations, ", x$n_dim,
      " variable(s)):\n\n", sep = "")

  print(matrix(x$bw, nrow = 1L,
               dimnames = list("Bandwidth(s):", names(x$bw))))

  cat("\nFormula: ", paste(deparse(x$formula), collapse = "\n"), sep = "")
  cat("\nObjective Function Value: ", format(x$lscv), sep = "")
  cat("\nNumber of Function Evaluations: ", format(x$n_feval), sep = "")

  cat("\n")
  if (x$n_conti > 0L) {
    cat("\nNo. Continuous Explanatory Vars.: ", x$n_conti, sep = "")
  }
  if (x$n_unord > 0L) {
    cat("\nNo. Unordered Categorical Explanatory Vars.: ", x$n_unord, sep = "")
  }
  if (x$n_ord > 0L) {
    cat("\nNo. Ordered Categorical Explanatory Vars.: ", x$n_ord, sep = "")
  }

  cat("\n\n")
  invisible(x)
}

# Rule-of-thumb starting values:  x = a continuous variable, n = number of
# observations, p = number of continuous variables.
rot_c <- function(x, n, p) 1.059224 * stats::sd(x) * (n^(-1 / (4 + p)))
rot_d <- function(n, p) 1.059224 * (n^(-2 / (4 + p)))

# 'na.action' as in model.frame(): a function or the name of a function. NULL
# gives the default, na.omit.
.llscoef_check_na_action <- function(na.action) {
  if (is.null(na.action)) {
    return(stats::na.omit)
  }

  if (!is.function(na.action) &&
      !(is.character(na.action) && length(na.action) == 1L &&
        !is.na(na.action))) {
    stop("'na.action' must be a function or the name of a function, such as ",
         "na.omit or na.exclude", call. = FALSE)
  }

  na.action
}

# y, X = [1, x] and Z as spscsfa() extracts them from 'y ~ x | z' or
# 'y ~ x | z | q' (the third part is not used here).
.llscoef_model_data <- function(formula, data, na.action = stats::na.omit) {

  parsed_formula <- tryCatch(
    Formula::Formula(formula),
    error = function(e) {
      stop("Could not parse 'formula' with Formula::Formula(): ", e$message,
           call. = FALSE)
    }
  )

  f_length <- length(parsed_formula)

  if (f_length[1L] != 1L) {
    stop("Invalid formula: exactly one response variable is required on the ",
         "LHS. Use 'y ~ x | z' or 'y ~ x | z | q'.", call. = FALSE)
  }

  if (f_length[2L] < 2L || f_length[2L] > 3L) {
    stop("Invalid formula: expected 2 or 3 RHS parts separated by '|'. ",
         "Use 'y ~ x | z' or 'y ~ x | z | q'.", call. = FALSE)
  }

  x_formula <- stats::formula(parsed_formula, lhs = 0, rhs = 1)

  if (!is.null(data) && !is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Formula's model.frame() reads '|' as a separator of formula parts and
  # deals with observations with missing values in any part from 'na.action'.
  mf <- if (is.null(data)) {
    stats::model.frame(parsed_formula, na.action = na.action)
  } else {
    stats::model.frame(parsed_formula, data = data, na.action = na.action)
  }

  y <- stats::model.response(mf)

  if (!is.numeric(y)) {
    stop("Response variable (y) must be numeric", call. = FALSE)
  }

  x <- stats::model.matrix(stats::terms(x_formula), data = mf)

  if (ncol(x) > 0L && colnames(x)[1L] == "(Intercept)") {
    x <- x[, -1L, drop = FALSE]
  }

  Z <- Formula::model.part(parsed_formula, data = mf, rhs = 2)
  attr(Z, "terms") <- NULL

  # na.pass leaves missing values in the data, which the estimation cannot use
  if (anyNA(y) || anyNA(x) || any(vapply(Z, anyNA, logical(1)))) {
    stop("missing values remain after 'na.action': use na.omit, na.exclude ",
         "or na.fail", call. = FALSE)
  }

  # Rows dropped because of missing values, as row numbers in 'data', and the
  # 'na.action' object of the model frame (class "omit" or "exclude"), which
  # napredict() and naresid() use to pad results
  na.obj <- attr(mf, "na.action")
  omit <- as.vector(na.obj)
  if (is.null(omit)) {
    omit <- integer(0)
  }

  list(X = cbind(rep(1, length(y)), x), Y = y, Z = Z, omit = omit,
       na.action = na.obj)
}

# Upper limit of the bandwidth of each column of Z, the same rule as in
# ss.ll.mixed.fast.progress() for npksum()'s default kernels: Inf (continuous),
# (c - 1) / c (unordered, Aitchison-Aitken, c = number of levels), 1 (ordered).
.llscoef_bw_max <- function(Z) {
  Z.is.ord <- vapply(Z, is.ordered, logical(1))
  Z.is.uno <- vapply(Z, is.factor, logical(1)) & !Z.is.ord
  bw.max <- rep(Inf, ncol(Z))
  bw.max[Z.is.ord] <- 1
  if (any(Z.is.uno)) {
    n.lev <- vapply(Z[Z.is.uno], nlevels, integer(1))
    bw.max[Z.is.uno] <- (n.lev - 1) / n.lev
  }
  bw.max
}

# Supplied bandwidths: right length and admissible
.llscoef_check_bw <- function(bws, Z, bw.max) {
  bws <- as.numeric(bws)

  if (length(bws) != ncol(Z)) {
    stop("length of 'bws' (", length(bws), ") does not match the number of ",
         "z variables (", ncol(Z), ")", call. = FALSE)
  }

  ok <- vapply(
    seq_along(bws),
    function(j) {
      is.finite(bws[j]) && (
        if (is.finite(bw.max[j])) {
          bws[j] >= 0 && bws[j] <= bw.max[j]
        } else {
          bws[j] > 0
        }
      )
    },
    logical(1)
  )

  if (!all(ok)) {
    stop("inadmissible bandwidth for: ",
         paste(names(Z)[!ok], collapse = ", "),
         " (continuous > 0; unordered in [0, (c - 1) / c]; ordered in [0, 1])",
         call. = FALSE)
  }

  bws
}

# 'max.search': the maximum number of Nelder-Mead searches, a positive integer
.llscoef_check_max_search <- function(max.search) {
  if (!is.numeric(max.search) || length(max.search) != 1L || is.na(max.search) ||
      !is.finite(max.search) || max.search < 1 || max.search != floor(max.search) ||
      max.search > .Machine$integer.max) {
    stop("'max.search' must be a positive integer", call. = FALSE)
  }

  as.integer(max.search)
}

# 'grid_scan': strictly increasing, finite, positive multipliers for the
# continuous starting bandwidths in the coarse scan.
.llscoef_check_grid_scan <- function(grid_scan) {
  if (!is.numeric(grid_scan) || !is.null(dim(grid_scan)) ||
      length(grid_scan) < 1L || any(!is.finite(grid_scan)) ||
      any(grid_scan <= 0) ||
      (length(grid_scan) > 1L && any(diff(grid_scan) <= 0))) {
    stop("'grid_scan' must be a strictly increasing numeric vector of finite ",
         "positive values", call. = FALSE)
  }

  as.numeric(grid_scan)
}

# Nelder-Mead search(es) for the LSCV bandwidths.
#   crit          the criterion, a function of the vector of bandwidths
#   stval         starting values
#   is.cont       which of the bandwidths are those of continuous variables
#   max.search    maximum number of searches (a positive integer)
#   grid_scan     strictly increasing scan multipliers
#   warn.handler  calling handler for the warnings of optim()
# Returns a list: par, value and convergence of the search that is kept, and
# n_feval, the evaluations of the criterion in the scan and in all searches;
# scan_factor, the multiplier used to initialize the search with the best fit.
.llscoef_search <- function(crit, stval, is.cont, max.search, grid_scan,
                            warn.handler) {

  big <- sqrt(.Machine$double.xmax)

  run <- function(start) {
    withCallingHandlers(
      stats::optim(start, crit, method = "Nelder-Mead"),
      warning = warn.handler
    )
  }

  best <- run(stval)
  best_scan_factor <- 1
  n_feval <- as.integer(best$counts[["function"]])

  if (max.search > 1L && any(is.cont)) {

    starts <- lapply(grid_scan, function(f) {
      b <- stval
      b[is.cont] <- f * stval[is.cont]
      b
    })

    value <- vapply(starts, crit, numeric(1))
    n_feval <- n_feval + length(value)

    # Inadmissible points and failed fits take the penalty value
    value[!is.finite(value) | value >= big] <- Inf

    # Local minima of the scan
    n.scan <- length(value)
    below.left <- value <= c(Inf, value[-n.scan])
    below.right <- value <= c(value[-1L], Inf)
    cand <- which(is.finite(value) & below.left & below.right)
    cand <- cand[order(value[cand])]

    n.extra <- 0L

    for (k in cand) {
      if (n.extra >= max.search - 1L) {
        break
      }

      # Only a point below what has been found can lead to something better;
      # this also leaves out the scan point that is the starting value
      if (!(value[k] < best$value)) {
        next
      }

      fit.k <- run(starts[[k]])
      n_feval <- n_feval + as.integer(fit.k$counts[["function"]])
      n.extra <- n.extra + 1L

      if (fit.k$value < best$value) {
        best <- fit.k
        best_scan_factor <- grid_scan[k]
      }
    }
  }

  list(par = best$par, value = best$value, convergence = best$convergence,
       n_feval = n_feval, scan_factor = best_scan_factor)
}

.spscsfa_progress <- local({
  env <- new.env(parent = emptyenv())
  env$depth <- 0L
  env$active <- FALSE
  env$visible <- FALSE
  env$rendered <- FALSE
  env$started <- 0
  env$last_emit <- 0
  env$grace <- 1
  env$interval <- 2
  env$ridging <- FALSE
  env$bar <- NULL
  env$aborted <- FALSE
  env
})

.spscsfa_progress_now <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

.spscsfa_progress_option <- function(name, default) {
  val <- suppressWarnings(as.numeric(getOption(name, default))[1L])
  if (!is.finite(val) || val < 0) default else val
}

.spscsfa_progress_enabled <- function() {
  interactive() &&
    isTRUE(getOption("spscsfa.messages", getOption("np.messages", TRUE))) &&
    isTRUE(getOption("progress_enabled", TRUE))
}

.spscsfa_progress_width <- function() {
  width <- suppressWarnings(as.integer(getOption("width", 80L))[1L])
  if (is.na(width)) {
    width <- 80L
  }
  max(20L, width - 4L)
}

.spscsfa_progress_begin <- function() {
  st <- .spscsfa_progress
  st$depth <- st$depth + 1L
  if (st$depth > 1L) {
    return(invisible(NULL))
  }

  now <- .spscsfa_progress_now()
  st$grace <- .spscsfa_progress_option("spscsfa.progress.start.grace.sec", 1)
  st$interval <- .spscsfa_progress_option("spscsfa.progress.interval.sec", 2)
  st$visible <- .spscsfa_progress_enabled()
  st$rendered <- FALSE
  st$started <- now
  st$last_emit <- now - st$interval
  st$ridging <- FALSE
  st$bar <- NULL
  st$aborted <- FALSE
  st$active <- TRUE

  if (st$visible) {
    st$bar <- tryCatch(
      progress::progress_bar$new(
        format = paste0(
          "{spscsfa} Bandwidth selection (elapsed: :sps_elapsed:f.value:note)"
        ),
        total = NA,
        width = .spscsfa_progress_width(),
        clear = TRUE,
        show_after = st$grace,
        message_class = "packageStartupMessage"
      ),
      error = function(e) NULL
    )

    if (is.null(st$bar)) {
      st$visible <- FALSE
    } else {
      tryCatch(
        {
          st$bar$tick(0, tokens = list(sps_elapsed = "0.0s", f.value = "", note = ""))
          st$rendered <- st$grace <= 0
        },
        error = function(e) {
          st$bar <- NULL
          st$visible <- FALSE
          st$rendered <- FALSE
        }
      )
    }
  }

  invisible(NULL)
}

.spscsfa_progress_end <- function() {
  st <- .spscsfa_progress
  if (st$depth > 0L) {
    st$depth <- st$depth - 1L
  }
  if (st$depth > 0L) {
    return(invisible(NULL))
  }

  if (!isTRUE(st$aborted)) {
    .spscsfa_progress_terminate()
  }
  st$bar <- NULL
  st$active <- FALSE
  st$visible <- FALSE
  st$rendered <- FALSE
  st$ridging <- FALSE
  st$aborted <- FALSE
  invisible(NULL)
}

.spscsfa_progress_terminate <- function() {
  st <- .spscsfa_progress
  if (!is.null(st$bar)) {
    try(st$bar$terminate(), silent = TRUE)
  }
  invisible(NULL)
}

# Search stopped by an error or interrupt: keep the last line and start a new one.
.spscsfa_progress_abort <- function() {
  st <- .spscsfa_progress
  if (!isTRUE(st$active) || isTRUE(st$aborted)) {
    return(invisible(NULL))
  }

  st$aborted <- TRUE
  if (isTRUE(st$rendered)) {
    try({
      cat("\n", file = stderr())
      flush(stderr())
    }, silent = TRUE)
  }
  st$rendered <- FALSE
  invisible(NULL)
}

# Called once per criterion evaluation. Ridging is reported at the next refresh.
.spscsfa_progress_tick <- function(ridging = FALSE, f.value = NA_real_) {
  st <- .spscsfa_progress
  if (!st$active || !st$visible || is.null(st$bar) || isTRUE(st$aborted)) {
    return(invisible(NULL))
  }

  tryCatch({
    if (isTRUE(ridging)) {
      st$ridging <- TRUE
    }

    now <- .spscsfa_progress_now()
    if ((now - st$started) < st$grace ||
        (now - st$last_emit) < st$interval) {
      return(invisible(NULL))
    }

    elapsed <- formatC(max(0, now - st$started), digits = 1L, format = "f")
    f.value_text <- if (is.finite(f.value)) {
      paste0(", f.value: ", formatC(f.value, digits = 9L, format = "fg"))
    } else {
      ""
    }

    note <- ""
    if (st$ridging) {
      line <- paste0(
        "{spscsfa} Bandwidth selection (elapsed: ", elapsed, "s", f.value_text
      )
      for (candidate in c(
        ", singular or near-singular local system, ridging",
        ", ridging"
      )) {
        if (nchar(paste0(line, candidate, ")"), type = "width") <=
            .spscsfa_progress_width()) {
          note <- candidate
          break
        }
      }
    }

    st$bar$tick(
      0,
      tokens = list(sps_elapsed = paste0(elapsed, "s"), f.value = f.value_text, 
                    note = note)
    )
    st$rendered <- TRUE
    st$ridging <- FALSE
    st$last_emit <- now
  }, error = function(e) {
    st$visible <- FALSE
    if (!is.null(st$bar)) {
      try(st$bar$terminate(), silent = TRUE)
    }
    st$bar <- NULL
    st$rendered <- FALSE
  })

  invisible(NULL)
}

.spscsfa_progress_select_bandwidth <- function(expr) {
  .spscsfa_progress_begin()
  on.exit(.spscsfa_progress_end(), add = TRUE)
  withCallingHandlers(
    expr,
    error = function(e) .spscsfa_progress_abort(),
    interrupt = function(i) .spscsfa_progress_abort()
  )
}

# ==============================================================================
# LSCV criterion
# ==============================================================================

ss.ll.mixed.fast.progress <- function(bw, X, Y, Z) {

  ukertype <- "aitchisonaitken"
  okertype <- "liracine"

  # Returned for inadmissible bandwidths and for failed fits.
  max.penalty <- sqrt(.Machine$double.xmax)

  # ============================================================
  # For evaluation point i, S_j^(i) = L_i T_j, where T_j is the design
  # built from X_j and the globally centered continuous variables Z_gj,
  # and L_i is a small (nc x nc) matrix that restores centering at Z_i.
  # Hence
  #     A_i = sum_{j != i} K_ij S_j^(i) S_j^(i)' = L_i A~_i L_i'
  #     B_i = sum_{j != i} K_ij S_j^(i) Y_j      = L_i B~_i
  # where A~_i and B~_i are the corresponding cross-products formed from
  # T_j. One leave-one-out npksum() call supplies them for every i.
  # ============================================================

  # ----------------------------------------------------------
  # Admissible bandwidths:
  #   continuous        finite and > 0
  #   unordered factor  0 <= bw <= (c - 1) / c   aitchisonaitken, c = nlevels
  #                     0 <= bw <= 1             liracine
  #   ordered factor    0 <= bw <= 1             (all ordered kernels)
  # ----------------------------------------------------------

  Z <- as.data.frame(Z)

  if (length(bw) != ncol(Z)) {
    stop("length of bandwidth vector does not match number of columns of 'Z'")
  }

  Z.is.ord <- vapply(Z, is.ordered, logical(1))
  Z.is.uno <- vapply(Z, is.factor, logical(1)) & !Z.is.ord
  Z.is.cat <- Z.is.ord | Z.is.uno

  if (!all(is.finite(bw)) || any(bw[!Z.is.cat] <= 0)) {
    return(max.penalty)
  }

  if (any(Z.is.cat)) {
    bw.max <- rep(1, ncol(Z))
    if (ukertype == "aitchisonaitken") {
      n.lev <- vapply(Z[Z.is.uno], nlevels, integer(1))
      bw.max[Z.is.uno] <- (n.lev - 1) / n.lev
    }
    if (any(bw[Z.is.cat] < 0 | bw[Z.is.cat] > bw.max[Z.is.cat])) {
      return(max.penalty)
    }
  }

  # Count only evaluated points. 'ridged' tells the progress
  # line whether a singular local system had to be ridged in this evaluation.
  ridged <- FALSE
  f.value <- NA_real_
  on.exit(.spscsfa_progress_tick(ridging = ridged, f.value = f.value), add = TRUE)
  
  X <- as.matrix(X)
  y <- drop(as.matrix(Y))

  n.train <- nrow(Z)
  k <- ncol(X)

  Z.col.numeric <- vapply(Z, is.numeric, logical(1))
  Z_c <- as.matrix(Z[, Z.col.numeric, drop = FALSE])
  p <- ncol(Z_c)

  # ----------------------------------------------------------
  # Globally centered local-linear design
  #
  # T = [X, X[,1]*Z_g, ..., X[,k]*Z_g]   (X-major ordering)
  #
  # Z_g is Z_c centered once at its global mean, giving one common design
  # matrix for all evaluation points. The transformation L_i below shifts
  # the local-linear terms to the evaluation point Z_i and therefore
  # recovers the exact observation-specific centering used by the
  # pointwise implementation.
  # ----------------------------------------------------------

  if (p > 0) {
    Z_g <- scale(Z_c, center = TRUE, scale = FALSE)
    XZ.col <- rep(seq_len(k), each = p)
    Z.col <- rep(seq_len(p), times = k)
    Tm <- cbind(
      X,
      X[, XZ.col, drop = FALSE] * Z_g[, Z.col, drop = FALSE]
    )
  } else {
    Tm <- X
  }

  nc <- ncol(Tm)

  # ----------------------------------------------------------
  # ONE kernel-weighted cross-product call.
  # Leave-one-out is performed by npksum().
  #
  # Using the same matrix W = [T, y] as tydat and weights gives, for
  # each i, the symmetric (nc+1) x (nc+1) matrix
  #     sum_{j != i} K_ij W_j W_j'  =  [ A~_i   B~_i ]
  #                                    [ B~_i'    *  ]
  # ----------------------------------------------------------

  W <- unname(cbind(Tm, y))

  ks <- np::npksum(
    bws = bw,
    txdat = Z,
    tydat = W,
    weights = W,
    leave.one.out = TRUE,
    ukertype = ukertype,
    okertype = okertype
  )$ksum

  # Non-finite kernel sums
  if (!all(is.finite(ks))) {
    return(max.penalty)
  }

  # ----------------------------------------------------------
  # Observation-specific centering:  S_j^(i) = L_i T_j
  #
  # L_i = [ I_k                 0     ]
  #       [ -(I_k (x) g_i)    I_{kp}  ],      g_i = Z_g[i, ]
  #
  # Only the lower-left block changes with evaluation point i.
  # ----------------------------------------------------------

  idx <- seq_len(nc)
  L <- diag(nc)

  if (p > 0) {
    L.pos <- cbind(k + seq_len(k * p), XZ.col)
    L.val <- -Z_g[, Z.col, drop = FALSE]
  }

  Y.hat <- numeric(n.train)

  for (i in seq_len(n.train)) {

    A <- ks[idx, idx, i]
    B <- ks[idx, nc + 1L, i]

    if (p > 0) {
      L[L.pos] <- L.val[i, ]
      A <- tcrossprod(L %*% A, L)
      B <- as.vector(L %*% B)
    }

    # -------------------------
    # No singularity problems
    # -------------------------

    A.inv <- tryCatch(
      suppressWarnings(as.matrix(solve(A))),
      error = function(e) NULL
    )

    if (is.null(A.inv)) {

      # --------------------------------------------------------
      # Ordinary diagonal ridging, matching sscm.ll.mixed() in llscoef.R:
      # add ridge * I to A, leave B unchanged, and increase ridge by
      # 1/n.train until the system can be solved.  No warning is issued
      # during bandwidth selection; the progress line notes ridging.
      # --------------------------------------------------------

      ridged <- TRUE
      ridging <- .ordinary_diagonal_ridging(
        A = A,
        n.train = n.train,
        n.diag = nc,
        suppress.warnings = TRUE
      )
      A.inv <- ridging$A.inv
      if (is.null(A.inv)) {
        return(max.penalty)
      }
    }

    theta <- as.vector(A.inv %*% B)

    Y.hat[i] <- sum(X[i, ] * theta[seq_len(k)])
  }

  if (!all(is.finite(Y.hat))) {
    return(max.penalty)
  }

  cv <- sum((y - Y.hat)^2) / length(y)

  if (!is.finite(cv)) {
    return(max.penalty)
  }
  f.value <- cv
  return(cv)
}

# Internal helper for ordinary diagonal ridging.
.ordinary_diagonal_ridging <- function(A, n.train, n.diag,
                                       suppress.warnings = FALSE,
                                       eval = NULL) {
  ridge_step <- 1 / n.train
  ridge <- 0
  A.inv <- NULL

  while (is.null(A.inv)) {
    ridge <- ridge + ridge_step

    if (ridge > 1e6) {
      break
    }

    A.inv <- tryCatch({
      if (suppress.warnings) {
        suppressWarnings(solve(A + diag(ridge, n.diag)))
      } else {
        solve(A + diag(ridge, n.diag))
      }
    }, error = function(e) NULL)
  }

  if (!is.null(A.inv) && !is.null(eval)) {
    warning(
      "Ridged at the evaluation point ", eval,
      ", ridge = ", signif(ridge, 6),
      call. = FALSE
    )
  }

  list(A.inv = A.inv, ridge = ridge)
}
