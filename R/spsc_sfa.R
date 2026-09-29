#' Semiparametric Smooth-Coefficient Stochastic Frontier Analysis
#'
#' Estimates a semiparametric smooth-coefficient stochastic frontier model
#' \eqn{y_i = \alpha(Z_i) + X_i' \beta(Z_i) + v_i - u(Z_i)}
#' following Sun and Kumbhakar (2013). The model allows the coefficients of the
#' parametric part to vary smoothly with a set of nonparametric variables.
#' Technical efficiency and inefficiency are computed using the Battese and
#' Coelli (1988) and Jondrow et al. (1982) methods, respectively.
#' Confidence intervals for technical efficiency are computed using the
#' approach of Horrace and Schmidt (1996).
#'
#' @details
#' The stochastic production frontier model is specified as
#' \deqn{
#' y_i = \alpha(Z_i) + X_i' \beta(Z_i) + v_i - u(Z_i),
#' }{
#' y_i = alpha(Z_i) + X_i' beta(Z_i) + v_i - u(Z_i)
#' }
#' where \eqn{y_i} can be viewed as the log of output, \eqn{X_i} can be viewed as a
#' vector of the log of inputs, and \eqn{Z_i} is a vector of covariates.
#' \eqn{\alpha(\cdot)} is the intercept, and \eqn{\beta(\cdot)} is the slope vector.
#' Both \eqn{\alpha(\cdot)} and \eqn{\beta(\cdot)} are unknown functions of \eqn{Z_i}.
#' The noise term \eqn{v_i \sim iid N(0, \sigma_v^2)}.
#' The inefficiency term is \eqn{u(Z_i) = \sigma_u(Z_i) \eta_i}, where
#' \eqn{\eta_i \sim iid N^+(0,1)}.
#' The random variables \eqn{\eta_i} and \eqn{v_i} are assumed to be independent of
#' each other and independent of \eqn{X} and \eqn{Z}.
#'
#' To guarantee positivity of \eqn{\sigma_u(Z_i)},
#' \deqn{
#' \sigma_u(Z_i) = \exp(\delta_0 + \delta_1' Z_i).
#' }{
#' sigma_u(Z_i) = exp(delta_0 + delta_1' Z_i).
#' }
#' These imply that
#' \deqn{
#' E(u_i \mid Z_i)
#' = \sigma_u(Z_i) E(\eta_i \mid Z_i)
#' = \sqrt{2/\pi}\,\sigma_u(Z_i)
#' = \sqrt{2/\pi}\,\exp(\delta_0 + \delta_1' Z_i)
#' }{
#' E(u_i | Z_i)
#' = sqrt(2/pi) * sigma_u(Z_i)
#' = sqrt(2/pi) * exp(delta_0 + delta_1' Z_i)
#' }
#' and
#' \eqn{
#' \sigma^2_u(Z_i) = \exp[2(\delta_0 + \delta_1' Z_i)].
#' }
#'
#' The estimation involves two steps. In the first step, estimate
#' \eqn{\theta(Z_i)} (\code{Intercept (theta)}), \eqn{\beta(Z_i)}, and
#' \eqn{\varepsilon_i} from
#' \deqn{
#' y_i = \theta(Z_i) + X_i'\beta(Z_i) + \varepsilon_i,
#' }{
#' y_i = theta(Z_i) + X_i' beta(Z_i) + varepsilon_i,
#' }
#' where
#' \eqn{
#' \theta(Z_i) = \alpha(Z_i)-E(u_i\mid Z_i)
#' } and 
#' \eqn{
#' \varepsilon_i = v_i-[u(Z_i)-E(u_i\mid Z_i)]
#' = E(u_i\mid Z_i) + v_i-u(Z_i).
#' }
#'
#' Then in the second step, use the estimated \eqn{\varepsilon_i} to estimate
#' \eqn{\delta_0}, \eqn{\delta_1} and \eqn{\sigma^2_v} from
#' \deqn{
#' \varepsilon_i =
#' \sqrt{2/\pi}\,\exp(\delta_0+\delta_1'Z_i)
#' + v_i
#' - \exp(\delta_0+\delta_1'Z_i)\eta_i.
#' }{
#' varepsilon_i = sqrt(2/pi) * exp(delta_0 + delta_1' Z_i) + v_i
#' - exp(delta_0 + delta_1' Z_i) * eta_i.
#' }
#'
#' Finally, compute
#' \eqn{E(u_i \mid Z_i) = \sqrt{2/\pi}\,\exp(\delta_0 + \delta_1' Z_i)}
#' from estimated \eqn{\delta_0} and \eqn{\delta_1}. The intercept in the original
#' stochastic production frontier \eqn{\alpha(Z_i)}
#' (\code{Intercept (alpha)}) is recovered as
#' \deqn{
#' \alpha(Z_i)
#' =\theta(Z_i)+E(u_i\mid Z_i).
#' }{
#' alpha(Z_i) = theta(Z_i) + E(u_i|Z_i).
#' }
#'
#' See Sun and Kumbhakar (2013, page 306) for more details.
#'
#' @param formula A formula of the form \code{y ~ x | z} or
#'   \code{y ~ x | z | q},
#'   where \code{y} is the dependent variable,
#'   \code{x} are parametric inputs,
#'   \code{z} are non-parametric variables determining smooth coefficients,
#'   and \code{q} (optional) are inefficiency determinants.
#'   When the third formula part \code{q} is omitted, the \code{z} variables also 
#'   serve as inefficiency determinants.
#'
#' @param data Optional data frame containing all variables used in \code{formula}.
#'   If supplied, variables are extracted from this data frame.
#'   If not supplied, variables are taken from the environment of \code{formula}.
#'
#' @param ezdat Optional data frame containing evaluation points for the
#'   non-parametric variables \code{z}. If provided, smooth coefficients
#'   are estimated at these points; otherwise, the training data points
#'   (from \code{z}) are used.
#'
#' @param bandwidth Optional bandwidth object (either from \code{npscoefbw} or
#'   a numeric vector) for the smooth-coefficient estimation. If \code{NULL}
#'   (default), bandwidth is selected by least-squares cross-validation (LSCV).
#'
#' @param sig_level Significance level for confidence intervals of technical
#'   efficiency (default 0.05).
#'
#' @param optim_method Character string specifying the optimizer:
#'   \code{"optim"} (uses Nelder-Mead), \code{"optim_BFGS"} (uses BFGS,
#'   the default), or \code{"nlm"}.
#'
#' @param optim_control List of control parameters passed to the optimizer.
#'
#' @param stval Starting values for the MLE. Use \code{"automatic"} (default)
#'   for internally computed starting values, or supply a numeric vector 
#'   in which the order of its elements is \eqn{\delta_0}, \eqn{\delta_1} (if any), and \eqn{\sigma^2_v}.
#'   See Details.
#'
#' @param penalty To ensure positivity of estimated \eqn{\sigma^2_v}, penalty
#'   applied when \eqn{\sigma^2_v} becomes negative during a numerical search
#'   (default 1e8).
#'
#' @param verbose Logical; if \code{TRUE}, progress and intermediate results
#'   are printed.
#'
#' @return An object of class \code{"spsc_sfa"} containing smooth coefficients, MLE
#'   results, technical efficiency, JLMS inefficiency, model summary, data
#'   information, and the parsed formula.
#'   
#'   A list of class \code{"spsc_sfa"} containing:
#'   \item{smooth_coefficients}{Smooth coefficient estimates, residuals,
#'        bandwidth, and metadata.}
#'   \item{mle_results}{MLE parameter estimates, standard errors, and
#'        log-likelihood.}
#'   \item{technical_efficiency}{Technical efficiency scores (BC method)
#'        and confidence intervals.}
#'   \item{inefficiency}{Inefficiency estimates (JLMS method).}
#'   \item{model_summary}{Summary information about the model specification.}
#'   \item{formula_parsed}{Parsed formula components.}
#'   \item{data_info}{Information about data cleaning and observations used.}
#' 
#' @note
#' Results can be summarized using \code{\link{summary.spsc_sfa}}.
#' Extraction functions for smooth coefficients, technical efficiency and
#' inefficiency are available through \code{\link{smooth_coef}},
#' \code{\link{te}}, and \code{\link{inefficiency}}.
#' 
#' @author
#' Kai Sun <ksun1@shu.edu.cn>
#' 
#' @references
#' George E. Battese and Tim J. Coelli (1988).
#' Prediction of firm-level technical efficiencies with a generalized frontier
#' production function and panel data.
#' \emph{Journal of Econometrics}, 38(3), 387--399.
#'
#' William C. Horrace and Peter Schmidt (1996).
#' Confidence statements for efficiency estimates from stochastic frontier models.
#' \emph{Journal of Productivity Analysis}, 7(2), 257--282.
#'
#' James Jondrow, C.A. Knox Lovell, Ivan S. Materov, and Peter Schmidt (1982).
#' On the estimation of technical inefficiency in the stochastic frontier
#' production function model.
#' \emph{Journal of Econometrics}, 19(2--3), 233--238.
#'
#' Kai Sun and Subal C. Kumbhakar (2013).
#' Semiparametric smooth-coefficient stochastic frontier model.
#' \emph{Economics Letters}, 120(2), 305--309.
#' 
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' z <- runif(n)
#' ez <- seq(0, 1, length = 1000)
#' q <- rnorm(n)
#' y <- cos(z) + sin(z) * x + rnorm(n, 0, 0.1) -
#'      exp(0.1 + 0.1 * z) * abs(rnorm(n, 0, 1))
#' dat <- data.frame(output = y, input = x, covariate = z, ineffvar = q)
#'
#' # covariate enters the nonparametric part and is also treated as inefficiency determinant
#' result <- spsc_sfa(output ~ input | covariate, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#' result$smooth_coefficients
#' result$mle_results
#' result$technical_efficiency
#' result$inefficiency
#' \donttest{
#' # This gives the same results: the inefficiency variables default to the nonparametric part
#' result <- spsc_sfa(output ~ input | covariate | NULL, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # Use evaluation data (ezdat) for estimating the smooth coefficients
#' result <- spsc_sfa(output ~ input | covariate, data = dat,
#'                    ezdat = data.frame(covariate = ez),
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#' nrow(result$smooth_coefficients$theta_beta)
#'
#' # Three-part formula: nonparametric part differs from inefficiency determinants
#' result <- spsc_sfa(output ~ input | covariate | ineffvar, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # Customized bandwidth for the nonparametric part
#' result <- spsc_sfa(output ~ input | covariate | ineffvar, data = dat,
#'                    optim_method = "optim_BFGS", bandwidth = c(0.1), verbose = TRUE)
#'
#' # Restrict no inefficiency determinants (intercept only)
#' result <- spsc_sfa(output ~ input | covariate | 1, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # More variables
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' z1 <- runif(n)
#' z2 <- runif(n)
#' y <- cos(z1 + z2) + sin(z1 - z2) * x1 + sin(z1 + z2) * x2 +
#'      rnorm(n, 0, 0.1) - exp(0.1 + 0.1 * z1 + 0.1 * z2) * abs(rnorm(n, 0, 1))
#' dat <- data.frame(y, x1, x2, z1, z2)
#'
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2, data = dat)
#'
#' # z3 is an unordered categorical variable
#' # spsc_sfa() detects if z3 is a factor
#' dat$z3 <- factor(rep(sample(c("A", "B", "C", "D", "E")), 20))
#'
#' # z1, z2, z3 enter the nonparametric part and are also inefficiency determinants
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2 + z3, data = dat)
#'
#' # z3 is preserved as an unordered categorical variable in the nonparametric part
#' result$smooth_coefficients$bandwidth
#' 
#' # z3 is converted to dummy variables as inefficiency determinants
#' result$mle_results$delta_table
#'
#' # Same as above (inefficiency defaults to nonparametric part)
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2 + z3 | NULL, data = dat)
#'
#' # Three-part formula: use z3 only as the inefficiency determinant
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2 + z3 | z3, data = dat)
#'
#' # Use external variables as inefficiency determinants
#' dat$q1 <- rnorm(n)
#' dat$q2 <- rnorm(n)
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2 + z3 | q1 + q2, data = dat)
#'
#' # No inefficiency determinants (intercept only)
#' result <- spsc_sfa(y ~ x1 + x2 | z1 + z2 + z3 | 1, data = dat, verbose = TRUE)
#'
#' # Application with real data (ricephil)
#' \dontrun{
#' 
#' data(ricephil)
#' 
#' ricephil$FARMERCODE <- as.factor(ricephil$FARMERCODE)
#' 
#' # this example shows categorical variable handling
#' result <- spsc_sfa(PROD ~ AREA + LABOR + NPK | EDYRS + FARMERCODE | FARMERCODE,
#'                    data = ricephil, optim_method = "optim_BFGS", verbose = TRUE)
#' 
#' # FARMERCODE is preserved as an unordered categorical variable in the nonparametric part
#' result$smooth_coefficients$bandwidth
#' 
#' # FARMERCODE is converted to dummy variables as inefficiency determinants
#' result$mle_results$delta_table
#' }
#' }
#' 
#' @import Formula
#' @import np
#' @import stats
#' @export
spsc_sfa <- function(formula, 
                     data = NULL,
                     ezdat = NULL,
                     sig_level = 0.05,
                     optim_method = "optim_BFGS",
                     optim_control = list(maxit = 5000),
                     stval = "automatic",
                     bandwidth = NULL,
                     penalty = 1e8, 
                     verbose = TRUE) {
  
  # ========== INPUT VALIDATION ==========
  
  if (verbose) cat("\n========== INPUT VALIDATION ==========\n")
  
  # ========== PARSE FORMULA ==========
  
  if (verbose) cat("Parsing formula...\n")
  
  
  # Parse formula using Formula package.
  #
  # Formula::lhs() and Formula::rhs() are not exported by the Formula
  # package. Use the exported formula() method to extract each component.
  #
  # Supported syntax:
  #   y ~ x | z
  #   y ~ x | z | q
  #
  # where x = parametric regressors, z = nonparametric variables, and
  # w = inefficiency determinants.
  parsed_formula <- tryCatch(
    Formula::Formula(formula),
    error = function(e) {
      stop("Could not parse 'formula' with Formula::Formula(): ", e$message)
    }
  )

  # Formula length gives the number of LHS and RHS parts.
  # For example, length(Formula(y ~ x | z | q)) == c(1, 3).
  f_length <- length(parsed_formula)

  if (f_length[1L] != 1L) {
    stop(
      "Invalid formula: exactly one response variable is required on the LHS. ",
      "Use 'y ~ x | z' or 'y ~ x | z | q'."
    )
  }

  n_rhs <- f_length[2L]
  if (n_rhs < 2L || n_rhs > 3L) {
    stop(
      "Invalid formula: expected 2 or 3 RHS parts separated by '|'. ",
      "Use 'y ~ x | z' or 'y ~ x | z | q'."
    )
  }

  # Extract standard R formula objects for each Formula component.
  # lhs = 0 removes the LHS; rhs = 0 removes the RHS.
  y_part_formula <- tryCatch(
    stats::formula(parsed_formula, lhs = 1, rhs = 0),
    error = function(e) {
      stop("Cannot extract response component from formula: ", e$message)
    }
  )

  x_formula <- tryCatch(
    stats::formula(parsed_formula, lhs = 0, rhs = 1),
    error = function(e) {
      stop(
        "Cannot extract first RHS component (parametric variables x): ",
        e$message
      )
    }
  )

  z_formula <- tryCatch(
    stats::formula(parsed_formula, lhs = 0, rhs = 2),
    error = function(e) {
      stop(
        "Cannot extract second RHS component (nonparametric variables z): ",
        e$message
      )
    }
  )

  inefficiency_formula <- NULL
  if (n_rhs == 3L) {
    inefficiency_formula <- tryCatch(
      stats::formula(parsed_formula, lhs = 0, rhs = 3),
      error = function(e) {
        stop(
          "Cannot extract third RHS component (inefficiency determinants): ",
          e$message
        )
      }
    )
  }

  # Keep the response formula under the name y_formula for compatibility
  # with the rest of the function and the returned formula_parsed object.
  y_formula <- y_part_formula

  if (verbose) {
    cat("Formula parsed successfully\n")
    cat("  Response (y): ", deparse(y_formula[[2]]), "\n", sep = "")
    cat("  Parametric (x): ", deparse(x_formula), "\n", sep = "")
    cat("  Nonparametric (z): ", deparse(z_formula), "\n", sep = "")
    if (!is.null(inefficiency_formula)) {
      cat("  Inefficiency determinants: ", deparse(inefficiency_formula), "\n", sep = "")
    }
  }
  
  # ========== DATA PREPARATION ==========
  
  if (verbose) cat("\n========== DATA PREPARATION AND CLEANING ==========\n")
  
  # Coerce data to data frame if supplied
  data_supplied <- !is.null(data)
  if (data_supplied) {
    if (!is.data.frame(data)) {
      tryCatch({
        data <- as.data.frame(data)
      }, error = function(e) {
        stop("Could not coerce 'data' to data.frame: ", e$message)
      })
    }
    if (verbose) cat("Data frame supplied (", nrow(data), " rows, ", ncol(data), " columns)\n", sep = "")
  }
  
  # Use Formula's model.frame method so that '|' is interpreted as a
  # formula-part separator rather than as base-R logical OR.
  # This also removes observations with missing values in any formula part.
  if (data_supplied) {
    clean_data <- tryCatch({
      stats::model.frame(parsed_formula, data = data, na.action = na.omit)
    }, error = function(e) {
      stop(
        "Error in data cleaning with Formula::model.frame(): ",
        e$message
      )
    })
  } else {
    clean_data <- tryCatch({
      stats::model.frame(parsed_formula, na.action = na.omit)
    }, error = function(e) {
      stop(
        "Error in data cleaning with Formula::model.frame(): ",
        e$message
      )
    })
  }

  n_total <- if (data_supplied) nrow(data) else NA
  n_used <- nrow(clean_data)
  n_removed <- n_total - n_used
  
  if (verbose) {
    cat("After removing observations with missing values:\n")
    if (!is.na(n_total)) {
      cat("  Total observations: ", n_total, "\n", sep = "")
      cat("  Observations removed: ", n_removed, "\n", sep = "")
    }
    cat("  Observations used: ", n_used, "\n", sep = "")
  }
  
  # ========== EXTRACT DATA MATRICES ==========
  
  if (verbose) cat("\nExtracting data matrices...\n")
  
  # Extract y (response variable) from the explicit one-response formula.
  y_extracted <- tryCatch({
    y_frame <- stats::model.frame(y_part_formula, data = clean_data)
    stats::model.response(y_frame)
  }, error = function(e) {
    stop("Could not extract response variable from formula: ", e$message)
  })

  if (!is.numeric(y_extracted)) {
    stop("Response variable (y) must be numeric")
  }
  
  if (verbose) cat("  y: ", length(y_extracted), " observations\n", sep = "")
  
  # Extract x (parametric variables) - as matrix
  x_terms <- attr(stats::model.frame(x_formula, data = clean_data), "terms")
  x_extracted <- stats::model.matrix(x_terms, data = clean_data)
  
  # Remove intercept if present (we'll add it later if needed)
  if (colnames(x_extracted)[1] == "(Intercept)") {
    x_extracted <- x_extracted[, -1, drop = FALSE]
  }
  
  if (verbose) cat("  x: ", nrow(x_extracted), " observations, ", ncol(x_extracted), " variables\n", sep = "")
  
  # Extract z (nonparametric variables) - as data frame
  z_extracted <- stats::model.frame(z_formula, data = clean_data)
  
  if (verbose) cat("  z: ", nrow(z_extracted), " observations, ", ncol(z_extracted), " variables\n", sep = "")
  
  # Extract inefficiency variables from formula's third part if present,
  # otherwise default to z variables.
  inefficiency_vars_data <- NULL
  if (!is.null(inefficiency_formula)) {
    inefficiency_vars_data <- stats::model.frame(inefficiency_formula, data = clean_data)
    if (verbose) cat("  inefficiency variables: ", nrow(inefficiency_vars_data), " observations, ", 
                     ncol(inefficiency_vars_data), " variables\n", sep = "")
  } else {
    # Default: use z variables for inefficiency (as is, with factors)
    inefficiency_vars_data <- z_extracted
    if (verbose) cat("  inefficiency variables (default from z): ", nrow(inefficiency_vars_data), 
                     " observations, ", ncol(inefficiency_vars_data), " variables\n", sep = "")
  }
  
  # ========== EVALUATION DATA (ezdat) ==========
  
  if (verbose) cat("\nPreparing evaluation data...\n")
  
  ezdat_extracted <- NULL
  if (!is.null(ezdat)) {
    if (!is.data.frame(ezdat)) {
      tryCatch({
        ezdat <- as.data.frame(ezdat)
      }, error = function(e) {
        stop("Could not coerce 'ezdat' to data.frame: ", e$message)
      })
    }
    ezdat_extracted <- tryCatch({
      stats::model.frame(z_formula, data = ezdat)
    }, error = function(e) {
      stop("Error extracting z variables from ezdat: ", e$message)
    })
    if (verbose) cat("  Evaluation data provided: ", nrow(ezdat_extracted), " evaluation points\n", sep = "")
  } else {
    ezdat_extracted <- z_extracted
    if (verbose) cat("  Using estimation data for evaluation (same as z data)\n")
  }
  
  # ========== DATA VALIDATION ==========

  # y_extracted is a vector, so nrow(y_extracted) is NULL.
  # Compare x and z row counts to length(y) instead.
  n_y <- length(y_extracted)
  n_x <- nrow(x_extracted)
  n_z <- nrow(z_extracted)

  if (is.null(n_x) || is.null(n_z)) {
    stop(
      "Internal error: x or z extraction did not return a row-oriented object"
    )
  }

  if (!isTRUE(n_x == n_y) || !isTRUE(n_z == n_y)) {
    stop(
      "Number of observations differ across y, x, and z variables after cleaning: ",
      "length(y) = ", n_y,
      ", nrow(x) = ", n_x,
      ", nrow(z) = ", n_z
    )
  }

  # ========== PREPARE OBJECTS FOR ESTIMATION ==========
  
  y <- y_extracted
  x_tz <- x_extracted
  z_tz <- z_extracted
  z_ez <- ezdat_extracted
  n_obs <- length(y)
  x_vars <- colnames(x_tz)
  z_var_names <- colnames(z_tz)
  
  # Determine factor status of z variables
  z_is_factor <- sapply(z_tz, is.factor)
  z_is_ordered <- sapply(z_tz, is.ordered)
  
  # Store parsed formula information
  formula_parsed <- list(
    original_formula = formula,
    y_formula = y_formula,
    x_formula = x_formula,
    z_formula = z_formula,
    inefficiency_formula = inefficiency_formula
  )
  
  # Store data information
  data_info <- list(
    data_supplied = data_supplied,
    n_total = n_total,
    n_used = n_used,
    n_removed = n_removed,
    variables_in_formula = all.vars(parsed_formula),
    y_name = deparse(y_formula[[2]]),
    x_names = colnames(x_extracted),
    z_names = colnames(z_extracted),
    inefficiency_names = if (!is.null(inefficiency_vars_data)) colnames(inefficiency_vars_data) else NULL
  )
  
  if (verbose) cat("\n========== DATA PREPARATION COMPLETE ==========\n\n")
  
  # ========== STEP 1: SMOOTH-COEFFICIENT ESTIMATION ==========
  
  if (verbose) cat("\n========== STEP 1: SMOOTH-COEFFICIENT ESTIMATION ==========\n")
  
  # Build data for npscoef (use training data for bandwidth selection)
  data_np_tz <- data.frame(y = y)
  for (i in 1:ncol(x_tz)) data_np_tz[[x_vars[i]]] <- x_tz[, i]
  for (i in 1:ncol(z_tz)) data_np_tz[[z_var_names[i]]] <- z_tz[, i]
  
  # Build formula for npscoef: y ~ x_vars | z_vars
  x_part <- paste(x_vars, collapse = " + ")
  z_part <- paste(z_var_names, collapse = " + ")
  formula_npscoef_str <- paste("y ~", x_part, "|", z_part)
  formula_npscoef <- as.formula(formula_npscoef_str)
  
  if (verbose) {
    cat("Building nonparametric model formula...\n")
    cat("Formula:", as.character(formula_npscoef)[2], "~", 
        as.character(formula_npscoef)[3], "\n")
  }
  
  # Bandwidth selection using training data, unless supplied by the user
  bandwidth_user_supplied <- !is.null(bandwidth)

  if (is.null(bandwidth)) {
    if (verbose) cat("Selecting bandwidth using LSCV (on training data)...\n")
    bw <- np::npscoefbw(formula = formula_npscoef, 
                    data = data_np_tz,
                    bwmethod = "cv.ls")
    if (verbose) cat("Bandwidth selected successfully by LSCV.\n")
  } else {
    bw <- np::npscoefbw(formula = formula_npscoef, 
                    bws = bandwidth,
                    data = data_np_tz,
                    bandwidth.compute = FALSE)
    if (verbose) cat("Using user-supplied bandwidth.\n")
  }
  
  # Estimate smooth coefficients
  if (verbose) cat("Estimating smooth coefficients...\n")
  
  # Check if evaluation data is different from training
  use_ezdat <- !is.null(ezdat) && !identical(z_ez, z_tz)
  
  if (use_ezdat) {

    # Build a complete evaluation data frame with all variables.
    # y and x are set to placeholders (means of training data) because
    # they are not used for estimating the smooth coefficients; only z matters.
    # The number of rows equals the number of evaluation points.
    n_eval <- nrow(z_ez)
    
    # Start with y (placeholders)
    data_np_ez <- data.frame(y = rep(mean(y), n_eval))
    
    # Add x columns with placeholders (means of each x variable)
    for (i in 1:ncol(x_tz)) {
      data_np_ez[[x_vars[i]]] <- rep(mean(x_tz[, i]), n_eval)
    }
    
    # Add z columns with the actual evaluation values
    for (i in 1:ncol(z_ez)) {
      data_np_ez[[z_var_names[i]]] <- z_ez[, i]
    }
    
    model_np <- np::npscoef(bws = bw, 
                        betas = TRUE, 
                        residuals = TRUE, 
                        data = data_np_tz,   # training data for estimation
                        newdata = data_np_ez)  # evaluation points
    if (verbose) cat("Smooth coefficients estimated on evaluation data.\n")
  } else {
    model_np <- np::npscoef(bws = bw, 
                        betas = TRUE, 
                        residuals = TRUE, 
                        data = data_np_tz)
    if (verbose) cat("Smooth coefficients estimated on training data.\n")
  }
  
  # Extract smooth coefficients
  coef_np <- coef(model_np)
  resid_np <- residuals(model_np)

  # when evaluation data are supplied, coef_np is evaluated at ezdat points,
  # which may differ in number from the training observations. Estimate the
  # smooth coefficients at the training observations separately so that
  # theta(Z_i) can be combined elementwise with Eu_mle to recover alpha(Z_i)
  if (use_ezdat) {
    model_np_tz <- np::npscoef(
      bws = bw,
      betas = TRUE,
      residuals = FALSE,
      data = data_np_tz
    )
    coef_np_tz <- coef(model_np_tz)
  } else {
    coef_np_tz <- coef_np
  }
  
  # Label the first smooth coefficient as the intercept theta(z)
  colnames(coef_np)[1] <- "Intercept (theta)"
  colnames(coef_np_tz)[1] <- "Intercept (theta)"
  
  if (verbose) {
    cat("\nSmooth Coefficients Summary:\n")
    print(summary(coef_np))
  }
  
  # ========== STEP 2: MAXIMUM LIKELIHOOD ESTIMATION ==========
  
  if (verbose) cat("\n========== STEP 2: MAXIMUM LIKELIHOOD ESTIMATION ==========\n")
  
  # Prepare Z matrix for MLE from inefficiency_vars_data
  # Convert categorical variables to binary (treatment contrasts)
  
  z_cat_warning_vars <- c()
  
  # Process inefficiency variables (convert categoricals to binary)
  ineff_var_names <- colnames(inefficiency_vars_data)
  ineff_is_cat <- sapply(inefficiency_vars_data, function(col) {
    is.factor(col) | is.character(col)
  })
  
  # Prepare Z matrix for MLE
  Z_ml_list <- list()
  for (i in seq_along(ineff_var_names)) {
    var_name <- ineff_var_names[i]
    is_cat <- ineff_is_cat[i]
    var_data <- inefficiency_vars_data[[var_name]]
    
    if (is_cat) {
      if (verbose) cat("Converting categorical variable '", var_name, "' to binary variables...\n", sep = "")
      z_cat_warning_vars <- c(z_cat_warning_vars, var_name)
      if (is.factor(var_data)) {
        levels_cat <- levels(var_data)
      } else {
        levels_cat <- unique(as.character(var_data))
        levels_cat <- sort(levels_cat)
      }
      for (k in 2:length(levels_cat)) {
        binary_var <- as.numeric(as.character(var_data) == levels_cat[k])
        Z_ml_list[[paste0(var_name, "_", levels_cat[k])]] <- binary_var
      }
    } else {
      Z_ml_list[[var_name]] <- var_data
    }
  }
  
  if (length(Z_ml_list) > 0) {
    Z_ml <- do.call(cbind, Z_ml_list)
  } else {
    Z_ml <- matrix(0, nrow = n_obs, ncol = 0)
  }
  Z_ml <- as.matrix(Z_ml)   # ensure matrix type
  n_z_ml <- ncol(Z_ml)
  
  # Issue warning if categorical variables were converted
  if (length(z_cat_warning_vars) > 0) {
    warning_msg <- paste(
      "NOTE: The following categorical variables in the inefficiency equation were converted to binary variables:\n",
      paste("  -", z_cat_warning_vars, collapse = "\n"),
      "\n"
    )
    if (verbose) cat("\n", warning_msg, "\n")
  } else {
    warning_msg <- NULL
  }
  
  # ========== MLE WITH HETEROGENEOUS INEFFICIENCY VARIANCE ==========
  
  resid_mle <- resid_np
  
  # Negative log-likelihood function (half-normal)
  llk_het <- function(param, y_resid, z_ml, n_z) {

    N <- length(y_resid)

    if (n_z == 0) {
      # No inefficiency variables: delta is scalar (intercept only)
      delta_0 <- param[1]  # Intercept only
      sigma2_v <- param[2]
      
      # Ensure positivity
      if (sigma2_v <= 1e-6) return(penalty)
      
      # sigma_u = exp(delta_0)
      s1 <- rep(exp(delta_0), N)
      s2 <- s1^2
      
    } else {
      # Has inefficiency variables
      z_ml <- cbind(1, z_ml)  # Add intercept
      n_z_full <- ncol(z_ml)
      
      # Extract parameters
      delta <- param[1:n_z_full]
      sigma2_v <- param[n_z_full + 1]
      
      # Ensure positivity
      if (sigma2_v <= 1e-6) return(penalty)
      
      # sigma_u = exp(delta'*z)
      s1 <- exp(z_ml %*% delta)
      s2 <- s1^2
    }
    
    if (any(s1 <= 0) | any(s2 <= 0)) return(penalty)
    
    sigma2 <- sigma2_v + s2
    sigma <- sqrt(sigma2)
    sigma_v <- sqrt(sigma2_v)
    lambda <- s1 / sigma_v
    
    # Adjustment for half-normal
    epsilon_star <- y_resid - s1 * sqrt(2 / pi)
    epsilon_star2 <- epsilon_star^2
    
    # Negative log-likelihood
    nll <- - N/2 * log(2/pi) + 0.5 * sum(log(s2 + sigma2_v) - 
                       2 * log(pnorm(-epsilon_star * lambda / sigma)) + 
                       epsilon_star2 / sigma2)
    
    return(nll)
  }
  
  # Initial values
  init_delta <- rep(0, n_z_ml + 1)  # n_z_ml parameters + intercept
  init_sigma2_v <- var(resid_mle) * 0.5
  
  if(any(stval=="automatic")) {
      initial_params <- c(init_delta, init_sigma2_v)
  } else {
      initial_params <- stval
  }
  
  if (verbose) {
    cat("Optimization method:", optim_method, "\n")
    if (n_z_ml == 0) {
      cat("Delta structure: scalar (no inefficiency variables)\n")
    } else {
      cat("Delta structure: intercept + ", n_z_ml, " variable(s)\n", sep = "")
    }

    if (any(stval=="automatic")) {
      cat("Initial sigma2_v:", round(init_sigma2_v, 6), "\n")
    } else {
      cat("Initial sigma2_v:", stval[length(stval)], "\n")
    }

    cat("Optimization starting...\n")
  }
  
  # Optimize
  if (optim_method == "optim") {
    optim_result <- optim(initial_params, 
                          llk_het,
                          method = "Nelder-Mead",
                          hessian = TRUE,
                          control = optim_control,
                          y_resid = resid_mle,
                          z_ml = Z_ml,
                          n_z = n_z_ml)
    
    if (optim_result$convergence != 0) {
      warning("Optimization did not converge. Convergence code:", optim_result$convergence)
    }
    
    param_mle <- optim_result$par
    loglik_value <- -optim_result$value
    optim_obj <- optim_result
    
  } else if (optim_method == "optim_BFGS") {
    optim_result <- optim(initial_params, 
                          llk_het,
                          method = "BFGS",
                          hessian = TRUE,
                          control = optim_control,
                          y_resid = resid_mle,
                          z_ml = Z_ml,
                          n_z = n_z_ml)
    
    if (optim_result$convergence != 0) {
      warning("Optimization did not converge. Convergence code:", optim_result$convergence)
    }
    
    param_mle <- optim_result$par
    loglik_value <- -optim_result$value
    optim_obj <- optim_result
    
    } else if (optim_method == "nlm") {
    nlm_result <- nlm(f = function(param) llk_het(param, resid_mle, Z_ml, n_z_ml),
                      p = initial_params,
                      hessian = TRUE,
                      iterlim = if (!is.null(optim_control$maxit)) optim_control$maxit else 100)
    
    param_mle <- nlm_result$estimate
    loglik_value <- -nlm_result$minimum
    optim_obj <- nlm_result
    
    if (verbose) cat("NLM iterations:", nlm_result$iterations, "\n")
    
  } else {
    stop("optim_method must be 'optim', 'optim_BFGS' or 'nlm'")
  }
  
# After optimization
if (optim_method %in% c("optim", "optim_BFGS")) {
    hessian <- optim_result$hessian
} else if (optim_method == "nlm") {
    hessian <- nlm_result$hessian
}

# Try to invert; if singular, use regularized Hessian
se_mle <- tryCatch({
    sqrt(diag(solve(hessian)))
}, error = function(e) {
    warning("Hessian is computationally singular. Using regularized Hessian (adding 1e-6 to diagonal) for standard errors.")
    hessian_reg <- hessian + diag(1e-6, nrow(hessian))
    # In case regularization still fails, return NA
    tryCatch({
        sqrt(diag(solve(hessian_reg)))
    }, error = function(e2) {
        warning("Regularized Hessian also singular. Standard errors set to NA.")
        rep(NA, length(param_mle))
    })
})

  # Extract MLE parameters
  n_delta <- n_z_ml + 1  # intercept + Z variables
  delta_mle <- param_mle[1:n_delta]
  delta_se <- se_mle[1:n_delta]
  sigma2_v_mle <- param_mle[n_delta + 1]
  sigma_v_mle <- sqrt(sigma2_v_mle)
  sigma2_v_se <- se_mle[n_delta + 1]
  
  # Calculate sigma_u for each observation
  if (n_z_ml == 0) {
    # Scalar delta
    delta_0 <- delta_mle[1]
    s1_mle <- rep(exp(delta_0), n_obs)
  } else {
    Z_ml_full <- cbind(1, Z_ml)
    s1_mle <- exp(Z_ml_full %*% delta_mle)
  }
  
  s2_mle <- s1_mle^2

  # Conditional mean of inefficiency and intercept alpha(z)
  Eu_mle <- s1_mle * sqrt(2 / pi)
  alpha <- as.matrix(coef_np_tz[, 1] + Eu_mle)
  colnames(alpha) <- "Intercept (alpha)"

  if (verbose) {
    cat("\nIntercept (alpha) Summary:\n")
    print(summary(alpha))
  }
  
  # Adjusted residuals
  epsilon_star_mle <- resid_mle - s1_mle * sqrt(2 / pi)
  
  # ========== BUILD MLE RESULTS TABLE ==========
  
  # Delta coefficients table
  if (n_z_ml == 0) {
    delta_names <- c("Intercept (delta_0)")
  } else {
    delta_names <- c("Intercept", paste0("z", 1:n_z_ml))
    if (!is.null(colnames(Z_ml)) & length(colnames(Z_ml)) == n_z_ml) {
      delta_names <- c("Intercept", colnames(Z_ml))
    }
  }
  
  t_stat_delta <- delta_mle / delta_se
  p_val_delta <- 2 * (1 - pnorm(abs(t_stat_delta)))
  
  mle_table_delta <- data.frame(
    Coefficient = delta_names,
    Estimate = delta_mle,
    Std.Error = delta_se,
    t.value = t_stat_delta,
    p.value = p_val_delta,
    row.names = NULL
  )
  
  # Sigma2_v and Log-likelihood
  t_stat_sigma2_v <- sigma2_v_mle / sigma2_v_se
  p_val_sigma2_v <- 2 * (1 - pnorm(abs(t_stat_sigma2_v)))
  
  if (n_z_ml == 0) {
    formula_desc <- "exp(delta_0) (scalar)"
  } else {
    formula_desc <- "exp(delta_0+delta_1'z)"
  }
  
  mle_summary_table <- data.frame(
    Parameter = c("sigma2_v", "sigma_u formula", "Log-Likelihood"),
    Estimate = c(sigma2_v_mle, formula_desc, loglik_value),
    Std.Error = c(sigma2_v_se, "", ""),
    t.value = c(t_stat_sigma2_v, "", ""),
    p.value = c(p_val_sigma2_v, "", ""),
    row.names = NULL
  )
  
  if (verbose) {
    cat("\n========== MLE RESULTS ==========\n")
    if (n_z_ml == 0) {
      cat("\nInefficiency Parameters (sigma_u = exp(delta_0) (scalar)):\n")
    } else {
      cat("\nInefficiency Determinants (sigma_u = exp(delta_0+delta_1'z)):\n")
    }
    print(mle_table_delta)
    cat("\nVariance Parameters:\n")
    print(mle_summary_table)
  }
  
  # ========== TECHNICAL EFFICIENCY CALCULATION (BC METHOD) ==========
  
  if (verbose) cat("\n========== TECHNICAL EFFICIENCY (BC METHOD) ==========\n")
  
  BC_efficiency <- function(epsilon, sigma2_v, sigma2_u, sig_level = 0.05) {
    sigma_v <- sqrt(sigma2_v)
    sigma_u <- sqrt(sigma2_u)
    sigma2 <- sigma2_v + sigma2_u
    sigma <- sqrt(sigma2)
    sigma_star <- sigma_v * sigma_u / sigma
    mu_star <- -epsilon * sigma2_u / sigma2
    
    cdf1 <- pnorm(sigma_star - (mu_star / sigma_star))
    cdf2 <- pnorm(-mu_star / sigma_star)
    
    TE <- as.vector(((1 - cdf1) / (1 - cdf2)) * exp(-mu_star + 0.5 * sigma_star^2))
    
    Z_L <- qnorm(1 - (sig_level / 2) * (1 - cdf2))
    Z_U <- qnorm(1 - (1 - sig_level / 2) * (1 - cdf2))
    
    TE_L <- as.vector(exp(-mu_star - Z_L * sigma_star))
    TE_U <- as.vector(exp(-mu_star - Z_U * sigma_star))
    
    list(TE = TE, TE_L = TE_L, TE_U = TE_U)
  }
  
  te_bc <- BC_efficiency(epsilon = epsilon_star_mle, 
                         sigma2_v = sigma2_v_mle, 
                         sigma2_u = s2_mle, 
                         sig_level = sig_level)
  
  TE <- te_bc$TE
  TE_L <- te_bc$TE_L
  TE_U <- te_bc$TE_U
  
  if (verbose) {
    cat("Technical Efficiency Summary:\n")
    cat("Mean TE:", round(mean(TE), 6), "\n")
    cat("Std Dev TE:", round(sd(TE), 6), "\n")
    cat("Min TE:", round(min(TE), 6), "\n")
    cat("Max TE:", round(max(TE), 6), "\n")
  }
  
  # ========== INEFFICIENCY CALCULATION (JLMS METHOD) ==========
  
  if (verbose) cat("\n========== INEFFICIENCY (JLMS METHOD) ==========\n")
  
  JLMS_inefficiency <- function(epsilon, sigma_v, sigma_u, mu = 0) {
    sigma2_v <- sigma_v^2
    sigma2_u <- sigma_u^2
    sigma2 <- sigma2_v + sigma2_u
    sigma <- sqrt(sigma2)
    sigma_star <- sigma_v * sigma_u / sigma
    mu_tilde <- (-sigma2_u * epsilon + mu * sigma2_v) / sigma2
    m <- mu_tilde / sigma_star
    
    cdf <- pnorm(m)
    pdf <- dnorm(m)
    
    as.vector(mu_tilde + sigma_star * pdf / cdf)
  }
  
  inefficiency <- JLMS_inefficiency(epsilon = epsilon_star_mle,
                                    sigma_v = sigma_v_mle,
                                    sigma_u = s1_mle,
                                    mu = 0)
  
  if (verbose) {
    cat("Inefficiency Summary:\n")
    cat("Mean Inefficiency:", round(mean(inefficiency), 6), "\n")
    cat("Std Dev Inefficiency:", round(sd(inefficiency), 6), "\n")
    cat("Min Inefficiency:", round(min(inefficiency), 6), "\n")
    cat("Max Inefficiency:", round(max(inefficiency), 6), "\n")
  }
  
  # ========== COMPILE RESULTS ==========
  
  results <- list(
    # Step 1 Results
    smooth_coefficients = list(
      theta_beta = coef_np,
      alpha = alpha,
      residuals = resid_np,
      bandwidth = bw,
      bandwidth_user_supplied = bandwidth_user_supplied,
      formula = formula_npscoef,
      n_training = nrow(data_np_tz),
      n_evaluation = if (use_ezdat) nrow(data_np_ez) else nrow(data_np_tz),
      used_evaluation_data = use_ezdat
    ),
    
    # Step 2 Results
    mle_results = list(
      delta_table = mle_table_delta,
      summary_table = mle_summary_table,
      Eu_mle = Eu_mle,
      delta_coef = delta_mle,
      delta_se = delta_se,
      n_z = n_z_ml,
      sigma2_v = sigma2_v_mle,
      sigma2_v_se = sigma2_v_se,
      log_likelihood = loglik_value,
      optim_method = optim_method,
      optim_result = optim_obj
    ),
    
    # Efficiency Results
    technical_efficiency = list(
      TE = TE,
      TE_lower = TE_L,
      TE_upper = TE_U,
      sig_level = sig_level,
      summary = c(
        Mean = mean(TE),
        SD = sd(TE),
        Min = min(TE),
        Max = max(TE),
        Median = median(TE)
      )
    ),
    
    inefficiency = list(
      u_JLMS = inefficiency,
      summary = c(
        Mean = mean(inefficiency),
        SD = sd(inefficiency),
        Min = min(inefficiency),
        Max = max(inefficiency),
        Median = median(inefficiency)
      )
    ),
    
    # Model Info
    model_summary = list(
      formula = formula,
      n_obs = n_obs,
      n_params_x = length(x_vars),
      n_params_z = length(z_var_names),
      n_cat_vars_z = sum(z_is_factor),
      x_names = x_vars,
      z_names = z_var_names,
      z_is_factor = z_is_factor,
      z_is_ordered = z_is_ordered,
      inefficiency_vars_na = FALSE,  # no longer used
      inefficiency_vars_null = is.null(inefficiency_formula),  # TRUE if no third part
      n_inefficiency_vars = n_z_ml,
      categorical_conversion_warning = warning_msg,
      sig_level = sig_level,
      data_supplied = data_supplied,
      bandwidth_user_supplied = bandwidth_user_supplied
    ),
    
    # Data used information
    data_info = data_info,
    formula_parsed = formula_parsed
  )
  
  class(results) <- "spsc_sfa"
  
  if (verbose) cat("\n========== ANALYSIS COMPLETE ==========\n")
  
  return(results)
}
