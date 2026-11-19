#' Semiparametric Smooth-Coefficient Stochastic Frontier Analysis
#'
#' Estimates a semiparametric smooth-coefficient stochastic frontier model
#' \eqn{y_i = \alpha(Z_i) + X_i' \beta(Z_i) + v_i - u(Z_i)}
#' following Sun and Kumbhakar (2013). The model allows the coefficients of the
#' parametric part to vary smoothly with a set of nonparametric variables.
#' Technical efficiency and inefficiency are computed using the Battese and
#' Coelli (1988) (BC) and Jondrow et al. (1982) (JLMS) methods, respectively.
#' Confidence intervals for technical efficiency are computed using the
#' approach of Horrace and Schmidt (1996). Marginal effects of the inefficiency 
#' determinants on the JLMS inefficiency are computed following Kumbhakar and Sun (2013).
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
#' The smooth coefficients of the first step can be estimated by either a local-constant
#' (LC) or a local-linear (LL) kernel estimator, as chosen by \code{spscreg}.
#' With \code{spscreg = "lc"} (default), the bandwidths are selected by
#' \code{np::npscoefbw} and the smooth coefficients are estimated by
#' \code{np::npscoef} with default arguments. With \code{spscreg = "ll"}, the
#' bandwidths are selected by \code{\link{llscoefbw}} and the smooth
#' coefficients are estimated by \code{\link{llscoef}}. If \code{bandwidth} is
#' \code{NULL} (default), the bandwidths are selected by least-squares
#' cross-validation (LSCV) for the chosen estimator. The residuals of the chosen
#' estimator are used in the second step, the maximum likelihood estimation.
#'
#' If \code{gradients = TRUE}, the gradients of the first-step smooth
#' coefficients with respect to the continuous variables in \code{z} are
#' estimated as well: by \code{\link{lcscoefgrad}} for the local-constant
#' estimator and by \code{\link{llscoef}} for the local-linear estimator. Since
#' \eqn{\alpha(Z_i) = \theta(Z_i) + E(u_i \mid Z_i)}, the gradient of
#' \eqn{\alpha(Z_i)} with respect to the \eqn{l}th continuous \eqn{Z} variable is
#' \deqn{
#' \partial\alpha(Z_i)/\partial z_l = \partial\theta(Z_i)/\partial z_l
#' + \sqrt{2/\pi}\,\exp(\delta_0+\delta_1'Z_i)\,\delta_{1l},
#' }{
#' d alpha(Z_i) / d z_l = d theta(Z_i) / d z_l
#' + sqrt(2/pi) * exp(delta_0 + delta_1' Z_i) * delta_1l,
#' }
#' where \eqn{\delta_{1l}} is the coefficient of that variable in the
#' inefficiency equation (zero if the variable is not an inefficiency
#' determinant) and \eqn{\partial\theta(Z_i)/\partial z_l} is evaluated at the
#' training data.
#' 
#' The analytic gradients of the local-constant smooth coefficient estimator
#' with respect to the continuous \code{z} variables, returned in
#' \code{smooth_coefficients_grad} if \code{gradients = TRUE} and
#' \code{spscreg = "lc"}, are computed following Geng and Sun (2019); see
#' \code{\link{lcscoefgrad}}. The marginal effects of the inefficiency determinants on the JLMS
#' inefficiency estimates, returned as \code{inefficiency_grad} if
#' \code{gradients = TRUE}, are computed following Kumbhakar and Sun (2013).
#'
#' The JLMS inefficiency estimator is the conditional mean
#' \eqn{E(u_i \mid \varepsilon_i) = \tilde\mu_i + \sigma_{*i}\, g_i}, where
#' \eqn{\tilde\mu_i = -\varepsilon_i \sigma_{ui}^2 / \sigma_i^2},
#' \eqn{\sigma_{*i} = \sigma_{ui}\sigma_v / \sigma_i},
#' \eqn{\sigma_i^2 = \sigma_{ui}^2 + \sigma_v^2},
#' \eqn{m_i = \tilde\mu_i / \sigma_{*i}} and
#' \eqn{g_i = \phi(m_i) / \Phi(m_i)}. Since \eqn{\sigma_{ui} = \exp(\delta_0 +
#' \delta_1' Z_i)} depends on the inefficiency determinants, the first-order
#' analytic derivative of the JLMS estimator with respect to the \eqn{l}-th
#' inefficiency determinant (holding \eqn{\varepsilon_i} fixed) is
#' \deqn{
#' \frac{\partial E(u_i \mid \varepsilon_i)}{\partial z_{li}}
#' = \frac{\delta_{1l}}{\sigma_i^2}
#' \Bigl[\sigma_v^2 \sigma_{*i} \bigl(g_i (1 + m_i^2) + m_i g_i^2\bigr)
#' - 2 \sigma_{*i}^2 \varepsilon_i (1 - g_i^2 - m_i g_i)\Bigr].
#' }{
#' dE(u_i | e_i) / dz_li = (delta_1l / sigma_i^2) * [ sigma_v^2 * sigma_*i *
#' (g_i (1 + m_i^2) + m_i g_i^2) - 2 * sigma_*i^2 * e_i * (1 - g_i^2 - m_i g_i) ].
#' }
#' This is the \eqn{\sigma_{ui}} channel of Equation (9) in Kumbhakar and Sun
#' (2013), where the coefficient on \eqn{z_{li}} in \eqn{\sigma_{ui}} was called
#' \eqn{\gamma_l} (it is \eqn{\delta_{1l}} here, as in Sun and Kumbhakar, 2013);
#' the other two channels vanish because \eqn{\mu_i = 0} and \eqn{\sigma_v} is
#' constant. If \code{gradients = TRUE}, the derivatives are returned in
#' \code{inefficiency_grad}. A column of \code{inefficiency_grad} that belongs
#' to a dummy variable created from a categorical inefficiency determinant is
#' the analytic partial derivative with the dummy treated as continuous.
#'
#' @param formula A formula of the form \code{y ~ x | z} or
#'   \code{y ~ x | z | q},
#'   where \code{y} is the dependent variable,
#'   \code{x} are parametric inputs,
#'   \code{z} are nonparametric variables determining smooth coefficients,
#'   and \code{q} (optional) are inefficiency determinants.
#'   When the third formula part \code{q} is omitted, the \code{z} variables also 
#'   serve as inefficiency determinants.
#'
#' @param data Optional data frame containing all variables used in \code{formula}.
#'   If supplied, variables are extracted from this data frame.
#'   If not supplied, variables are taken from the environment of \code{formula}.
#'
#' @param ezdat Optional data frame containing evaluation points for the
#'   nonparametric variables \code{z}. If provided, smooth coefficients
#'   are estimated at these points; otherwise, the training data points
#'   (from \code{z}) are used.
#'
#' @param spscreg Character string; the estimator of the smooth coefficients
#'   in the first step. \code{"lc"} (default) is the local-constant estimator
#'   (\code{np::npscoefbw} and \code{np::npscoef} both being made with default 
#'   \code{regtype = "lc"} and default kernels, along with their default 
#'   \code{bwtype = "fixed"} and \code{bwscaling = FALSE}) and \code{"ll"} is the 
#'   local-linear estimator (\code{\link{llscoefbw}} and \code{\link{llscoef}}).
#'
#' @param bandwidth Optional bandwidths for the smooth-coefficient estimation. 
#'   If \code{NULL} (default and recommended), bandwidths are selected by least-squares 
#'   cross-validation (LSCV): for the local-constant estimator if \code{spscreg = "lc"}, and for 
#'   the local-linear estimator if \code{spscreg = "ll"}. Or a numeric vector with one 
#'   bandwidth for each variable in \code{z}, in the order of the formula, of the kind 
#'   that suits \code{spscreg} (local-constant bandwidths for \code{"lc"}, local-linear 
#'   ones for \code{"ll"}). A named vector, such as the one returned by \code{llscoefbw} 
#'   with \code{LSCV = FALSE}, is matched to the \code{z} variables by
#'   name if its names are those of the \code{z} variables, and is taken in the
#'   order of the formula otherwise. Or a bandwidth object: for \code{spscreg = "lc"} from
#'   \code{np::npscoefbw} with \code{regtype = "lc"} (its default), and for
#'   \code{spscreg = "ll"} from \code{np::npscoefbw} with
#'   \code{regtype = "ll"}, with both using fixed bandwidths (\code{bwtype = "fixed"},
#'   \code{bwscaling = FALSE}) and the default kernels; or from
#'   \code{\link{llscoefbw}}. An object of the other kind is an error, as the
#'   bandwidths of the two estimators are not interchangeable; a vector of the
#'   other kind cannot be told apart and is not suitable.
#'
#' @param gradients Logical; if \code{TRUE}, the gradients of the smooth
#'   coefficients with respect to the continuous \code{z} variables (local-constant
#'   or local-linear, as chosen by \code{spscreg}) and the gradients of the JLMS
#'   inefficiency estimates are computed and returned. If \code{FALSE}
#'   (default), no gradient is computed or returned.
#'
#' @param verbose Logical; if \code{TRUE}, progress and intermediate results
#'   are printed.
#'
#' @param sig_level Significance level for confidence intervals of technical
#'   efficiency (default 0.05).
#'
#' @param optim_method Character string specifying the optimizer for MLE in the second step:
#'   \code{"optim"} (uses Nelder-Mead), \code{"optim_BFGS"} (uses BFGS,
#'   the default), or \code{"nlm"}.
#'
#' @param optim_control List of control parameters passed to the optimizer for MLE in the second step.
#'
#' @param stval Starting values for the MLE in the second step. Use \code{"automatic"} (default)
#'   for internally computed starting values, or supply a numeric vector 
#'   in which the order of its elements is \eqn{\delta_0}, \eqn{\delta_1} (if any), and \eqn{\sigma^2_v}.
#'   See Details.
#'
#' @param penalty To ensure positivity of estimated \eqn{\sigma^2_v}, penalty
#'   applied when \eqn{\sigma^2_v} becomes negative during a numerical search
#'   (default 1e8).
#'
#' @return A list of class \code{"spscsfa"} containing:
#'   \item{smooth_coefficients}{Smooth coefficient estimates (local-constant if
#'        \code{spscreg = "lc"}, local-linear if \code{spscreg = "ll"}),
#'        residuals, bandwidth, and metadata.}
#'   \item{smooth_coefficients_grad}{If \code{gradients = TRUE}. A list
#'        whose first element, \code{theta_beta_grad}, is itself a list with one
#'        matrix for each continuous \code{z} variable (named after it): the
#'        gradients of \eqn{\theta(z)} and \eqn{\beta(z)} with respect to that
#'        variable, at the evaluation points, with columns
#'        \code{d(Intercept (theta))/d(z)}, \code{d(x1)/d(z)}, \ldots. Like \code{x1} 
#'        used in a regression table, the \code{x1} here refers to the smooth coefficient
#'        of \code{x1}. Its second element, \code{alpha_grad}, is a matrix with one column for
#'        each continuous \code{z} variable, named \code{d(Intercept (alpha))/d(z)}: 
#'        the gradients of \eqn{\alpha(Z_i)} at the training data.}
#'   \item{mle_results}{MLE parameter estimates, standard errors, and
#'        log-likelihood.}
#'   \item{technical_efficiency}{Technical efficiency scores (BC method)
#'        and confidence intervals.}
#'   \item{inefficiency}{Inefficiency estimates (JLMS method).}
#'   \item{inefficiency_grad}{If \code{gradients = TRUE}. A matrix of the
#'        first-order analytic derivatives of the JLMS inefficiency estimates
#'        with respect to the inefficiency determinants (Kumbhakar and Sun,
#'        2013). It has one row per training observation and one column per
#'        inefficiency determinant used in the maximum likelihood step, named
#'        after the columns of that matrix (categorical variables appear as
#'        their binary variables). It is \code{NULL} when there are no
#'        inefficiency determinants (intercept only).}
#'   \item{model_summary}{Summary information about the model specification.}
#'   \item{formula_parsed}{Parsed formula components.}
#'   \item{data_info}{Information about data cleaning and observations used.}
#'   \item{model_frame}{The estimation sample: the model frame, without missing
#'        values, of the variables used in \code{formula}. It is kept so that
#'        \code{\link{plot.spscsfa}} can re-evaluate the smooth coefficients
#'        on a grid of points.}
#' 
#' @note
#' Results can be summarized using \code{\link{summary.spscsfa}}.
#' Extraction functions for smooth coefficients, technical efficiency and
#' inefficiency are available through \code{\link{smoothcoef}},
#' \code{\link{te}}, and \code{\link{inefficiency}}, respectively. If
#' \code{gradients = TRUE}, the gradients are extracted with
#' \code{\link{grad}}.
#' The smooth coefficients can be plotted against the continuous nonparametric
#' variables using \code{\link{plot.spscsfa}}.
#'
#' \code{spsc_sfa} is an alias of \code{\link{spscsfa}}, \code{smooth_coef} is an alias of \code{\link{smoothcoef}},
#' and the methods of the class \code{"spscsfa"} also handle fitted models and
#' summaries of the class \code{"spsc_sfa"}.
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
#' Xin Geng and Kai Sun (2019).
#' Gradient estimation of the local-constant semiparametric smooth coefficient
#' model.
#' \emph{Economics Letters}, 185, 108684.
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
#' Subal C. Kumbhakar and Kai Sun (2013).
#' Derivation of marginal effects of determinants of technical inefficiency.
#' \emph{Economics Letters}, 120(2), 249--253.
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
#' result <- spscsfa(output ~ input | covariate, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#' result$smooth_coefficients
#' result$mle_results
#' result$technical_efficiency
#' result$inefficiency
#' 
#' # This gives the same results: the inefficiency variables default to the nonparametric part
#' result <- spscsfa(output ~ input | covariate | NULL, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # Use evaluation data (ezdat) for estimating the smooth coefficients
#' result <- spscsfa(output ~ input | covariate, data = dat,
#'                    ezdat = data.frame(covariate = ez),
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#' nrow(result$smooth_coefficients$theta_beta)
#' 
#' # Gradients of the smooth coefficients, of alpha, and of the inefficiency
#' result <- spscsfa(output ~ input | covariate, data = dat, gradients = TRUE, verbose = TRUE)
#' result$smooth_coefficients_grad$theta_beta_grad
#' head(result$smooth_coefficients_grad$alpha_grad)
#' head(result$inefficiency_grad)
#' grad(result, zvar = "covariate")
#'
#' # Three-part formula: nonparametric part differs from inefficiency determinants
#' result <- spscsfa(output ~ input | covariate | ineffvar, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # Customized bandwidth for the nonparametric part
#' result <- spscsfa(output ~ input | covariate | ineffvar, data = dat,
#'                    optim_method = "optim_BFGS", bandwidth = c(0.1), verbose = TRUE)
#'
#' # Restrict no inefficiency determinants (intercept only)
#' result <- spscsfa(output ~ input | covariate | 1, data = dat,
#'                    optim_method = "optim_BFGS", verbose = TRUE)
#'
#' # Local-linear smooth coefficients (the default is spscreg = "lc")
#' result <- spscsfa(output ~ input | covariate, data = dat, spscreg = "ll", verbose = TRUE)
#' result$smooth_coefficients$theta_beta
#'
#' # Gradients with local-linear smooth coefficients
#' result <- spscsfa(output ~ input | covariate, data = dat,
#'                    spscreg = "ll", gradients = TRUE, verbose = TRUE)
#' result
#'
#' # User-provided bandwidth for the local-linear estimator: llscoefbw() with
#' # LSCV = FALSE returns a named vector
#' bw <- llscoefbw(output ~ input | covariate, bws = 0.2, data = dat, LSCV = FALSE)
#' result <- spscsfa(output ~ input | covariate, data = dat,
#'                    spscreg = "ll", bandwidth = bw, verbose = TRUE)
#' \donttest{
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
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2, data = dat, spscreg = "ll", gradients = TRUE)
#'
#' # z3 is an unordered categorical variable
#' # spscsfa() detects if z3 is a factor
#' dat$z3 <- factor(rep(sample(c("A", "B", "C", "D", "E")), 20))
#'
#' # z1, z2, z3 enter the nonparametric part and are also inefficiency determinants
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2 + z3, data = dat)
#'
#' # z3 is preserved as an unordered categorical variable in the nonparametric part
#' result$smooth_coefficients$bandwidth
#' 
#' # z3 is converted to dummy variables as inefficiency determinants
#' result$mle_results$delta_table
#'
#' # Same as above (inefficiency defaults to nonparametric part)
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2 + z3 | NULL, data = dat)
#'
#' # Three-part formula: use z3 only as the inefficiency determinant
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2 + z3 | z3, data = dat)
#'
#' # Use external variables as inefficiency determinants
#' dat$q1 <- rnorm(n)
#' dat$q2 <- rnorm(n)
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2 + z3 | q1 + q2, data = dat)
#'
#' # No inefficiency determinants (intercept only)
#' result <- spscsfa(y ~ x1 + x2 | z1 + z2 + z3 | 1, data = dat, verbose = TRUE)
#'
#' # Application with real data (ricephil)
#' \dontrun{
#' 
#' data(ricephil)
#' 
#' ricephil$FARMERCODE <- as.factor(ricephil$FARMERCODE)
#' 
#' # this example shows categorical variable handling
#' result <- spscsfa(PROD ~ AREA + LABOR + NPK | EDYRS + FARMERCODE | FARMERCODE,
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
spscsfa <- function(formula, 
                     data = NULL,
                     ezdat = NULL,
                     spscreg = c("lc", "ll"),
                     bandwidth = NULL,
                     gradients = FALSE,
                     verbose = TRUE,
                     sig_level = 0.05,
                     optim_method = "optim_BFGS",
                     optim_control = list(maxit = 5000),
                     stval = "automatic",
                     penalty = 1e8) {
  
  # ========== INPUT VALIDATION ==========
  
  if (verbose) cat("\n========== INPUT VALIDATION ==========\n")

  # Estimator of the smooth coefficients: "lc" (local-constant, default) or
  # "ll" (local-linear)
  spscreg <- match.arg(spscreg)

  # The bandwidths of the two estimators are not interchangeable.
  if (inherits(bandwidth, "scbandwidth")) {
    if (spscreg == "ll") {
      .llscoef_np_bw(bandwidth, "bandwidth")
    } else {
      .lcscoef_np_bw(bandwidth, "bandwidth")
    }
  }

  if (spscreg == "lc" && inherits(bandwidth, "llsmoothcoefficientbw")) {
    stop("'bandwidth' is an object returned by llscoefbw(), whose bandwidths ",
         "are for the local-linear estimator: with spscreg = \"lc\" use an ",
         "object returned by npscoefbw() with regtype = \"lc\" or a numeric ",
         "vector of bandwidths", call. = FALSE)
  }

  if (!is.logical(gradients) || length(gradients) != 1L || is.na(gradients)) {
    stop("'gradients' must be TRUE or FALSE")
  }

  if (verbose) {
    cat("Smooth-coefficient estimator: ",
        if (spscreg == "ll") "local-linear (LL)" else "local-constant (LC)",
        "\n", sep = "")
    cat("Gradients: ", if (gradients) "requested" else "not requested",
        "\n", sep = "")
  }
  
  # ========== VALIDATE SIGNIFICANCE LEVEL ==========
  
  if (!is.numeric(sig_level) ||
      length(sig_level) != 1L ||
      !is.finite(sig_level) ||
      sig_level <= 0 ||
      sig_level >= 1) {
    
    stop(
      "'sig_level' must be a single finite number strictly between 0 and 1",
      call. = FALSE
    )
  }
  
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
  # q = inefficiency determinants.
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
  
  # The parts of the formula are read from the model frame built above and are
  # not evaluated again. Its columns are named as the terms are written in
  # 'formula'.

  # Extract y (response variable) from the model frame.
  y_extracted <- tryCatch({
    stats::model.response(clean_data)
  }, error = function(e) {
    stop("Could not extract response variable from formula: ", e$message)
  })

  if (!is.numeric(y_extracted)) {
    stop("Response variable (y) must be numeric")
  }
  
  if (verbose) cat("  y: ", length(y_extracted), " observations\n", sep = "")
  
  # Extract x (parametric variables) - as matrix
  x_extracted <- stats::model.matrix(stats::terms(x_formula), data = clean_data)
  
  # Remove intercept if present (we'll add it later if needed)
  if (colnames(x_extracted)[1] == "(Intercept)") {
    x_extracted <- x_extracted[, -1, drop = FALSE]
  }
  
  if (verbose) cat("  x: ", nrow(x_extracted), " observations, ", ncol(x_extracted), " variables\n", sep = "")
  
  # Extract z (nonparametric variables) - as data frame
  z_extracted <- Formula::model.part(parsed_formula, data = clean_data, rhs = 2)
  attr(z_extracted, "terms") <- NULL
  
  if (verbose) cat("  z: ", nrow(z_extracted), " observations, ", ncol(z_extracted), " variables\n", sep = "")
  
  # Extract inefficiency variables from formula's third part if present,
  # otherwise default to z variables.
  inefficiency_vars_data <- NULL
  if (!is.null(inefficiency_formula)) {
    inefficiency_vars_data <- Formula::model.part(parsed_formula,
                                                  data = clean_data, rhs = 3)
    attr(inefficiency_vars_data, "terms") <- NULL
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
  
  # Build syntactically named data for the first-step estimator; bandwidth
  # selection and estimation use the training observations.
  np_names <- make.names(c("y", x_vars, z_var_names), unique = TRUE)
  np_x <- np_names[1L + seq_along(x_vars)]
  np_z <- np_names[1L + length(x_vars) + seq_along(z_var_names)]

  data_np_tz <- data.frame(y = y)
  for (i in 1:ncol(x_tz)) data_np_tz[[np_x[i]]] <- x_tz[, i]
  for (i in 1:ncol(z_tz)) data_np_tz[[np_z[i]]] <- z_tz[, i]
  
  # Build the first-step formula: y ~ x_vars | z_vars
  x_part <- paste(np_x, collapse = " + ")
  z_part <- paste(np_z, collapse = " + ")
  formula_npscoef_str <- paste("y ~", x_part, "|", z_part)
  formula_npscoef <- as.formula(formula_npscoef_str)
  
  if (verbose) {
    cat("Building nonparametric model formula...\n")
    cat("Formula:", as.character(formula_npscoef)[2], "~", 
        as.character(formula_npscoef)[3], "\n")
  }
  
  bandwidth_user_supplied <- !is.null(bandwidth)

  # Evaluation data is different from training data
  use_ezdat <- !is.null(ezdat) && !identical(z_ez, z_tz)

  # Names of the smooth coefficients: the first one is the intercept theta(z)
  coef_names <- c("Intercept (theta)", x_vars)

  # Gradients are taken with respect to the continuous z variables only
  z_is_cont <- vapply(z_tz, is.numeric, logical(1))
  z_cont_names <- z_var_names[z_is_cont]
  n_z_cont <- length(z_cont_names)

  if (gradients && n_z_cont == 0L) {
    warning("'gradients = TRUE', but there is no continuous z variable: ",
            "the gradients of the smooth coefficients are not available",
            call. = FALSE)
  }

  # Gradients of theta and beta at the evaluation points (grad_tb) and at the
  # training data (grad_tb_tz, used for alpha): lists with one matrix for each
  # continuous z variable
  grad_tb <- NULL
  grad_tb_tz <- NULL

  # Evaluation data frame with all variables, if it differs from the training
  # data. y and x are set to placeholders (means of training data) because they
  # are not used for estimating the smooth coefficients; only z matters. The
  # number of rows equals the number of evaluation points.
  data_np_ez <- NULL

  if (use_ezdat) {
    n_eval <- nrow(z_ez)

    data_np_ez <- data.frame(y = rep(mean(y), n_eval))

    for (i in 1:ncol(x_tz)) {
      data_np_ez[[np_x[i]]] <- rep(mean(x_tz[, i]), n_eval)
    }

    # Add z columns with the actual evaluation values
    for (i in 1:ncol(z_ez)) {
      data_np_ez[[np_z[i]]] <- z_ez[, i]
    }
  }

  if (spscreg == "lc") {

    # ---------- local-constant smooth coefficients (np::npscoef) ----------

    # Bandwidth selection using training data
    if (is.null(bandwidth)) {
      if (verbose) cat("Selecting bandwidth using LSCV for the local-constant estimator (on training data)...\n")
      bw <- np::npscoefbw(formula = formula_npscoef, 
                      data = data_np_tz,
                      bwmethod = "cv.ls")
      if (verbose) cat("Bandwidth selected successfully by LSCV.\n")
    } else {
      # A bandwidth object of npscoefbw is used as it is; a numeric vector is
      # matched to the z variables
      if (!inherits(bandwidth, "scbandwidth")) {
        bandwidth <- .spsc_bw_vector(bandwidth, z_var_names, np_z)
      }
      bw <- np::npscoefbw(formula = formula_npscoef, 
                      bws = bandwidth,
                      data = data_np_tz,
                      bandwidth.compute = FALSE)
      if (verbose) cat("Using user-supplied bandwidth.\n")
    }
    
    # Estimate smooth coefficients
    if (verbose) cat("Estimating smooth coefficients...\n")
    
    if (use_ezdat) {
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

    # Print the smooth coefficient summary first
    colnames(coef_np) <- coef_names
    if (verbose) {
      cat("\nSmooth Coefficients Summary:\n")
      print(summary(coef_np))
    }

    # Gradients of the local-constant smooth coefficients: lcscoefgrad() returns
    # one matrix for each continuous z variable, with the gradients of
    # beta_1, ..., beta_p (the intercept first) with respect to that variable
    if (gradients && n_z_cont > 0L) {
      if (verbose) cat("\nEstimating gradients of the local-constant smooth coefficients...\n")
      grad_tb <- lcscoefgrad(bws = bw,
                             txdat = x_tz,
                             tydat = y,
                             tzdat = z_tz,
                             ezdat = if (use_ezdat) z_ez else NULL)
      grad_tb_tz <- if (use_ezdat) {
        lcscoefgrad(bws = bw, txdat = x_tz, tydat = y, tzdat = z_tz)
      } else {
        grad_tb
      }
    }

  } else {

    # ---------- local-linear smooth coefficients (llscoefbw, llscoef) ----------

    # Bandwidth selection using training data
    if (is.null(bandwidth)) {
      if (verbose) cat("Selecting bandwidth using LSCV for the local-linear estimator (on training data)...\n")
      bw <- llscoefbw(formula_npscoef, data = data_np_tz)
      if (verbose) cat("Bandwidth selected successfully by LSCV.\n")
    } else {
      # llscoefbw() with LSCV = FALSE checks the bandwidths and
      # returns them as a vector named after the z variables. llscoef() needs
      # a formula to be given such a vector.
      bw <- llscoefbw(formula_npscoef,
                      bws = .spsc_bw_vector(bandwidth, z_var_names, np_z),
                      data = data_np_tz,
                      LSCV = FALSE)
      if (verbose) cat("Using user-supplied bandwidth.\n")
    }

    # Estimate smooth coefficients (and their gradients). A bandwidth object
    # from llscoefbw() keeps its formula and llscoef() does not accept another
    # one with it; a vector of bandwidths (supplied by the user) needs it.
    ll_formula <- if (inherits(bw, "llsmoothcoefficientbw")) NULL else formula_npscoef

    if (verbose) cat("Estimating smooth coefficients...\n")

    fit_ll <- llscoef(bw,
                      formula = ll_formula,
                      data = data_np_tz,       # training data for estimation
                      newdata = data_np_ez,    # evaluation points (or NULL)
                      betas = TRUE,
                      residuals = TRUE,
                      gradients = gradients && n_z_cont > 0L)

    if (verbose) {
      if (use_ezdat) {
        cat("Smooth coefficients estimated on evaluation data.\n")
      } else {
        cat("Smooth coefficients estimated on training data.\n")
      }
    }

    coef_np <- coef(fit_ll)
    resid_np <- residuals(fit_ll)

    # Smooth coefficients (and gradients) at the training observations, for alpha
    if (use_ezdat) {
      fit_ll_tz <- llscoef(bw,
                           formula = ll_formula,
                           data = data_np_tz,
                           betas = TRUE,
                           gradients = gradients && n_z_cont > 0L)
      coef_np_tz <- coef(fit_ll_tz)
    } else {
      fit_ll_tz <- fit_ll
      coef_np_tz <- coef_np
    }

    # llscoef() returns the gradients in a matrix with the order
    # d beta_1 / d z_1, ..., d beta_1 / d z_q, d beta_2 / d z_1, ...
    # rearrange them in one matrix for each continuous z variable, as lcscoefgrad() does
    if (gradients && n_z_cont > 0L) {
      grad_tb <- .spsc_ll_grad_list(fit_ll$grad, length(coef_names), z_cont_names)
      grad_tb_tz <- if (use_ezdat) {
        .spsc_ll_grad_list(fit_ll_tz$grad, length(coef_names), z_cont_names)
      } else {
        grad_tb
      }
    }
  }
  
  # Label the first smooth coefficient as the intercept theta(z), and the others
  # with the names of the x variables
  colnames(coef_np) <- coef_names
  colnames(coef_np_tz) <- coef_names

  # Name the gradients
  if (!is.null(grad_tb)) {
    grad_tb <- .spsc_name_grad_list(grad_tb, coef_names)
    grad_tb_tz <- .spsc_name_grad_list(grad_tb_tz, coef_names)
  }
  
  if (verbose) {
    # Local-linear gradients are estimated together with the coefficients
    if (spscreg == "ll") {
      cat("\nSmooth Coefficients Summary:\n")
      print(summary(coef_np))
    }

    # gradients of theta and beta, one z variable at a time
    for (zn in names(grad_tb)) {
      cat("\nGradients of the smooth coefficients with respect to ", zn,
          ":\n", sep = "")
      print(summary(grad_tb[[zn]]))
    }
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
  # columns of the Z matrix that each inefficiency determinant gives: one for a
  # continuous variable, the binary variables for a categorical one
  ineff_columns <- list()
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
        binary_name <- paste0(var_name, "_", levels_cat[k])
        Z_ml_list[[binary_name]] <- binary_var
        ineff_columns[[var_name]] <- c(ineff_columns[[var_name]], binary_name)
      }
    } else {
      Z_ml_list[[var_name]] <- var_data
      ineff_columns[[var_name]] <- var_name
    }
  }
  
  if (length(Z_ml_list) > 0) {
    Z_ml <- do.call(cbind, Z_ml_list)
  } else {
    Z_ml <- matrix(0, nrow = n_obs, ncol = 0)
  }
  Z_ml <- as.matrix(Z_ml)   # ensure matrix type
  n_z_ml <- ncol(Z_ml)
  
  # Prepare a note to report if categorical variables were converted
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
      message(
        "Warning: Optimization did not converge. Convergence code: ",
        optim_result$convergence
      )
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
      message(
        "Warning: Optimization did not converge. Convergence code: ",
        optim_result$convergence
      )
    }
    
    param_mle <- optim_result$par
    loglik_value <- -optim_result$value
    optim_obj <- optim_result
    
    } else if (optim_method == "nlm") {
    nlm_result <- nlm(f = function(param) llk_het(param, resid_mle, Z_ml, n_z_ml),
                      p = initial_params,
                      hessian = TRUE,
                      iterlim = if (!is.null(optim_control$maxit)) optim_control$maxit else 100)
    
    # Report only nlm termination codes that indicate a convergence problem.
    # Codes 1 and 2 indicate that the current iterate is probably a solution.
    if (nlm_result$code == 3L) {
      message("Warning: nlm stopped with code 3.")
    } else if (nlm_result$code == 4L) {
      message("Warning: nlm stopped with code 4.")
    } else if (nlm_result$code == 5L) {
      message("Warning: nlm stopped with code 5.")
    }
    
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
    message(
      "Warning: Hessian is computationally singular. Using regularized Hessian ",
      "(adding 1e-6 to diagonal) for standard errors."
    )
    hessian_reg <- hessian + diag(1e-6, nrow(hessian))
    # In case regularization still fails, return NA
    tryCatch({
        sqrt(diag(solve(hessian_reg)))
    }, error = function(e2) {
        message("Warning: Regularized Hessian also singular. Standard errors set to NA.")
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

  # Gradient of alpha(z) = theta(z) + E(u|z), E(u|z) = sqrt(2/pi) exp(delta_0 +
  # delta_1'z), with respect to each continuous z variable (training data):
  #   d alpha / d z_l = d theta / d z_l + E(u|z) * delta_1l,
  # where delta_1l is the coefficient of that variable in the inefficiency equation.
  alpha_grad <- NULL

  if (gradients) {
    alpha_grad <- matrix(NA_real_, nrow = n_obs, ncol = n_z_cont)
    # with no continuous z variable alpha_grad has no columns
    if (n_z_cont > 0L) {
      colnames(alpha_grad) <- paste0("d(Intercept (alpha))/d(", z_cont_names, ")")
    }

    for (l in seq_len(n_z_cont)) {
      zn <- z_cont_names[l]
      delta_l <- if (n_z_ml > 0L && zn %in% colnames(Z_ml)) {
        delta_mle[1L + match(zn, colnames(Z_ml))]
      } else {
        0
      }
      alpha_grad[, l] <- grad_tb_tz[[zn]][, 1L] + as.vector(Eu_mle) * delta_l
    }

    if (verbose && n_z_cont > 0L) {
      cat("\nIntercept (alpha) Gradients Summary:\n")
      print(summary(alpha_grad))
    }
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
  
  # JLMS estimator E(u|epsilon) and its first-order analytic derivatives
  # (Kumbhakar and Sun, 2013, Economics Letters 120, 249-253, Equation 9).
  #
  # u ~ Truncated Normal(mu, sigma_u^2) and v ~ Normal(0, sigma_v^2), with
  #   mu      = c0 + delta' z_delta
  #   sigma_u = exp(c1 + gamma' z_gamma)
  #   sigma_v = exp(c2 + rho' z_rho)
  # l = 0 returns E(u|epsilon); l = j returns the derivative of E(u|epsilon)
  # with respect to the j-th variable. delta, gamma and rho hold the slope
  # coefficients of mu, sigma_u and sigma_v, and any of them
  # may be NA when the corresponding function does not depend on z.
  # In the semiparametric smooth-coefficient model u is half-normal (mu = 0,
  # delta = NA) and sigma_v is constant (rho = NA), so only gamma, the
  # coefficients on z in sigma_u (delta_1 in the notation of Sun and 
  # Kumbhakar, 2013, Economics Letters, 120, 305–309), is used.
  JLMS_inefficiency <- function(l = 0, epsilon, sigma_v, sigma_u, mu = 0,
                                delta = NA, gamma = NA, rho = NA) {
    sigma2_v <- sigma_v^2
    sigma2_u <- sigma_u^2
    sigma2 <- sigma2_v + sigma2_u
    sigma <- sqrt(sigma2)
    sigma_star <- sigma_v * sigma_u / sigma
    mu_tilde <- (-sigma2_u * epsilon + mu * sigma2_v) / sigma2
    m <- mu_tilde / sigma_star
    
    cdf <- pnorm(m)
    pdf <- dnorm(m)
    g <- pdf / cdf
    
    if (l == 0) {
      
      # E(u|epsilon)
      as.vector(mu_tilde + sigma_star * g)
      
    } else {
      
      m2 <- m^2
      g2 <- g^2
      sigma2_star <- sigma_star^2
      
      # channel through mu
      comp1 <- delta[l] * (sigma2_v * (1 - g * m - g2) / sigma2)
      
      # channel through sigma_u
      comp2 <- gamma[l] * (sigma2_v * sigma_star * (g * (1 + m2) + m * g2) -
                             2 * sigma2_star * (epsilon + mu) * (1 - g2 - m * g)) / sigma2
      
      # channel through sigma_v
      comp3 <- rho[l] * (sigma2_u * sigma_star * (g * (1 + m2) + m * g2) +
                           2 * sigma2_star * (epsilon + mu) * (1 - g2 - m * g)) / sigma2
      
      # derivative of E(u|epsilon) with respect to the l-th variable
      as.vector(rowSums(cbind(comp1, comp2, comp3), na.rm = TRUE))
      
    }
  }
  
  inefficiency <- JLMS_inefficiency(l = 0,
                                    epsilon = epsilon_star_mle,
                                    sigma_v = sigma_v_mle,
                                    sigma_u = s1_mle,
                                    mu = 0)
  
  # First-order analytic derivatives of the JLMS inefficiency estimates with
  # respect to each column of Z_ml: one column per inefficiency variable. They
  # are computed only if requested ('gradients = TRUE'), and are NULL if there
  # is no inefficiency variable
  inefficiency_grad <- NULL

  if (gradients && n_z_ml > 0) {
    inefficiency_grad <- matrix(NA_real_, nrow = n_obs, ncol = n_z_ml,
                                dimnames = list(NULL, colnames(Z_ml)))
    for (l in seq_len(n_z_ml)) {
      inefficiency_grad[, l] <- JLMS_inefficiency(l = l,
                                                  epsilon = epsilon_star_mle,
                                                  sigma_v = sigma_v_mle,
                                                  sigma_u = s1_mle,
                                                  mu = 0,
                                                  delta = NA,
                                                  gamma = delta_mle[-1],
                                                  rho = NA)
    }
  }
  
  if (verbose) {
    cat("Inefficiency Summary:\n")
    cat("Mean Inefficiency:", round(mean(inefficiency), 6), "\n")
    cat("Std Dev Inefficiency:", round(sd(inefficiency), 6), "\n")
    cat("Min Inefficiency:", round(min(inefficiency), 6), "\n")
    cat("Max Inefficiency:", round(max(inefficiency), 6), "\n")
    if (gradients && !is.null(inefficiency_grad)) {
      cat("\nInefficiency Gradients Summary:\n")
      print(summary(inefficiency_grad))
    }
  }
  
  # ========== COMPILE RESULTS ==========
  
  # Estimation sample (the model frame without its attributes). It is kept in
  # the returned object so that methods such as plot() can re-evaluate the
  # smooth coefficients on a grid of points through 'ezdat'.
  model_frame <- as.data.frame(clean_data)
  attr(model_frame, "terms") <- NULL
  attr(model_frame, "na.action") <- NULL

  results <- list(
    # Step 1 Results
    smooth_coefficients = list(
      theta_beta = coef_np,
      alpha = alpha,
      residuals = resid_np,
      bandwidth = bw,
      bandwidth_user_supplied = bandwidth_user_supplied,
      spscreg = spscreg,
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
      inefficiency_columns = ineff_columns,
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
      bandwidth_user_supplied = bandwidth_user_supplied,
      spscreg = spscreg,
      gradients = gradients,
      z_cont_names = z_cont_names
    ),
    
    # Data used information
    data_info = data_info,
    formula_parsed = formula_parsed,
    
    # Estimation sample (used by plot.spscsfa)
    model_frame = model_frame
  )
  
  # Gradients are returned only if requested: smooth_coefficients_grad right
  # after smooth_coefficients, and inefficiency_grad right after inefficiency
  if (gradients) {
    if (is.null(grad_tb)) {
      grad_tb <- stats::setNames(list(), character(0))
    }

    results <- append(
      results,
      list(smooth_coefficients_grad = list(theta_beta_grad = grad_tb,
                                           alpha_grad = alpha_grad)),
      after = match("smooth_coefficients", names(results))
    )

    results <- append(
      results,
      list(inefficiency_grad = inefficiency_grad),
      after = match("inefficiency", names(results))
    )
  }

  class(results) <- "spscsfa"
  
  if (verbose) cat("\n========== ANALYSIS COMPLETE ==========\n")
  
  return(results)
}


#' @rdname spscsfa
#' @usage NULL
#' @export
spsc_sfa <- spscsfa


## ---------------------------------------------------------------------------
## Internal helpers (not exported, not documented)
## ---------------------------------------------------------------------------

## Regression type of a bandwidth object of npscoefbw(): "lc" (local-constant,
## the default of npscoefbw(), "ll" (local-linear) or "lp" (local-polynomial).
.spsc_np_regtype <- function(bw_obj) {
  regtype <- bw_obj[["regtype"]]

  if (is.null(regtype)) "lc" else as.character(regtype)[1L]
}

## Check the bandwidth type, scaling, and kernel settings recorded by
## npscoefbw(). Missing components are treated as the corresponding defaults.
.spsc_validate_np_bw_settings <- function(bw_obj, arg, estimator) {
  settings <- list(
    list("type",      "bwtype",    "fixed"),
    list("scaling",   "bwscaling", FALSE),
    list("ckertype",  "ckertype",  "gaussian"),
    list("ckerorder", "ckerorder", 2),
    list("ckerbound", "ckerbound", "none"),
    list("ukertype",  "ukertype",  "aitchisonaitken"),
    list("okertype",  "okertype",  "liracine")
  )

  for (set in settings) {
    value <- bw_obj[[set[[1L]]]]

    if (!is.null(value) &&
        !identical(as.character(value)[1L], as.character(set[[3L]]))) {
      stop("'", arg, "' is an object returned by npscoefbw() with ",
           set[[2L]], " = ", deparse(value)[1L], ": the ", estimator,
           " estimator needs ", set[[2L]], " = ",
           deparse(set[[3L]])[1L],
           " (it uses fixed bandwidths and the default kernels of ",
           "npscoefbw())", call. = FALSE)
    }
  }

  invisible(NULL)
}

## An npscoefbw() object supplies LC bandwidths only when it records
## regtype = "lc" (or the LC default), fixed bandwidths, and the default kernels 
## used by spscsfa(spscreg = "lc"). Returns its numeric bandwidth vector.
.lcscoef_np_bw <- function(bw_obj, arg) {
  regtype <- .spsc_np_regtype(bw_obj)

  if (!identical(regtype, "lc")) {
    stop("'", arg, "' is an object returned by npscoefbw() with regtype = \"",
         regtype, "\"",
         if (identical(regtype, "ll")) {
           ", whose bandwidths are for the local-linear estimator"
         } else if (identical(regtype, "lp")) {
           ", whose bandwidths are for the local-polynomial estimator"
         },
         ": the local-constant estimator needs regtype = \"lc\"",
         call. = FALSE)
  }

  .spsc_validate_np_bw_settings(bw_obj, arg, "local-constant")

  bw <- bw_obj[["bw"]]
  if (!is.numeric(bw)) {
    stop("'", arg, "' does not hold numeric bandwidths", call. = FALSE)
  }

  bw
}

## Bandwidths given by the user as a numeric vector with one value for each z
## variable, in the order of the formula. 'bandwidth' can be a numeric vector,
## named or not, or an object with a component 'bw' (the bandwidth object of
## llscoefbw(), or of npscoefbw() with regtype = "ll", for the local-linear
## estimator). A named vector is matched by name to the z variables, either as
## they are written in the formula ('z_names') or as np names them
## ('np_z_names'); otherwise it is taken in the order of the formula.
.spsc_bw_vector <- function(bandwidth, z_names, np_z_names) {

  bw <- bandwidth

  if (is.list(bw) && !is.null(bw[["bw"]])) {
    bw <- bw[["bw"]]
  }

  if (!is.numeric(bw)) {
    stop("'bandwidth' must be a numeric vector of bandwidths or a bandwidth ",
         "object", call. = FALSE)
  }

  if (length(bw) != length(z_names)) {
    stop("'bandwidth' has ", length(bw), " value(s) but there are ",
         length(z_names), " z variable(s): ", paste(z_names, collapse = ", "),
         call. = FALSE)
  }

  nm <- names(bw)

  if (!is.null(nm) && !anyNA(nm) && !anyDuplicated(nm)) {
    if (setequal(nm, z_names)) {
      bw <- bw[z_names]
    } else if (setequal(nm, np_z_names)) {
      bw <- bw[np_z_names]
    }
  }

  as.numeric(bw)
}

## The gradients of llscoef(), a matrix with one row for each evaluation point
## and the columns d beta_1 / d z_1, ..., d beta_1 / d z_q, d beta_2 / d z_1,
## ..., d beta_2 / d z_q, ..., as a list with one matrix for each continuous 
## z variable: the columns of the matrix for z_l are d beta_1 / d z_l, ..., 
## d beta_p / d z_l, the layout of lcscoefgrad().
.spsc_ll_grad_list <- function(G, n_coef, z_cont_names) {

  q <- length(z_cont_names)

  if (is.null(G) || ncol(G) != n_coef * q) {
    stop("internal error: the gradients of llscoef() do not have the expected ",
         "dimensions", call. = FALSE)
  }

  out <- lapply(seq_len(q), function(l) {
    G[, (seq_len(n_coef) - 1L) * q + l, drop = FALSE]
  })

  names(out) <- z_cont_names
  out
}

## Column names d(coefficient)/d(z variable); the coefficients are named as  
## in the smooth coefficient matrix: "Intercept (theta)", then the x variables.
.spsc_name_grad_list <- function(grad_list, coef_names) {

  for (zn in names(grad_list)) {
    if (ncol(grad_list[[zn]]) != length(coef_names)) {
      stop("internal error: the gradients do not have one column for each ",
           "smooth coefficient", call. = FALSE)
    }
    dimnames(grad_list[[zn]]) <- list(
      NULL,
      paste0("d(", coef_names, ")/d(", zn, ")")
    )
  }

  grad_list
}
