#' Extract technical-efficiency scores
#'
#' @param x An object of class `spsc_sfa`.
#' @return A numeric vector of technical-efficiency estimates.
#' @export
te <- function(x) {
  if (!inherits(x, "spsc_sfa")) stop("x must be of class 'spsc_sfa'")
  x$technical_efficiency$TE
}

#' Extract JLMS inefficiency estimates
#'
#' @param x An object of class `spsc_sfa`.
#' @return A numeric vector of JLMS conditional inefficiency estimates.
#' @export
inefficiency <- function(x) {
  if (!inherits(x, "spsc_sfa")) stop("x must be of class 'spsc_sfa'")
  x$inefficiency$u_JLMS
}

#' Extract smooth coefficients
#'
#' @param x An object of class `spsc_sfa`.
#' @return The estimated smooth-coefficient matrix.
#' @export
smooth_coef <- function(x) {
  if (!inherits(x, "spsc_sfa")) stop("x must be of class 'spsc_sfa'")
  list(
    theta_beta = x$smooth_coefficients$theta_beta,
    alpha = x$smooth_coefficients$alpha
  )
}
