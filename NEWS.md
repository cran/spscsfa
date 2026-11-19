# spscsfa 0.2.0

All changes since version 0.1.0.

## New arguments of `spscsfa()`

* `spscreg` selects the estimator of the smooth coefficients in the first
  step: `"lc"` (default) for the local-constant estimator, and `"ll"` for the
  local-linear estimator.
  * With `bandwidth = NULL` (default) the bandwidths are selected by
    least-squares cross-validation (LSCV) for the chosen estimator.
  * The residuals of the chosen estimator are used in the second step, the
    maximum likelihood estimation.
  * `smooth_coefficients$theta_beta` and `smooth_coefficients$alpha` are local
    constant or local-linear estimates, as chosen. The choice is recorded in
    `smooth_coefficients$spscreg` and `model_summary$spscreg`.
* `gradients` (default `FALSE`) asks for gradients. If `FALSE`, gradient estimates 
  are not returned. The returned `model_summary$gradients` is set to 
  `FALSE`. If `TRUE`, `spscsfa()` returns two more components:
  * `smooth_coefficients_grad`, a list with two elements:
    * `theta_beta_grad`, a list with one matrix for each continuous `z`
      variable (named after it), with the gradients of the intercept
      `theta(z)` and of the slope coefficients `beta(z)` with respect to that
      variable, at the evaluation points (`ezdat` if supplied, the training data
      otherwise).
    * `alpha_grad`, a matrix at the training data with one column for each
      continuous `z` variable, named `d(Intercept (alpha))/d(z)`.
  * `inefficiency_grad`, the first-order analytic derivatives of the JLMS inefficiency 
    estimates with respect to the inefficiency determinants. 
    `NULL` if there is no determinant.
  * With `spscreg = "lc"` the gradients are those of the local-constant
    estimator, computed by `lcscoefgrad()`; with  `spscreg = "ll"` 
    they are those of the local-linear estimator, from `llscoef()`.
  * A warning is given if there is no continuous `z` variable.
  * With `verbose = TRUE`, summaries of the gradients are shown as they are 
    estimated.
* `bandwidth` accepts a numeric vector with one bandwidth for each `z`
  variable, or a bandwidth object. A named vector, such as the one returned by
  `llscoefbw(LSCV = FALSE)`, is matched to the `z` variables by
  name; otherwise the bandwidths are taken in the order of the formula.

## New functions

* `llscoefbw()` selects the LSCV bandwidths of the local-linear semiparametric
  smooth-coefficient (SPSC) estimator.
  * The bandwidth vector is named after the `z` variables of the formula.
  * It returns an object of class `"llsmoothcoefficientbw"` when the bandwidths 
    are computed, and the named vector with `LSCV = FALSE`.
  * It utilizes a search algorithm that combines Nelder-Mead search with a coarse
    grid scan.
  * The new argument `max.search` (default 2) is the maximum number of searches;
    `max.search = 1` gives the single search from the starting values, without the
    scan. `n_feval` counts the evaluations of the scan and of every search.
    `spscsfa(spscreg = "ll")` uses the default.
* `llscoef()` is the local-linear SPSC estimator, which works with `llscoefbw()`.
  * The bandwidths can be a bandwidth object, a formula (the bandwidths are
    then selected first) or a numeric vector (a formula is then needed).
  * It estimates the smooth coefficients at the training or evaluation data, 
    and returns residuals and, with `gradients = TRUE`, the gradients of the 
    coefficients with respect to the continuous `z` variables.
  * `coef()`, `fitted()`, `residuals()` and `print()` methods for the class
    `"llsmoothcoefficient"`.
* `lcscoefgrad()` computes the analytic gradients of the local-constant SPSC
  estimator with respect to the continuous `z` variables. 
  It returns a list with one matrix for each continuous `z` variable, with
  the gradients of the intercept and of the slope coefficients. 
  Evaluation points can be given in `ezdat`.
* `grad()` extracts the gradients of a model estimated with
  `gradients = TRUE`: a list with `theta_beta_grad`, `alpha_grad` and
  `inefficiency_grad`. With `zvar` it extracts those with respect to some
  variables, for example `grad(fit, zvar = "z1")`. If the variable is not a
  continuous `z` variable (a factor, say), only `inefficiency_grad` for its
  binary variables can be extracted.
* `spsc_sfa()` and `smooth_coef()` become aliases of `spscsfa()` and
    `smoothcoef()`, respectively.
* Fits (and summaries) saved with the class `"spsc_sfa"` dispatch to the
    new `plot`, `print`, and `summary` methods, and are accepted by `te()`,
    `inefficiency()` and `smoothcoef()`.

## New summary and plot methods

* The methods are now `summary.spscsfa()`, `print.spscsfa()`,
  `print.summary.spscsfa()` and `plot.spscsfa()`.
* `summary()` now returns an object of class `"summary.spscsfa"` rather
  than printing the results itself. Printing is handled by the new
  `print.summary.spscsfa()` method.
* The summary object is a list (model specification, step 1, step 2,
  technical efficiency and inefficiency), whose components include
  `coefficients` and `variance`, matrices with the columns `Estimate`,
  `Std. Error`, `z value` and `Pr(>|z|)`, `logLik`, `optim_method`, `te`,
  `inefficiency`, `spscreg` and, if gradients were estimated, `gradients`.
* The inefficiency coefficients and the variance parameter are printed by
  `printCoefmat()`, so `print()` takes `digits`, `signif.stars`, `signif.legend` 
  and other arguments of `print()`. `print()` of a fit shows the
  same as the summary and returns the fit invisibly.
* If the model was estimated with `gradients = TRUE`, `print()` and
  `summary()` also report the summary statistics of the gradients, whatever `verbose` is. 
  `summary()` stores them as `gradients`.
* `plot()` is a new method for fitted models. It plots each smooth coefficient
  (`theta(z)` and `beta(z)`) against one continuous `z` variable at a time, over
  a grid from `min(z)` to `max(z)`, evaluated through `ezdat`. The other
  continuous `z` variables are held at their averages, and discrete ones at
  their most frequent category and not plotted. It uses the estimator of the 
  fitted model (`spscreg`, local-constant or local-linear), re-uses its bandwidth, and 
  returns the plotted values invisibly. `plot()` works with `factor()`, `ordered()` and 
  transformed terms.

## Other changes

* Terms such as `factor(z2)`, `ordered(z2)` or `log(x1)` now work when their variables are
  only in `data`.
* New help pages: `grad`, `llscoefbw`, `llscoef`, `lcscoefgrad` and the
  methods `plot.spscsfa`, `print.summary.spscsfa`. The help pages of
  `spscsfa`, `summary`, `print` and `smoothcoef` are updated, and the old
  names are documented as aliases.
