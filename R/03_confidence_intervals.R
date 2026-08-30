#' Run an expression, muffling bbmle's profiling diagnostics and
#' identifying which parameter(s) they refer to
#'
#' \code{bbmle::profile()}/\code{bbmle::confint(..., method = "spline")}
#' can emit several low-level diagnostic warnings when a parameter's
#' profile is poorly behaved (flat, hits the iteration limit, or is
#' non-monotonic so the spline interpolation falls back to linear) - e.g.
#' \code{"stepsize effectively zero/flat profile (Linf)"} or
#' \code{"non-monotonic profile (Linf): reverting from spline to linear
#' approximation"}. Individually these are cryptic and easy to miss;
#' this helper muffles them and instead returns which parameter names
#' (matched by the literal \code{"(name)"} bbmle appends) were flagged, so
#' the caller can issue one clear, consolidated warning instead.
#'
#' @param expr Expression to evaluate (e.g. a \code{bbmle::confint(...)} call).
#' @param parameters Character vector of parameter names to look for.
#' @return List with \code{result} (the value of \code{expr}) and
#'   \code{flaky} (character vector of parameter names implicated in a
#'   muffled warning, possibly empty).
#' @keywords internal
.with_profile_warnings <- function(expr, parameters) {
  flaky <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      hit <- parameters[vapply(parameters, function(p) grepl(paste0("(", p, ")"), msg, fixed = TRUE),
                                logical(1))]
      if (length(hit)) flaky <<- union(flaky, hit)
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, flaky = flaky)
}

#' Confidence intervals via likelihood profile
#'
#' Computes confidence intervals for the parameters of an already-fitted
#' growth model, using the likelihood profile method
#' (\code{stats::confint.nls} when \code{method = "ls"}, or
#' \code{bbmle::confint} when \code{method = "mle"}). Likelihood profiles
#' tend to produce asymmetric intervals, which are more realistic than
#' those based on the quadratic (Wald) approximation, especially with
#' small samples or correlated parameters.
#'
#' Sometimes the profile for a specific parameter is poorly behaved (flat
#' or non-monotonic - typically because that parameter is weakly
#' identified by the data, e.g. \code{Linf} when no individuals are old
#' enough to approach the asymptote). \code{bbmle} normally reports this
#' as several cryptic low-level warnings naming the affected parameter
#' (e.g. \code{"flat profile (Linf)"}); \code{profile_ci()} catches these,
#' still returns \pkg{bbmle}'s best-effort CI, but replaces the raw
#' warnings with a single clear one and flags the affected row(s) via the
#' \code{reliable} column, so a suspect interval is easy to spot rather
#' than silently trusted.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param level Confidence level (default 0.95).
#' @param profile_method Only for \code{method = "mle"}: \code{"spline"}
#'   (full likelihood profile, recommended) or \code{"quad"} (quadratic /
#'   Wald approximation, faster).
#' @return Data frame with columns \code{parameter}, \code{estimate},
#'   \code{ci_lower}, \code{ci_upper}, \code{level}, \code{method}, and
#'   \code{reliable} (\code{FALSE} for any parameter whose profile hit a
#'   numerical problem - see Details), or \code{NULL} (with a warning) if
#'   the profile could not be computed at all (in which case
#'   \code{\link{bootstrap_ci}} is recommended instead).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy",
#'                    method = "mle")
#' profile_ci(fit, level = 0.95)
#' @export
profile_ci <- function(fit, level = 0.95, profile_method = c("spline", "quad")) {
  stopifnot(inherits(fit, "growth_fit"))
  profile_method <- match.arg(profile_method)
  parameters <- names(fit$coefficients)
  flaky_params <- character(0)

  if (fit$method == "ls") {
    ci <- try(suppressMessages(stats::confint(fit$fit_obj, level = level)), silent = TRUE)
    method_label <- "likelihood profile (nls)"
  } else {
    wrapped <- .with_profile_warnings(
      try(suppressMessages(bbmle::confint(fit$fit_obj, level = level, method = profile_method)),
          silent = TRUE),
      parameters
    )
    ci <- wrapped$result
    flaky_params <- wrapped$flaky
    method_label <- if (profile_method == "spline") {
      "likelihood profile (mle2)"
    } else {
      "quadratic / Wald approximation (mle2)"
    }
  }

  if (inherits(ci, "try-error")) {
    warning("Could not compute the likelihood-profile CI for model '",
            fit$model, "': ", conditionMessage(attr(ci, "condition")),
            ". Try bootstrap_ci() instead.")
    return(NULL)
  }
  if (!is.matrix(ci)) {
    # bbmle::confint(..., method = "spline") returns the *refit* mle2 object
    # (instead of a confint matrix) when profiling finds a better optimum,
    # i.e. when the original fit had not fully converged.
    warning("Could not compute the likelihood-profile CI for model '", fit$model,
            "': profiling found a better optimum, meaning the original fit had ",
            "not fully converged. Refit with different starting values (or ",
            "method = \"ls\"), or use bootstrap_ci() instead.")
    return(NULL)
  }
  if (is.null(dim(ci))) ci <- matrix(ci, nrow = 1, dimnames = list(parameters, names(ci)))
  ci <- ci[parameters, , drop = FALSE]

  if (length(flaky_params) > 0) {
    plural <- length(flaky_params) > 1
    msg <- sprintf(
      "profile_ci(): the likelihood profile for %s had numerical problems (flat and/or non-monotonic profile) while fitting model '%s' - treat %s confidence interval%s as unreliable. Consider bootstrap_ci() as a cross-check, or refit with different starting values.",
      paste(sprintf("'%s'", flaky_params), collapse = " and "),
      fit$model,
      if (plural) "their" else "its",
      if (plural) "s" else ""
    )
    warning(msg, call. = FALSE)
  }

  data.frame(
    parameter = parameters,
    estimate = as.numeric(fit$coefficients[parameters]),
    ci_lower = as.numeric(ci[, 1]),
    ci_upper = as.numeric(ci[, 2]),
    level = level,
    method = method_label,
    reliable = !(parameters %in% flaky_params),
    row.names = NULL
  )
}

#' Compute the raw likelihood-profile trace of a growth_fit object
#'
#' Low-level helper that returns the raw profile trace (parameter value
#' vs. signed profile statistic) used internally by \code{\link{plot_profile}}.
#' Most users will call \code{\link{plot_profile}} directly instead.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param which Optional character vector restricting the profile to a
#'   subset of parameters. By default, all model parameters (excluding
#'   \code{log_sigma} for MLE fits).
#' @param level Confidence level used for the critical value and (if
#'   \code{ci} is not supplied) the profile confidence interval (default 0.95).
#' @param ci Optional: an already-computed \code{\link{profile_ci}} result
#'   for this exact \code{fit}. When supplied, it is used as-is for the
#'   vertical CI reference lines instead of being recomputed internally -
#'   useful to guarantee that what gets plotted is exactly the object the
#'   user already computed and inspected (e.g. with a specific
#'   \code{profile_method}), rather than a silent second call to
#'   \code{profile_ci()} with default arguments. If \code{NULL} (default),
#'   computed internally exactly as before. \code{ci} is validated against
#'   \code{fit}: it must cover all of \code{fit}'s parameters by name, and
#'   (when a \code{profile_ci()} result, which stamps an \code{estimate}
#'   column) its point estimates must match \code{fit$coefficients} - this
#'   second check catches the case where \code{ci} was computed for a
#'   different model that happens to reuse the same parameter names (e.g.
#'   \code{von_bertalanffy} and \code{gompertz} both have \code{Linf},
#'   \code{K}, \code{t0}), which a name-only check would miss.
#' @return List with elements \code{data} (a tidy data frame with columns
#'   \code{parameter}, \code{focal}, \code{tau}, and \code{deviance} —
#'   the chi-squared profile statistic, \code{deviance = tau^2}),
#'   \code{critical_chisq} (the \eqn{\chi^2_1} critical value for the
#'   chosen confidence level, i.e. \code{qchisq(level, df = 1)}),
#'   \code{critical_tau} (its square root, the corresponding threshold on
#'   the signed \eqn{\tau} scale) and \code{ci} (the profile confidence
#'   intervals, or \code{NULL}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy",
#'                    method = "mle")
#' trace <- profile_trace(fit)
#' str(trace$data)
#' trace$critical_chisq
#' @export
profile_trace <- function(fit, which = NULL, level = 0.95, ci = NULL) {
  stopifnot(inherits(fit, "growth_fit"))
  parameters <- names(fit$coefficients)
  if (!is.null(which)) parameters <- intersect(parameters, which)
  flaky_params <- character(0)

  if (!is.null(ci)) {
    stopifnot(is.data.frame(ci), all(c("parameter", "ci_lower", "ci_upper") %in% names(ci)))
    if (!all(parameters %in% ci$parameter)) {
      stop("profile_trace(): 'ci' does not cover all parameters of 'fit' (fit: ",
           paste(parameters, collapse = ", "), "; ci: ",
           paste(unique(ci$parameter), collapse = ", "), "). It looks like 'ci' ",
           "was computed for a different model - pass a profile_ci() result ",
           "computed for this exact 'fit'.")
    }
    # Matching parameter *names* is not enough: two different models can use
    # the same parameter names (e.g. von_bertalanffy and gompertz both have
    # Linf, K, t0). Cross-check the point estimates too, since profile_ci()
    # always stamps its own fit's coefficients into the 'estimate' column -
    # a mismatch there means 'ci' almost certainly comes from a different fit.
    if ("estimate" %in% names(ci)) {
      ci_est <- ci$estimate[match(parameters, ci$parameter)]
      fit_est <- as.numeric(fit$coefficients[parameters])
      if (!isTRUE(all.equal(ci_est, fit_est, tolerance = 1e-6, check.attributes = FALSE))) {
        stop("profile_trace(): 'ci' shares parameter names with 'fit' (",
             paste(parameters, collapse = ", "), ") but its point estimates ",
             "do not match 'fit$coefficients' (fit: ",
             paste(sprintf("%.6g", fit_est), collapse = ", "), "; ci: ",
             paste(sprintf("%.6g", ci_est), collapse = ", "), "). This looks like ",
             "'ci' was computed for a different model fit that happens to reuse ",
             "the same parameter names - pass a profile_ci() result computed ",
             "for this exact 'fit'.")
      }
    }
    if ("level" %in% names(ci) && !isTRUE(all.equal(unique(ci$level), level))) {
      warning("profile_trace(): 'ci' was computed at level = ",
               paste(unique(ci$level), collapse = ", "), " but 'level' = ", level,
               " was requested here; the chi-squared reference line uses 'level', ",
               "while the plotted CI bounds come from 'ci' - the two may not line up.")
    }
  }

  if (fit$method == "ls") {
    prof <- try(stats::profile(fit$fit_obj), silent = TRUE)
    if (inherits(prof, "try-error")) {
      warning("Could not compute the likelihood profile for model '", fit$model,
               "': ", conditionMessage(attr(prof, "condition")))
      return(NULL)
    }
    rows <- lapply(parameters, function(pn) {
      pdf <- prof[[pn]]
      data.frame(parameter = pn, focal = pdf$par.vals[, pn], tau = pdf$tau)
    })
    trace_df <- do.call(rbind, rows)
  } else {
    wrapped <- .with_profile_warnings(try(bbmle::profile(fit$fit_obj), silent = TRUE), parameters)
    prof <- wrapped$result
    flaky_params <- wrapped$flaky
    if (inherits(prof, "try-error")) {
      warning("Could not compute the likelihood profile for model '", fit$model,
               "': ", conditionMessage(attr(prof, "condition")))
      return(NULL)
    }
    if (!inherits(prof, "profile.mle2")) {
      # bbmle::profile() returns the *refit* mle2 object (instead of a
      # profile.mle2 object) when it finds a better optimum, i.e. when the
      # original fit had not fully converged.
      warning("Could not compute the likelihood profile for model '", fit$model,
              "': profiling found a better optimum, meaning the original fit ",
              "had not fully converged. Refit with different starting values ",
              "(or method = \"ls\"), or use bootstrap_ci() instead.")
      return(NULL)
    }
    raw <- as.data.frame(prof)
    raw <- raw[raw$param %in% parameters, , drop = FALSE]
    trace_df <- data.frame(parameter = raw$param, focal = raw$focal, tau = raw$z)
  }

  # Profile deviance, D(theta) = 2 * (logLik_max - logLik(theta)) = tau^2,
  # asymptotically chi-squared with 1 df (Wilks' theorem); this is the
  # standard chi-squared profile-likelihood statistic and reference
  # distribution, used uniformly for both "ls" and "mle" fits.
  trace_df$deviance <- trace_df$tau^2
  critical_chisq <- stats::qchisq(level, df = 1)
  critical_tau <- sqrt(critical_chisq)

  if (is.null(ci)) {
    # suppressWarnings here: profile_ci() would otherwise re-run its own
    # profiling and (redundantly) re-warn about the same flaky parameter(s)
    # already detected and reported for this trace just above.
    ci <- try(suppressWarnings(profile_ci(fit, level = level)), silent = TRUE)
    if (inherits(ci, "try-error") || is.null(ci)) ci <- NULL
  }
  if (!is.null(ci)) {
    # keep only the parameters actually present in the trace, otherwise a
    # geom_vline layer built from 'ci' would introduce extra (empty) facets
    ci <- ci[ci$parameter %in% parameters, , drop = FALSE]
  }

  if (length(flaky_params) > 0) {
    plural <- length(flaky_params) > 1
    msg <- sprintf(
      "profile_trace(): the likelihood profile for %s had numerical problems (flat and/or non-monotonic profile) while fitting model '%s' - treat that region of the profile plot (and %s confidence interval%s) with caution. Consider bootstrap_ci() as a cross-check.",
      paste(sprintf("'%s'", flaky_params), collapse = " and "),
      fit$model,
      if (plural) "their" else "its",
      if (plural) "s" else ""
    )
    warning(msg, call. = FALSE)
  }

  list(data = trace_df, critical_chisq = critical_chisq, critical_tau = critical_tau,
       ci = ci, level = level)
}

#' Extract a growth_fit object's parameter standard errors
#'
#' Thin dispatcher over \code{vcov()}: for \code{method = "ls"} fits,
#' \code{stats::vcov()} already dispatches correctly on the underlying
#' \code{nls} object; for \code{method = "mle"} fits, the S3 method lives
#' in \pkg{bbmle} and must be called explicitly (\pkg{bbmle} is only
#' \emph{imported}, not attached, so the plain generic would not find it).
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @return Named numeric vector of standard errors, one per element of
#'   \code{fit$coefficients}.
#' @keywords internal
.growth_fit_se <- function(fit) {
  v <- if (fit$method == "ls") stats::vcov(fit$fit_obj) else bbmle::vcov(fit$fit_obj)
  parameters <- names(fit$coefficients)
  se <- sqrt(diag(v))[parameters]
  names(se) <- parameters
  se
}

#' Profile log-likelihood of a growth_fit object with two parameters fixed
#'
#' Refits \code{fit}'s model with the two named parameters in
#' \code{fixed_params} held constant, optimising over whatever parameters
#' remain (using exactly the same \code{.fit_core()} engine - and
#' therefore the same LS/MLE algorithm and error distribution - as the
#' original fit), and returns the resulting log-likelihood. If fixing the
#' two parameters leaves none free (a two-parameter model with both
#' parameters fixed), the concentrated normal/lognormal log-likelihood is
#' computed directly from the residuals, with no optimisation needed.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param fixed_params Named list of length 2: the two parameter values to
#'   hold fixed (names must be among \code{names(fit$coefficients)}).
#' @return A single numeric log-likelihood value, or \code{NA_real_} if the
#'   constrained refit failed to converge or produced a non-finite fit.
#' @keywords internal
.profile_loglik_fixed <- function(fit, fixed_params) {
  remaining <- setdiff(names(fit$coefficients), names(fixed_params))
  fn_wrapped <- function(t, ...) {
    free_args <- list(...)
    do.call(fit$fn, c(list(t = t), free_args, fixed_params))
  }
  if (length(remaining) == 0) {
    yhat <- fn_wrapped(t = fit$t)
    if (any(!is.finite(yhat))) return(NA_real_)
    resid <- if (fit$distribution == "normal") fit$y - yhat else log(fit$y) - log(pmax(yhat, 1e-8))
    n <- length(resid)
    rss <- sum(resid^2)
    if (!is.finite(rss) || rss <= 0) return(NA_real_)
    sigma2 <- rss / n
    return(-n / 2 * (log(2 * pi) + log(sigma2) + 1))
  }
  start <- as.list(fit$coefficients[remaining])
  core <- try(
    .fit_core(y = fit$y, predictors = list(t = fit$t), fn = fn_wrapped,
              parameters = remaining, method = fit$method,
              distribution = fit$distribution, start = start),
    silent = TRUE
  )
  if (inherits(core, "try-error")) return(NA_real_)
  ll <- if (core$method == "ls") as.numeric(stats::logLik(core$fit_obj)) else -core$fit_obj@min
  if (!is.finite(ll)) return(NA_real_)
  ll
}

#' Joint (bivariate) confidence region via a 2D likelihood profile
#'
#' For models such as von Bertalanffy, \code{Linf} and \code{K} are
#' typically strongly (negatively) correlated: many combinations of a
#' slightly smaller \code{Linf} with a slightly larger \code{K} fit the
#' data almost as well as the joint optimum. Reporting independent,
#' one-parameter profile intervals for each of the two (as
#' \code{\link{profile_ci}} does, each referenced to
#' \eqn{\chi^2_{1,\,1-\alpha}}) ignores this correlation and can be
#' misleading: the two marginal intervals, read together as a rectangle,
#' generally do \emph{not} match the shape of the actual joint region of
#' statistically-as-good-as-the-optimum parameter pairs.
#'
#' \code{profile_ci_bivariate()} instead profiles the pair \emph{jointly}:
#' it evaluates the profile log-likelihood over a 2D grid of
#' (\code{params[1]}, \code{params[2]}) values (refitting every other
#' parameter at each grid point - see \code{\link{.profile_loglik_fixed}}),
#' and compares the resulting profile deviance,
#' \eqn{D(\theta_1,\theta_2) = 2(\hat\ell - \ell(\theta_1,\theta_2))}, to
#' \eqn{\chi^2_{2,\,1-\alpha}} (\strong{2 degrees of freedom}, one per
#' jointly-profiled parameter - Wilks' theorem), rather than
#' \eqn{\chi^2_1}. The set of grid points with deviance at or below that
#' threshold is the joint confidence region, typically an elongated,
#' tilted ellipse-like shape along the ridge of correlation between the
#' two parameters - visualised with \code{\link{plot_profile_bivariate}}.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param params Character vector of exactly two parameter names (must be
#'   among \code{names(fit$coefficients)}) - the pair to profile jointly.
#' @param level Confidence level (default 0.95).
#' @param grid_n Number of grid points per parameter (the grid is
#'   \code{grid_n x grid_n}, so total refits are \code{grid_n^2}).
#' @param range_mult The grid for each parameter spans
#'   \code{estimate +/- range_mult * SE}, with \code{SE} taken from the
#'   fit's variance-covariance matrix (\code{vcov()}). Widen this if the
#'   confidence region touches the edge of the grid.
#' @return A list of class \code{"profile_ci_bivariate"} with elements:
#'   \code{data} (data frame with one row per grid point: \code{par1},
#'   \code{par2} (the actual parameter values, named after \code{params}),
#'   \code{loglik}, \code{deviance}, and \code{in_region}), \code{params},
#'   \code{level}, \code{critical_chisq} (\eqn{\chi^2_{2,\,1-\alpha}}),
#'   \code{estimate} (the joint MLE/LS point estimate of the pair),
#'   \code{model}, \code{label}, and \code{n_failed} (grid points where the
#'   constrained refit did not converge, excluded from \code{in_region}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy",
#'                    method = "ls")
#' region <- profile_ci_bivariate(fit, params = c("Linf", "K"), level = 0.95,
#'                                 grid_n = 15)
#' region
#' @export
profile_ci_bivariate <- function(fit, params, level = 0.95, grid_n = 25, range_mult = 3.5) {
  stopifnot(inherits(fit, "growth_fit"))
  fit_params <- names(fit$coefficients)
  if (length(params) != 2 || !all(params %in% fit_params)) {
    stop("profile_ci_bivariate(): 'params' must name exactly two parameters of 'fit' (",
         paste(fit_params, collapse = ", "), "), got: ", paste(params, collapse = ", "), ".")
  }
  se <- try(.growth_fit_se(fit), silent = TRUE)
  if (inherits(se, "try-error") || any(!is.finite(se[params]))) {
    stop("profile_ci_bivariate(): could not obtain standard errors for '", fit$model,
         "' (needed to size the grid) - vcov() failed. Try refitting, or check that the ",
         "model converged to a proper optimum.")
  }

  p1 <- params[1]; p2 <- params[2]
  est1 <- fit$coefficients[[p1]]; est2 <- fit$coefficients[[p2]]
  grid1 <- seq(est1 - range_mult * se[[p1]], est1 + range_mult * se[[p1]], length.out = grid_n)
  grid2 <- seq(est2 - range_mult * se[[p2]], est2 + range_mult * se[[p2]], length.out = grid_n)

  ll_full <- fit$logLik
  grid <- expand.grid(par1 = grid1, par2 = grid2)
  grid$loglik <- vapply(seq_len(nrow(grid)), function(i) {
    fixed <- stats::setNames(list(grid$par1[i], grid$par2[i]), c(p1, p2))
    .profile_loglik_fixed(fit, fixed)
  }, numeric(1))
  grid$deviance <- 2 * (ll_full - grid$loglik)

  critical_chisq <- stats::qchisq(level, df = 2)
  grid$in_region <- is.finite(grid$deviance) & grid$deviance <= critical_chisq
  n_failed <- sum(!is.finite(grid$loglik))
  if (n_failed > 0) {
    warning("profile_ci_bivariate(): the constrained refit did not converge at ", n_failed,
             " of ", nrow(grid), " grid points for model '", fit$model,
             "' - those points are excluded from the confidence region. Consider a smaller ",
             "'range_mult' or a coarser 'grid_n'.")
  }

  names(grid)[names(grid) == "par1"] <- p1
  names(grid)[names(grid) == "par2"] <- p2

  out <- list(
    data = grid,
    params = params,
    level = level,
    critical_chisq = critical_chisq,
    estimate = stats::setNames(c(est1, est2), params),
    model = fit$model,
    label = fit$label,
    n_failed = n_failed
  )
  class(out) <- "profile_ci_bivariate"
  out
}

#' Print a profile_ci_bivariate object
#'
#' @param x Object of class \code{"profile_ci_bivariate"}.
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' region <- profile_ci_bivariate(fit, params = c("Linf", "K"), grid_n = 15)
#' print(region)
#' @export
print.profile_ci_bivariate <- function(x, ...) {
  cat(sprintf("Joint (bivariate) likelihood-profile confidence region - %s model\n", x$label))
  cat(sprintf("Parameters: %s = %.5g, %s = %.5g\n",
              x$params[1], x$estimate[[1]], x$params[2], x$estimate[[2]]))
  cat(sprintf("Confidence level: %.0f%% | critical chi-squared (2 df): %.4f\n",
              x$level * 100, x$critical_chisq))
  cat(sprintf("Grid: %d x %d points | in region: %d | failed refits: %d\n",
              length(unique(x$data[[x$params[1]]])), length(unique(x$data[[x$params[2]]])),
              sum(x$data$in_region), x$n_failed))
  in_reg <- x$data[x$data$in_region, , drop = FALSE]
  if (nrow(in_reg) > 0) {
    cat(sprintf("Region range: %s in [%.5g, %.5g] | %s in [%.5g, %.5g]\n",
                x$params[1], min(in_reg[[x$params[1]]]), max(in_reg[[x$params[1]]]),
                x$params[2], min(in_reg[[x$params[2]]]), max(in_reg[[x$params[2]]])))
  }
  invisible(x)
}

#' Confidence intervals via bootstrap
#'
#' Computes confidence intervals (percentiles) for the parameters of a
#' growth model using bootstrap resampling, refitting the model \code{R}
#' times. Two schemes are available: case resampling (complete age-size
#' observations) or residual resampling (keeps the observed ages fixed
#' and resamples the residuals of the original fit).
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param R Number of bootstrap replicates.
#' @param type \code{"cases"} (resampling of age-size pairs, more robust
#'   to model misspecification) or \code{"residuals"} (resampling of
#'   residuals, assumes the model is well specified).
#' @param level Confidence level.
#' @param seed Random seed (\code{NULL} to not set one).
#' @param verbose If \code{TRUE}, reports how many replicates converged.
#' @return List of class \code{"bootstrap_ci_growth"} with: \code{table}
#'   (data frame with estimate and CI per parameter), \code{replicates}
#'   (matrix of bootstrap estimates), \code{R_requested}, \code{R_converged},
#'   \code{level}, and \code{type}.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' boot <- bootstrap_ci(fit, R = 100, type = "cases", seed = 1)
#' boot
#' @export
bootstrap_ci <- function(fit, R = 500, type = c("cases", "residuals"),
                          level = 0.95, seed = NULL, verbose = TRUE) {
  stopifnot(inherits(fit, "growth_fit"))
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  model <- fit$model
  method <- fit$method
  distribution <- fit$distribution
  parameters <- names(fit$coefficients)
  n <- fit$n
  start0 <- as.list(fit$coefficients)

  t <- fit$t; y <- fit$y
  yhat <- do.call(fit$fn, c(list(t = t), start0))
  residuals_ <- y - yhat

  reps <- matrix(NA_real_, nrow = R, ncol = length(parameters),
                  dimnames = list(NULL, parameters))

  for (i in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (type == "cases") {
      tb <- t[idx]; yb <- y[idx]
    } else {
      tb <- t; yb <- yhat + residuals_[idx]
    }
    dfb <- data.frame(age = tb, size = yb)
    fb <- try(fit_growth(dfb, age = "age", size = "size", model = model,
                          method = method, distribution = distribution, start = start0),
              silent = TRUE)
    if (!inherits(fb, "try-error") && !is.null(fb)) {
      reps[i, ] <- fb$coefficients[parameters]
    }
  }

  reps_ok <- reps[stats::complete.cases(reps), , drop = FALSE]
  if (verbose) {
    message(sprintf("Bootstrap (%s, %s): %d of %d replicates converged.",
                     type, toupper(method), nrow(reps_ok), R))
  }
  if (nrow(reps_ok) < 20) {
    warning("Very few bootstrap replicates converged; the CIs may be unreliable.")
  }

  alpha <- 1 - level
  tbl <- data.frame(
    parameter = parameters,
    estimate = as.numeric(fit$coefficients[parameters]),
    ci_lower = apply(reps_ok, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE),
    ci_upper = apply(reps_ok, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE),
    boot_se = apply(reps_ok, 2, stats::sd, na.rm = TRUE),
    row.names = NULL
  )

  out <- list(table = tbl, replicates = reps_ok, R_requested = R,
              R_converged = nrow(reps_ok), level = level, type = type,
              model = model, label = fit$label)
  class(out) <- "bootstrap_ci_growth"
  out
}

#' Print a bootstrap_ci_growth object
#'
#' Prints the bootstrap type and model, the number of converged
#' replicates out of those requested, the confidence level, and the
#' per-parameter confidence interval table.
#'
#' @param x Object of class \code{"bootstrap_ci_growth"} (from
#'   \code{\link{bootstrap_ci}}).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' boot <- bootstrap_ci(fit, R = 100, seed = 1)
#' print(boot)
#' @export
print.bootstrap_ci_growth <- function(x, ...) {
  cat(sprintf("Bootstrap (%s) for model %s\n", x$type, x$label))
  cat(sprintf("Converged replicates: %d / %d | Confidence level: %.0f%%\n",
              x$R_converged, x$R_requested, x$level * 100))
  print(x$table, row.names = FALSE)
  invisible(x)
}
