#' @importFrom stats coef nls predict logLik AIC BIC vcov confint sd approx dnorm lm as.formula
NULL

# ---- Automatic starting values ----------------------------------------------

#' Estimate starting values for a growth model
#'
#' Computes, heuristically, reasonable starting values for the
#' optimisation algorithm from a classic linearisation of the von
#' Bertalanffy model (regression of log(Linf0 - size) on age). These
#' values are reused, with small adjustments, as starting points for the
#' other models in the family. Users can always supply their own values
#' via the \code{start} argument of \code{\link{fit_growth}}.
#'
#' @param t Vector of ages.
#' @param y Vector of sizes (or weights).
#' @param model Model name (see \code{names(GROWTH_MODELS)}).
#' @return Named list of starting values for each model parameter.
#' @examples
#' data(growth_data)
#' initial_values(growth_data$age, growth_data$size, model = "von_bertalanffy")
#' initial_values(growth_data$age, growth_data$size, model = "gompertz")
#' @export
initial_values <- function(t, y, model) {
  stopifnot(model %in% names(GROWTH_MODELS))

  Linf0 <- max(y, na.rm = TRUE) * 1.05
  lin_data <- data.frame(t = t, y = y)
  lin_data <- lin_data[is.finite(lin_data$t) & is.finite(lin_data$y) &
                          lin_data$y < Linf0, , drop = FALSE]
  lin_data$logdiff <- log(Linf0 - lin_data$y)

  K0 <- 0.3
  t0_0 <- 0
  lin_fit <- try(stats::lm(logdiff ~ t, data = lin_data), silent = TRUE)
  if (!inherits(lin_fit, "try-error")) {
    slope <- unname(stats::coef(lin_fit)["t"])
    intercept <- unname(stats::coef(lin_fit)["(Intercept)"])
    if (is.finite(slope) && slope < 0) {
      K0 <- -slope
      t0_0 <- (intercept - log(Linf0)) / K0
    }
  }
  if (!is.finite(K0) || K0 <= 0) K0 <- 0.3
  if (!is.finite(t0_0)) t0_0 <- 0

  switch(model,
    von_bertalanffy  = list(Linf = Linf0, K = K0, t0 = t0_0),
    gompertz         = list(Linf = Linf0, K = K0, t0 = max(t0_0, 0)),
    logistic         = list(Linf = Linf0, K = K0, t0 = max(t0_0, 0)),
    richards         = list(Linf = Linf0, K = K0, t0 = t0_0, b = 1),
    schnute_richards = list(Linf = Linf0, b = 0.9, c = K0, d = 1, g = 1),
    vb_seasonal      = list(Linf = Linf0, K = K0, t0 = t0_0, C = 0.5, ts = 0.5),
    gompertz_laird   = {
      L0_0 <- max(min(y, na.rm = TRUE), 1e-3)
      list(L0 = L0_0, A = K0, k = K0)
    },
    schnute = ,
    schnute_case2 = ,
    schnute_case3 = ,
    schnute_case4 = {
      # y1_0/y2_0: starting values for the sizes at the reference ages
      # t1 = min(t)/t2 = max(t). Since t1 and t2 are themselves elements of
      # t, this is just the (mean) observed size at those exact ages - no
      # interpolation needed. (stats::approx(t, y, xout = t1) would give
      # the same answer here, but warns "collapsing to unique 'x' values"
      # whenever t has tied ages, which is common with real/simulated data;
      # direct subsetting avoids that spurious warning entirely.)
      t1 <- min(t, na.rm = TRUE); t2 <- max(t, na.rm = TRUE)
      y1_0 <- mean(y[t == t1], na.rm = TRUE)
      y2_0 <- mean(y[t == t2], na.rm = TRUE)
      switch(model,
        schnute       = list(a = K0, b = 1, y1 = y1_0, y2 = y2_0),
        schnute_case2 = list(a = K0, y1 = y1_0, y2 = y2_0),
        schnute_case3 = list(b = 1, y1 = y1_0, y2 = y2_0),
        schnute_case4 = list(y1 = y1_0, y2 = y2_0)
      )
    },
    gallucci_quinn = list(omega = Linf0 * K0, K = K0, t0 = t0_0),
    francis_vb = {
      # L1/L2/L3: starting values at the (fixed) reference ages t1 = min(t),
      # t3 = max(t), t2 = midpoint - evaluated from the classic VB
      # linearisation already computed above, rather than interpolated
      # from the raw data (t2 generally does not fall exactly on an
      # observed age).
      t1 <- min(t, na.rm = TRUE); t3 <- max(t, na.rm = TRUE); t2 <- (t1 + t3) / 2
      list(L1 = von_bertalanffy(t1, Linf0, K0, t0_0),
           L2 = von_bertalanffy(t2, Linf0, K0, t0_0),
           L3 = von_bertalanffy(t3, Linf0, K0, t0_0))
    },
    biphasic_growth = {
      T0 <- stats::median(t, na.rm = TRUE)
      list(h = Linf0 * K0, t0 = t0_0, T = T0, Linf = Linf0, K = K0)
    },
    persistence = {
      # Approaches a*t^b as t -> Inf, so a log(y) ~ log(t) regression
      # (excluding t <= 0, where the power function is undefined) gives a
      # reasonable starting point for a/b; c controls how fast the
      # exponent decays away from b at small ages, started at a fraction
      # of the age range.
      pos <- t > 0 & y > 0
      b0 <- 1; a0 <- Linf0 / max(t, na.rm = TRUE)
      lin_pow <- try(stats::lm(log(y[pos]) ~ log(t[pos])), silent = TRUE)
      if (!inherits(lin_pow, "try-error")) {
        slope <- unname(stats::coef(lin_pow)[2])
        intercept <- unname(stats::coef(lin_pow)[1])
        if (is.finite(slope) && slope > 0) b0 <- slope
        if (is.finite(intercept)) a0 <- exp(intercept)
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- 1
      c0 <- max(stats::median(t[pos], na.rm = TRUE) * 0.3, 1e-3)
      list(a = a0, b = b0, c = c0)
    },
    tanaka = {
      # Tanaka's model is algebraically an inverse-hyperbolic-sine curve,
      # L(t) = (1/sqrt(f))*asinh(sqrt(f/a)*(t-c)) + D; anchor the curve to
      # pass through the median observed (age, size) point using generic
      # scale parameters derived from the VB linearisation above.
      c0 <- t0_0
      f0 <- max(K0, 0.05)
      a0 <- max(stats::var(t, na.rm = TRUE), 1)
      t_med <- stats::median(t, na.rm = TRUE)
      y_med <- stats::median(y, na.rm = TRUE)
      d0 <- y_med - (1 / sqrt(f0)) * asinh(sqrt(f0 / a0) * (t_med - c0))
      list(a = a0, c = c0, d = d0, f = f0)
    },
    linear = {
      # L(t) = a + b*t is itself linear, so the ordinary least-squares fit
      # of y on t is already the exact LS solution - an ideal starting point.
      a0 <- mean(y, na.rm = TRUE); b0 <- 1
      lin_fit <- try(stats::lm(y ~ t), silent = TRUE)
      if (!inherits(lin_fit, "try-error")) {
        cf <- unname(stats::coef(lin_fit))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- cf[1]; b0 <- cf[2] }
      }
      list(a = a0, b = b0)
    },
    power = {
      # L(t) = a*t^b => log(y) ~ log(a) + b*log(t): a linear regression on
      # the log-log scale (excluding t <= 0, where the model is undefined).
      pos <- t > 0 & y > 0
      a0 <- 1; b0 <- 1
      lin_pow <- try(stats::lm(log(y[pos]) ~ log(t[pos])), silent = TRUE)
      if (!inherits(lin_pow, "try-error")) {
        cf <- unname(stats::coef(lin_pow))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- exp(cf[1]); b0 <- cf[2] }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- 1
      list(a = a0, b = b0)
    },
    exponential = {
      # L(t) = a*exp(b*t) => log(y) ~ log(a) + b*t: a linear regression on
      # the log scale (excluding non-positive sizes).
      pos <- y > 0
      a0 <- max(min(y, na.rm = TRUE), 1e-3); b0 <- 0.1
      lin_exp <- try(stats::lm(log(y[pos]) ~ t[pos]), silent = TRUE)
      if (!inherits(lin_exp, "try-error")) {
        cf <- unname(stats::coef(lin_exp))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- exp(cf[1]); b0 <- cf[2] }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- 1
      list(a = a0, b = b0)
    },
    logarithmic = {
      # L(t) = a + b*log(t): a linear regression of y on log(t) (excluding
      # t <= 0, where the model is undefined).
      pos <- t > 0
      a0 <- mean(y, na.rm = TRUE); b0 <- 1
      lin_log <- try(stats::lm(y[pos] ~ log(t[pos])), silent = TRUE)
      if (!inherits(lin_log, "try-error")) {
        cf <- unname(stats::coef(lin_log))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- cf[1]; b0 <- cf[2] }
      }
      list(a = a0, b = b0)
    },
    hyperbolic = {
      # L(t) = a*t/(b+t) -> 1/y = 1/a + (b/a)*(1/t) (Lineweaver-Burk-style
      # linearisation), fit on the reciprocal scale (excluding t <= 0, y <= 0).
      pos <- t > 0 & y > 0
      a0 <- Linf0; b0 <- max(stats::median(t, na.rm = TRUE), 1e-3)
      lin_hyp <- try(stats::lm(I(1 / y[pos]) ~ I(1 / t[pos])), silent = TRUE)
      if (!inherits(lin_hyp, "try-error")) {
        cf <- unname(stats::coef(lin_hyp))
        if (length(cf) == 2 && all(is.finite(cf)) && cf[1] > 0) {
          a0 <- 1 / cf[1]; b0 <- cf[2] / cf[1]
        }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- Linf0
      if (!is.finite(b0) || b0 <= 0) b0 <- max(stats::median(t, na.rm = TRUE), 1e-3)
      list(a = a0, b = b0)
    },
    ricker = {
      # L(t) = a*t*exp(-b*t) -> log(y/t) ~ log(a) - b*t (excluding t <= 0,
      # y <= 0); b is started positive (the characteristic rise-then-fall).
      pos <- t > 0 & y > 0
      a0 <- Linf0 / max(t, na.rm = TRUE); b0 <- 1 / max(t, na.rm = TRUE)
      lin_rk <- try(stats::lm(log(y[pos] / t[pos]) ~ t[pos]), silent = TRUE)
      if (!inherits(lin_rk, "try-error")) {
        cf <- unname(stats::coef(lin_rk))
        if (length(cf) == 2 && all(is.finite(cf))) {
          a0 <- exp(cf[1]); b0 <- -cf[2]
        }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- Linf0 / max(t, na.rm = TRUE)
      if (!is.finite(b0) || b0 <= 0) b0 <- 1 / max(t, na.rm = TRUE)
      list(a = a0, b = b0)
    },
    beverton_holt = {
      # Same Lineweaver-Burk-style linearisation as the hyperbolic model
      # above, adapted to L(t) = a*t/(1+b*t) -> 1/y = 1/(a*t) + b/a.
      pos <- t > 0 & y > 0
      a0 <- Linf0 / max(stats::median(t, na.rm = TRUE), 1e-3); b0 <- 1 / max(t, na.rm = TRUE)
      lin_bh <- try(stats::lm(I(1 / y[pos]) ~ I(1 / t[pos])), silent = TRUE)
      if (!inherits(lin_bh, "try-error")) {
        cf <- unname(stats::coef(lin_bh))
        if (length(cf) == 2 && all(is.finite(cf)) && cf[2] > 0) {
          a0 <- 1 / cf[2]; b0 <- cf[1] / cf[2]
        }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- Linf0 / max(stats::median(t, na.rm = TRUE), 1e-3)
      if (!is.finite(b0) || b0 <= 0) b0 <- 1 / max(t, na.rm = TRUE)
      list(a = a0, b = b0)
    },
    gamma = {
      # Starts from a plain power fit (g as the power exponent) with a
      # small decay rate b, since a full Gamma-shape linearisation is not
      # available in closed form.
      pos <- t > 0 & y > 0
      a0 <- 1; g0 <- 1
      lin_pow <- try(stats::lm(log(y[pos]) ~ log(t[pos])), silent = TRUE)
      if (!inherits(lin_pow, "try-error")) {
        cf <- unname(stats::coef(lin_pow))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- exp(cf[1]); g0 <- max(cf[2], 0.1) }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- 1
      b0 <- 0.05
      list(a = a0, b = b0, g = g0)
    },
    weibull = list(alpha = Linf0, kappa = K0, gamma = min(t0_0, 0), delta = 1),
    extended_power = {
      # Starts from the plain power fit (b as the large-t exponent), with
      # c = 0 (i.e. starts exactly at the plain power model).
      pos <- t > 0 & y > 0
      a0 <- 1; b0 <- 1
      lin_pow <- try(stats::lm(log(y[pos]) ~ log(t[pos])), silent = TRUE)
      if (!inherits(lin_pow, "try-error")) {
        cf <- unname(stats::coef(lin_pow))
        if (length(cf) == 2 && all(is.finite(cf))) { a0 <- exp(cf[1]); b0 <- cf[2] }
      }
      if (!is.finite(a0) || a0 <= 0) a0 <- 1
      list(a = a0, b = b0, c = 0)
    },
    johnson = {
      # t0 must stay below every observed age (k*(t-t0) != 0 for any data
      # point); if the VB-style linearisation above didn't already land
      # there, fall back to half an age-range below the minimum age.
      age_span <- max(t, na.rm = TRUE) - min(t, na.rm = TRUE)
      t0_guess <- if (is.finite(t0_0) && t0_0 < min(t, na.rm = TRUE)) {
        t0_0
      } else {
        min(t, na.rm = TRUE) - 0.5 * max(age_span, 1)
      }
      list(Linf = Linf0, k = K0, t0 = t0_guess)
    },
    .custom_initial_values(model, t, y)
  )
}

#' Starting values for a model added via add_growth_model()
#'
#' Fallback used by \code{\link{initial_values}} for any model not in its
#' built-in \code{switch()} - i.e. one registered at runtime with
#' \code{\link{add_growth_model}}. If a \code{start} heuristic was
#' registered for it (see \code{add_growth_model}'s \code{start}
#' argument), that is used; otherwise every parameter defaults to 1, with
#' a message pointing at the two ways to do better (a \code{start}
#' heuristic registered once at \code{add_growth_model()} time, or a
#' one-off \code{start} passed directly to \code{fit_growth()}).
#'
#' @param model Model name.
#' @param t Vector of ages.
#' @param y Vector of sizes.
#' @return Named list of starting values.
#' @keywords internal
.custom_initial_values <- function(model, t, y) {
  if (!(model %in% names(GROWTH_MODELS))) {
    stop("Unrecognised model: ", model)
  }
  spec <- GROWTH_MODELS[[model]]
  if (!is.null(spec$start)) {
    start <- if (is.function(spec$start)) spec$start(t = t, y = y) else spec$start
    missing_p <- setdiff(spec$parameters, names(start))
    if (length(missing_p) > 0) {
      stop("initial_values(): the 'start' heuristic registered for model '", model,
           "' did not return a value for parameter(s): ", paste(missing_p, collapse = ", "),
           ". Fix the heuristic passed to add_growth_model(), or pass start = list(...) ",
           "directly to fit_growth().")
    }
    return(start[spec$parameters])
  }
  .auto_multistart(fn = spec$fn, parameters = spec$parameters, t = t, y = y, model = model)
}

#' Automatic multi-start search for starting values
#'
#' Used by \code{\link{.custom_initial_values}} whenever a model added via
#' \code{\link{add_growth_model}} has no registered \code{start} heuristic.
#' A user-defined \code{fn} is a black box - the package has no built-in
#' knowledge of what a reasonable \code{Linf} or \code{K} would be for it,
#' the way it does for the built-in models' hand-written linearisations.
#' Rather than falling back to a single naive guess (every parameter = 1,
#' which converges poorly whenever the true parameter scales differ from
#' 1), this runs a bounded search over many candidate parameter
#' combinations - built from a magnitude ladder (dimensionless orders of
#' magnitude, for rate/shape-like parameters), the same ladder scaled to
#' the span of \code{y} (for size-like parameters) and to the span of
#' \code{t} (for time-scale-like parameters), plus a few age-anchored
#' candidates at \code{min(t)}/\code{max(t)} (for location-like
#' parameters, e.g. a \code{t0}-style reference age) - and scores each
#' candidate directly against the data (residual sum of squares from one
#' call to \code{fn}, no optimiser run), returning the best-scoring
#' combination. This means \code{start} in \code{\link{add_growth_model}}
#' is genuinely optional: a user does not need to already understand
#' their own model's mathematics (the very thing they may be turning to
#' \code{add_growth_model} to avoid) to get a working fit. A registered
#' \code{start} heuristic, or one passed directly to \code{fit_growth()},
#' is still useful for speed or for unusually difficult/multi-modal
#' models, and always takes priority over this search when supplied.
#'
#' The random component of the search (used once there are 3 or more
#' parameters, where a full factorial grid would be too large) is made
#' reproducible by temporarily seeding and then restoring R's random
#' number generator state, so it does not disturb the caller's RNG
#' stream (e.g. inside \code{\link{simulate_growth_data}} or
#' \code{\link{fit_growth_bayes}}).
#'
#' @param fn The model function, as stored in \code{GROWTH_MODELS[[model]]$fn}.
#' @param parameters Character vector of \code{fn}'s parameter names.
#' @param t Vector of ages.
#' @param y Vector of sizes.
#' @param model Model name (used only for the diagnostic message on total failure).
#' @param max_evals Upper bound on the number of candidate combinations scored.
#' @return Named list of starting values.
#' @keywords internal
.auto_multistart <- function(fn, parameters, t, y, model, max_evals = 4000) {
  p <- length(parameters)
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]

  y_scale <- max(abs(diff(range(y))), abs(mean(y)), 1)
  t_min <- min(t); t_max <- max(t)
  t_span <- max(t_max - t_min, 1)

  ladder <- c(-100, -10, -3, -1, -0.3, -0.1, -0.01, 0.01, 0.1, 0.3, 1, 3, 10, 100)
  candidates <- unique(c(
    ladder,
    ladder * y_scale,
    ladder * t_span,
    t_min, t_max, (t_min + t_max) / 2, t_min - 0.5 * t_span
  ))
  candidates <- candidates[is.finite(candidates)]

  score <- function(theta) {
    pred <- tryCatch(
      do.call(fn, c(list(t = t), stats::setNames(as.list(theta), parameters))),
      error = function(e) NA_real_
    )
    if (!is.numeric(pred) || length(pred) != length(t)) return(Inf)
    finite <- is.finite(pred)
    if (mean(finite) < 0.8) return(Inf)
    sum((y[finite] - pred[finite])^2)
  }

  best_theta <- stats::setNames(rep(1, p), parameters)
  best_rss <- score(best_theta)

  # Seed only the search's own randomness, then restore the caller's RNG
  # state exactly as it was - this function must not be observable via
  # side effects on later random draws (simulation, bootstrap, MCMC, ...).
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(20240610)

  if (p <= 2) {
    grid <- as.matrix(expand.grid(rep(list(candidates), p)))
    if (nrow(grid) > max_evals) {
      grid <- grid[sample.int(nrow(grid), max_evals), , drop = FALSE]
    }
  } else {
    n_random <- min(max_evals, 6000)
    grid <- matrix(sample(candidates, n_random * p, replace = TRUE), nrow = n_random, ncol = p)
    n_extra <- min(500, n_random)
    extra <- matrix(
      sample(c(-1, 1), n_extra * p, replace = TRUE) * 10^stats::runif(n_extra * p, -2, 3),
      nrow = n_extra, ncol = p
    )
    grid <- rbind(grid, extra)
  }

  for (i in seq_len(nrow(grid))) {
    theta <- grid[i, ]
    rss <- score(theta)
    if (is.finite(rss) && rss < best_rss) {
      best_rss <- rss
      best_theta <- stats::setNames(theta, parameters)
    }
  }

  if (!is.finite(best_rss)) {
    message("initial_values(): the automatic multi-start search for custom model '", model,
             "' could not find any starting point where fn() returns mostly finite values ",
             "over the observed data - defaulting every parameter to 1, which will likely ",
             "fail to converge. Consider passing start = list(...) directly to fit_growth(), ",
             "or registering a start heuristic via the 'start' argument of add_growth_model().")
    return(stats::setNames(as.list(rep(1, p)), parameters))
  }

  as.list(best_theta)
}

# ---- Shared reference-age fixing (Schnute cases, Francis) --------------------

#' @keywords internal
.build_fixed_fn <- function(spec, t) {
  if (isTRUE(spec$needs_t1_t2)) {
    t1_fixed <- min(t); t2_fixed <- max(t)
    base_fn <- spec$fn
    fn_local <- function(t, ...) {
      do.call(base_fn, c(list(t = t), list(...), list(t1 = t1_fixed, t2 = t2_fixed)))
    }
    return(list(fn_local = fn_local, fixed_extra = list(t1 = t1_fixed, t2 = t2_fixed)))
  }
  if (isTRUE(spec$needs_t1_t3)) {
    t1_fixed <- min(t); t3_fixed <- max(t); t2_fixed <- (t1_fixed + t3_fixed) / 2
    base_fn <- spec$fn
    fn_local <- function(t, ...) {
      do.call(base_fn, c(list(t = t), list(...),
                          list(t1 = t1_fixed, t2 = t2_fixed, t3 = t3_fixed)))
    }
    return(list(fn_local = fn_local,
                fixed_extra = list(t1 = t1_fixed, t2 = t2_fixed, t3 = t3_fixed)))
  }
  list(fn_local = spec$fn, fixed_extra = list())
}

# ---- Shared fitting core (LS and MLE) ----------------------------------------

#' @keywords internal
.fit_core <- function(y, predictors, fn, parameters, method, distribution,
                       start, control_ls = NULL, control_mle = NULL) {
  method <- match.arg(method, c("ls", "mle"))
  distribution <- match.arg(distribution, c("normal", "lognormal"))
  start <- start[parameters]

  if (method == "ls") {
    .model_fn <- fn
    pred_names <- names(predictors)
    call_args <- paste(sprintf("%s = %s", c(pred_names, parameters),
                                c(pred_names, parameters)), collapse = ", ")
    form_str <- paste0("y ~ .model_fn(", call_args, ")")
    model_formula <- stats::as.formula(form_str, env = environment())
    nls_data <- c(predictors, list(y = y))
    ctrl <- if (is.null(control_ls)) minpack.lm::nls.lm.control(maxiter = 1000) else control_ls
    fit <- minpack.lm::nlsLM(model_formula, start = start, data = nls_data,
                              control = ctrl)
    return(list(fit_obj = fit, method = "ls", distribution = "normal"))
  }

  # method == "mle"
  nll_vec <- function(p) {
    params <- as.list(p[parameters])
    yhat <- do.call(fn, c(predictors, params))
    log_sigma <- p[["log_sigma"]]
    sigma <- exp(log_sigma)
    resid <- if (distribution == "normal") y - yhat else log(y) - log(pmax(yhat, 1e-8))
    if (any(!is.finite(yhat)) || any(!is.finite(resid)) || !is.finite(sigma) || sigma <= 0) {
      return(1e10)
    }
    -sum(stats::dnorm(resid, mean = 0, sd = sigma, log = TRUE))
  }

  yhat0 <- do.call(fn, c(predictors, start))
  resid0 <- if (distribution == "normal") y - yhat0 else log(y) - log(pmax(yhat0, 1e-8))
  sigma0 <- stats::sd(resid0, na.rm = TRUE)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  start_vec <- c(unlist(start), log_sigma = log(sigma0))

  ctrl_mle <- if (is.null(control_mle)) list(maxit = 5000, reltol = 1e-10) else control_mle
  bbmle::parnames(nll_vec) <- names(start_vec)
  fit <- bbmle::mle2(minuslogl = nll_vec, start = start_vec, vecpar = TRUE,
                      method = "Nelder-Mead", control = ctrl_mle)
  list(fit_obj = fit, method = "mle", distribution = distribution)
}

# ---- Main interface: fit_growth ----------------------------------------------

#' Fit an age-size growth model
#'
#' Fits any of the models in \code{\link{GROWTH_MODELS}} (including
#' Schnute) to age-size data, either by nonlinear least squares
#' (\code{method = "ls"}, using \code{minpack.lm::nlsLM}, Levenberg-Marquardt
#' algorithm) or by maximum likelihood (\code{method = "mle"}, using
#' \code{bbmle::mle2}), assuming additive normal or multiplicative
#' lognormal errors.
#'
#' @param data Data frame with the data.
#' @param age Name (character) of the age column in \code{data}.
#' @param size Name (character) of the size/weight column in \code{data}.
#' @param model Name of the model to fit (see \code{names(GROWTH_MODELS)}).
#' @param method \code{"ls"} (least squares) or \code{"mle"} (maximum likelihood).
#' @param distribution For \code{method = "mle"}: \code{"normal"} (additive
#'   error) or \code{"lognormal"} (multiplicative error). Ignored if
#'   \code{method = "ls"}.
#' @param start Named list of starting values. If \code{NULL} (default), it
#'   is computed automatically with \code{\link{initial_values}}.
#' @param control_ls Control list for \code{minpack.lm::nls.lm.control}.
#' @param control_mle Control list for \code{optim} (via \code{bbmle::mle2}).
#' @return An object of class \code{"growth_fit"}.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' summary(fit)
#'
#' # Maximum likelihood, lognormal errors:
#' fit_mle <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy",
#'                        method = "mle", distribution = "lognormal")
#' fit_mle
#' @export
fit_growth <- function(data, age, size, model = "von_bertalanffy",
                        method = c("ls", "mle"),
                        distribution = c("normal", "lognormal"),
                        start = NULL, control_ls = NULL, control_mle = NULL) {
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  stopifnot(is.data.frame(data), age %in% names(data), size %in% names(data))
  stopifnot(model %in% names(GROWTH_MODELS))

  t <- data[[age]]
  y <- data[[size]]
  ok <- is.finite(t) & is.finite(y)
  if (any(!ok)) warning(sum(!ok), " observations with NA/Inf were excluded.")
  t <- t[ok]; y <- y[ok]

  spec <- GROWTH_MODELS[[model]]
  if (is.null(start)) start <- initial_values(t, y, model)

  built <- .build_fixed_fn(spec, t)
  fn_local <- built$fn_local
  fixed_extra <- built$fixed_extra

  core <- .fit_core(y = y, predictors = list(t = t), fn = fn_local,
                     parameters = spec$parameters, method = method,
                     distribution = distribution, start = start,
                     control_ls = control_ls, control_mle = control_mle)

  .build_growth_fit(core = core, model = model, spec = spec, t = t, y = y,
                     age = age, size = size, raw_data = data,
                     fixed_extra = fixed_extra, fn_local = fn_local)
}

#' @keywords internal
.build_growth_fit <- function(core, model, spec, t, y, age, size, raw_data,
                               fixed_extra, fn_local) {
  fit_obj <- core$fit_obj
  n <- length(y)
  if (core$method == "ls") {
    ll <- as.numeric(stats::logLik(fit_obj))
    coefs <- stats::coef(fit_obj)
    sigma <- summary(fit_obj)$sigma
  } else {
    ll <- -fit_obj@min
    coefs <- bbmle::coef(fit_obj)[spec$parameters]
    sigma <- exp(bbmle::coef(fit_obj)[["log_sigma"]])
  }
  # Number of parameters reported as "k": the growth model's own structural
  # parameters only (e.g. 3 for von Bertalanffy: Linf, K, t0), NOT counting
  # the residual error variance/sigma as an extra parameter. This matches
  # the convention most commonly used in the fish-growth model-selection
  # literature (e.g. Cerrato 1990; Katsanevakis 2006), and keeps "k" equal
  # to what a reader intuitively expects when they know a model's formula.
  # `stats::logLik.nls()`/`bbmle::coef()` would instead count sigma as an
  # additional parameter (k + 1); since every model here estimates exactly
  # one such extra parameter, that convention only shifts AIC/BIC by the
  # same constant for every model and never changes AIC/BIC rankings - it
  # only makes a (typically negligible) difference to AICc's small-sample
  # correction term, which depends on k non-linearly.
  k <- length(spec$parameters)
  aic <- -2 * ll + 2 * k
  aicc <- if (n - k - 1 > 0) aic + (2 * k * (k + 1)) / (n - k - 1) else NA_real_
  bic <- -2 * ll + k * log(n)

  out <- list(
    model = model,
    label = spec$label,
    method = core$method,
    distribution = core$distribution,
    fit_obj = fit_obj,
    coefficients = coefs,
    sigma = sigma,
    logLik = ll,
    k = k,
    n = n,
    aic = aic,
    aicc = aicc,
    bic = bic,
    t = t,
    y = y,
    age = age,
    size = size,
    data = raw_data,
    spec = spec,
    fixed_extra = fixed_extra,
    fn = fn_local
  )
  class(out) <- "growth_fit"
  out
}

# ---- S3 methods ----------------------------------------------------------------

#' Print a growth_fit object
#'
#' Prints a short summary of a fitted growth model: label/formula, fitting
#' method and error distribution, sample size, number of estimated
#' parameters, log-likelihood, information criteria (AIC/AICc/BIC), and
#' the point estimates of the coefficients.
#'
#' @param x Object of class \code{"growth_fit"} (from \code{\link{fit_growth}}).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' print(fit)   # or just `fit` at the console
#' @export
print.growth_fit <- function(x, ...) {
  cat(sprintf("Growth model: %s (%s)\n", x$label, x$model))
  cat(sprintf("Formula: %s\n", x$spec$formula))
  cat(sprintf("Method: %s | Error distribution: %s\n",
              toupper(x$method), x$distribution))
  cat(sprintf("n = %d | estimated parameters = %d | logLik = %.3f\n",
              x$n, x$k, x$logLik))
  cat(sprintf("AIC = %.3f | AICc = %.3f | BIC = %.3f\n", x$aic, x$aicc, x$bic))
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

#' Summarise a growth_fit object
#'
#' Prints everything \code{\link{print.growth_fit}} prints, plus the
#' residual standard deviation (\code{sigma}).
#'
#' @param object Object of class \code{"growth_fit"}.
#' @param ... Unused.
#' @return \code{object}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' summary(fit)
#' @export
summary.growth_fit <- function(object, ...) {
  print(object)
  cat(sprintf("Sigma (residual standard deviation): %.4f\n", object$sigma))
  invisible(object)
}

#' Extract coefficients from a growth_fit object
#'
#' @param object Object of class \code{"growth_fit"}.
#' @param ... Unused.
#' @return Named numeric vector of estimated parameters.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' coef(fit)
#' @export
coef.growth_fit <- function(object, ...) object$coefficients

#' Extract the log-likelihood from a growth_fit object
#'
#' @param object Object of class \code{"growth_fit"}.
#' @param ... Unused.
#' @return An object of class \code{"logLik"}, with the number of
#'   estimated parameters attached as its \code{"df"} attribute (see
#'   \code{\link[stats]{logLik}}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' logLik(fit)
#' @export
logLik.growth_fit <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$k
  class(val) <- "logLik"
  val
}

#' Extract (or recompute) the AIC of a growth_fit object
#'
#' @param object Object of class \code{"growth_fit"}.
#' @param ... Unused.
#' @param k Penalty per estimated parameter (default 2, the classic AIC).
#'   For \code{k = 2} the already-stored \code{object$aic} is returned;
#'   any other value recomputes \code{-2*logLik + k*(number of parameters)}.
#' @return Numeric AIC value.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' AIC(fit)
#' AIC(fit, k = log(nrow(growth_data)))   # BIC penalty, computed manually
#' @export
AIC.growth_fit <- function(object, ..., k = 2) {
  if (k == 2) object$aic else -2 * object$logLik + k * object$k
}

#' Predict sizes from a growth_fit object
#'
#' @param object Object of class \code{"growth_fit"}.
#' @param new_t Vector of ages (or times) at which to predict. If \code{NULL},
#'   a regular grid over the observed range is used.
#' @param n Number of grid points if \code{new_t} is \code{NULL}.
#' @param ... Unused.
#' @return Data frame with columns \code{t} and \code{y_pred}.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' predict(fit, new_t = c(1, 5, 10))
#' head(predict(fit))
#' @export
predict.growth_fit <- function(object, new_t = NULL, n = 200, ...) {
  if (is.null(new_t)) {
    new_t <- seq(min(object$t), max(object$t), length.out = n)
  }
  params <- as.list(object$coefficients)
  y_pred <- do.call(object$fn, c(list(t = new_t), params))
  data.frame(t = new_t, y_pred = as.numeric(y_pred))
}

#' Fit several growth models to the same data (multi-model)
#'
#' Fits, in a single call, a collection of growth models to the same
#' age-size data. Models that fail to converge are dropped with a
#' warning, rather than stopping the whole run.
#'
#' @param data Data frame with the data.
#' @param age Name of the age column.
#' @param size Name of the size column.
#' @param models Vector of model names to fit (by default, all available
#'   age-size models).
#' @param method \code{"ls"} or \code{"mle"}.
#' @param distribution \code{"normal"} or \code{"lognormal"} (only if
#'   \code{method = "mle"}).
#' @param ... Additional arguments passed to \code{\link{fit_growth}}.
#' @return Named list of \code{"growth_fit"} objects (class
#'   \code{"multimodel_growth"}), in the same order as \code{models}.
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' aic_table(mm)
#' @export
fit_multimodel <- function(data, age, size,
                            models = names(GROWTH_MODELS),
                            method = c("ls", "mle"),
                            distribution = c("normal", "lognormal"), ...) {
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  results <- list()
  for (m in models) {
    fit <- try(
      fit_growth(data, age, size, model = m, method = method,
                 distribution = distribution, ...),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      warning(sprintf("Model '%s' did not converge and was skipped: %s", m,
                       conditionMessage(attr(fit, "condition"))))
    } else {
      results[[m]] <- fit
    }
  }
  class(results) <- "multimodel_growth"
  results
}
