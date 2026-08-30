#' @importFrom stats rnorm dnorm cov quantile median var acf sd
NULL

# ---- Priors --------------------------------------------------------------

#' @keywords internal
.default_priors <- function(start, log_sigma0) {
  priors <- lapply(start, function(v) {
    centre <- v
    scale <- max(2 * abs(v), 1)
    function(x) stats::dnorm(x, mean = centre, sd = scale, log = TRUE)
  })
  priors[["log_sigma"]] <- (function(centre) {
    function(x) stats::dnorm(x, mean = centre, sd = 2, log = TRUE)
  })(log_sigma0)
  priors
}

# ---- Likelihood (pointwise, for WAIC) and log-posterior ------------------

#' @keywords internal
.loglik_pointwise_bayes <- function(p, y, predictors, fn, parameters, distribution) {
  log_sigma <- p[["log_sigma"]]
  sigma <- exp(log_sigma)
  if (!is.finite(sigma) || sigma <= 0) return(rep(-Inf, length(y)))
  params <- as.list(p[parameters])
  yhat <- do.call(fn, c(predictors, params))
  resid <- if (distribution == "normal") y - yhat else log(y) - log(pmax(yhat, 1e-8))
  if (any(!is.finite(yhat)) || any(!is.finite(resid))) return(rep(-Inf, length(y)))
  stats::dnorm(resid, mean = 0, sd = sigma, log = TRUE)
}

#' @keywords internal
.log_posterior_bayes <- function(p, y, predictors, fn, parameters, distribution, priors) {
  ll <- .loglik_pointwise_bayes(p, y, predictors, fn, parameters, distribution)
  if (any(!is.finite(ll))) return(-Inf)
  lp <- 0
  for (pn in names(priors)) {
    lp <- lp + priors[[pn]](p[[pn]])
    if (!is.finite(lp)) return(-Inf)
  }
  sum(ll) + lp
}

# ---- Adaptive Metropolis-Hastings sampler ---------------------------------

#' @keywords internal
.adaptive_metropolis <- function(log_post_fn, start_vec, n_iter, n_burnin,
                                  adapt_every = 100, target_accept = 0.234,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- length(start_vec)
  pnames <- names(start_vec)
  chain <- matrix(NA_real_, nrow = n_iter, ncol = d, dimnames = list(NULL, pnames))

  theta <- start_vec
  lp_curr <- log_post_fn(theta)
  tries <- 0
  while (!is.finite(lp_curr) && tries < 200) {
    theta <- start_vec + stats::rnorm(d, sd = pmax(abs(start_vec), 1) * 0.05)
    names(theta) <- pnames
    lp_curr <- log_post_fn(theta)
    tries <- tries + 1
  }
  if (!is.finite(lp_curr)) {
    stop("Could not find a finite starting log-posterior after ", tries,
         " attempts; check 'start'/'priors' for this model.")
  }

  base_scale <- pmax(abs(start_vec) * 0.1, 0.05)
  cov_mat <- diag(base_scale^2, d)
  scale_factor <- 2.38^2 / d
  n_accept_recent <- 0

  for (i in seq_len(n_iter)) {
    z <- stats::rnorm(d)
    step <- as.numeric(t(chol(scale_factor * cov_mat + diag(1e-10, d))) %*% z)
    prop <- theta + step
    names(prop) <- pnames
    lp_prop <- log_post_fn(prop)

    if (is.finite(lp_prop) && log(stats::runif(1)) < (lp_prop - lp_curr)) {
      theta <- prop
      lp_curr <- lp_prop
      n_accept_recent <- n_accept_recent + 1
    }
    chain[i, ] <- theta

    # Diminishing adaptation: only during burn-in, and only using the
    # burn-in history so far, so the sampler settles into a fixed proposal
    # once sampling ("post-burn-in") begins - required for a valid
    # (asymptotically correct) adaptive MCMC chain.
    if (i <= n_burnin && i %% adapt_every == 0) {
      hist_window <- chain[max(1, i - adapt_every + 1):i, , drop = FALSE]
      if (nrow(hist_window) > d) {
        emp_cov <- stats::cov(hist_window)
        if (all(is.finite(emp_cov))) cov_mat <- emp_cov + diag(1e-8, d)
      }
      accept_rate <- n_accept_recent / adapt_every
      scale_factor <- scale_factor * exp((accept_rate - target_accept) * 0.5)
      scale_factor <- max(min(scale_factor, 1e4), 1e-6)
      n_accept_recent <- 0
    }
  }

  list(chain = chain, n_accept_total = sum(diff(chain[, 1]) != 0, na.rm = TRUE))
}

# ---- Diagnostics -----------------------------------------------------------

#' @keywords internal
.rhat_gelman_rubin <- function(chains_list) {
  # chains_list: list of post-burn-in matrices (n_samples x d), same d/names.
  m <- length(chains_list)
  if (m < 2) return(stats::setNames(rep(NA_real_, ncol(chains_list[[1]])), colnames(chains_list[[1]])))
  n <- nrow(chains_list[[1]])
  pnames <- colnames(chains_list[[1]])
  vapply(pnames, function(pn) {
    vals <- lapply(chains_list, function(ch) ch[, pn])
    chain_means <- vapply(vals, mean, numeric(1))
    chain_vars <- vapply(vals, stats::var, numeric(1))
    W <- mean(chain_vars)
    B <- n * stats::var(chain_means)
    var_hat <- ((n - 1) / n) * W + B / n
    if (!is.finite(W) || W <= 0) return(NA_real_)
    sqrt(var_hat / W)
  }, numeric(1))
}

#' @keywords internal
.effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 10 || stats::sd(x) <= 0) return(NA_real_)
  ac <- try(stats::acf(x, plot = FALSE, lag.max = min(n - 1, 1000)), silent = TRUE)
  if (inherits(ac, "try-error")) return(NA_real_)
  rho <- as.numeric(ac$acf)[-1]
  pos <- which(rho <= 0)
  cutoff <- if (length(pos) > 0) pos[1] - 1 else length(rho)
  cutoff <- max(cutoff, 0)
  denom <- 1 + 2 * sum(rho[seq_len(cutoff)])
  ess <- n / max(denom, 1e-6)
  min(max(ess, 1), n)
}

# ---- WAIC and DIC -----------------------------------------------------------

#' @keywords internal
.compute_waic <- function(log_lik_matrix) {
  # log_lik_matrix: S (samples) x N (observations)
  S <- nrow(log_lik_matrix); N <- ncol(log_lik_matrix)
  lppd_i <- apply(log_lik_matrix, 2, function(col) {
    m <- max(col)
    m + log(mean(exp(col - m)))
  })
  p_waic_i <- apply(log_lik_matrix, 2, stats::var)
  elpd_i <- lppd_i - p_waic_i
  lppd <- sum(lppd_i)
  p_waic <- sum(p_waic_i)
  waic <- -2 * (lppd - p_waic)
  se_waic <- sqrt(N * stats::var(elpd_i)) * 2
  list(waic = waic, elpd_waic = lppd - p_waic, p_waic = p_waic,
       se_waic = se_waic, pointwise = elpd_i)
}

#' @keywords internal
.compute_dic <- function(log_lik_matrix, log_lik_at_mean) {
  deviance_samples <- -2 * rowSums(log_lik_matrix)
  D_bar <- mean(deviance_samples)
  D_hat <- -2 * sum(log_lik_at_mean)
  p_D <- D_bar - D_hat
  list(dic = D_bar + p_D, p_D = p_D, D_bar = D_bar)
}

# ---- Main interface: fit_growth_bayes --------------------------------------

#' Fit an age-size growth model by Bayesian inference (MCMC)
#'
#' Fits any of the models in \code{\link{GROWTH_MODELS}} using a
#' self-contained adaptive Metropolis-Hastings sampler (Haario et al.
#' 2001) - no external MCMC engine (Stan/JAGS/NIMBLE) is required, keeping
#' the package installable as a single source tarball, consistent with
#' the rest of MMgrowthR. Several independent chains are run so that
#' convergence can be checked (Gelman-Rubin \eqn{\hat R}), and posterior
#' credible intervals, WAIC, and DIC are computed automatically.
#'
#' \strong{Priors}: by default, each model parameter gets a weakly
#' informative Normal prior centred on the automatic starting value
#' (\code{\link{initial_values}}) with a generous standard deviation
#' (\code{max(2*|start|, 1)}); the residual log-standard-deviation
#' (\code{log_sigma}) gets \code{Normal(log(sigma0), 2)}. Override any of
#' these via \code{priors}, a named list of functions taking a parameter
#' value and returning its log-density (e.g.
#' \code{priors = list(Linf = function(x) dnorm(x, 80, 15, log = TRUE))});
#' parameters not named in \code{priors} keep their default.
#'
#' \strong{Likelihood}: identical construction (additive normal or
#' log-scale normal for \code{distribution = "lognormal"}) to
#' \code{\link{fit_growth}(method = "mle")}, so that WAIC/DIC from this
#' function are directly comparable in spirit to the AIC/AICc/BIC
#' produced by \code{\link{aic_table}} for the same data.
#'
#' @param data Data frame with the data.
#' @param age Name (character) of the age column in \code{data}.
#' @param size Name (character) of the size/weight column in \code{data}.
#' @param model Name of the model to fit (see \code{names(GROWTH_MODELS)}).
#' @param distribution \code{"normal"} (additive error) or
#'   \code{"lognormal"} (multiplicative error).
#' @param priors Optional named list overriding the default priors (see Details).
#' @param start Named list of starting values / prior centres. If
#'   \code{NULL} (default), computed automatically with \code{\link{initial_values}}.
#' @param n_iter Number of MCMC iterations per chain (including burn-in).
#' @param n_burnin Number of initial iterations discarded as burn-in
#'   (also the adaptation period - the proposal is frozen afterwards).
#' @param n_chains Number of independent chains (\code{>= 2} recommended,
#'   needed to compute \eqn{\hat R}).
#' @param thin Keep only every \code{thin}-th post-burn-in sample.
#' @param level Credible interval level (default 0.95, equal-tailed/percentile).
#' @param seed Optional random seed (chains use \code{seed + chain_index - 1}).
#' @param verbose If \code{TRUE} (default), reports acceptance rates and \eqn{\hat R}.
#' @return An object of class \code{"growth_fit_bayes"} with the posterior
#'   samples (\code{samples}), a per-parameter summary table
#'   (\code{summary}: mean, median, sd, credible interval), \code{waic},
#'   \code{dic}, \code{rhat}, \code{ess}, and the elements needed by
#'   \code{\link{predict.growth_fit_bayes}}/\code{\link{plot_growth_fit_bayes}}.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' summary(fit_b)
#' }
#' @export
fit_growth_bayes <- function(data, age, size, model = "von_bertalanffy",
                              distribution = c("normal", "lognormal"),
                              priors = NULL, start = NULL,
                              n_iter = 8000, n_burnin = 2000, n_chains = 3,
                              thin = 1, level = 0.95, seed = NULL, verbose = TRUE) {
  distribution <- match.arg(distribution)
  stopifnot(is.data.frame(data), age %in% names(data), size %in% names(data))
  stopifnot(model %in% names(GROWTH_MODELS))
  stopifnot(n_chains >= 1, n_iter > n_burnin)

  t <- data[[age]]; y <- data[[size]]
  ok <- is.finite(t) & is.finite(y)
  if (any(!ok)) warning(sum(!ok), " observations with NA/Inf were excluded.")
  t <- t[ok]; y <- y[ok]

  spec <- GROWTH_MODELS[[model]]
  if (is.null(start)) start <- initial_values(t, y, model)
  built <- .build_fixed_fn(spec, t)
  fn_local <- built$fn_local
  fixed_extra <- built$fixed_extra
  parameters <- spec$parameters

  yhat0 <- do.call(fn_local, c(list(t = t), start))
  resid0 <- if (distribution == "normal") y - yhat0 else log(y) - log(pmax(yhat0, 1e-8))
  sigma0 <- stats::sd(resid0, na.rm = TRUE)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  start_vec <- c(unlist(start[parameters]), log_sigma = log(sigma0))

  default_priors <- .default_priors(start_vec, log(sigma0))
  if (!is.null(priors)) {
    for (pn in names(priors)) default_priors[[pn]] <- priors[[pn]]
  }

  log_post_fn <- function(p) {
    .log_posterior_bayes(p, y = y, predictors = list(t = t), fn = fn_local,
                          parameters = parameters, distribution = distribution,
                          priors = default_priors)
  }

  chains_raw <- vector("list", n_chains)
  accept_rates <- numeric(n_chains)
  for (ch in seq_len(n_chains)) {
    ch_seed <- if (!is.null(seed)) seed + ch - 1 else NULL
    jitter <- if (ch == 1) rep(0, length(start_vec)) else
      stats::rnorm(length(start_vec), sd = pmax(abs(start_vec), 1) * 0.1)
    start_ch <- start_vec + jitter
    names(start_ch) <- names(start_vec)
    res <- .adaptive_metropolis(log_post_fn, start_ch, n_iter = n_iter,
                                 n_burnin = n_burnin, seed = ch_seed)
    chains_raw[[ch]] <- res$chain
    accept_rates[ch] <- res$n_accept_total / n_iter
    if (isTRUE(verbose)) {
      message(sprintf("Chain %d/%d: acceptance rate = %.2f", ch, n_chains, accept_rates[ch]))
    }
  }

  keep_idx <- seq(n_burnin + 1, n_iter, by = thin)
  chains_post <- lapply(chains_raw, function(ch) ch[keep_idx, , drop = FALSE])
  rhat <- .rhat_gelman_rubin(chains_post)
  samples <- do.call(rbind, chains_post)
  ess <- apply(samples, 2, .effective_sample_size)

  if (isTRUE(verbose) && n_chains >= 2) {
    message("Rhat (Gelman-Rubin, target ~1.00-1.05): ",
            paste(sprintf("%s = %.3f", names(rhat), rhat), collapse = ", "))
  }

  # Pointwise log-likelihood matrix (posterior samples x observations) for WAIC/DIC.
  log_lik_matrix <- t(apply(samples, 1, function(p) {
    .loglik_pointwise_bayes(p, y = y, predictors = list(t = t), fn = fn_local,
                             parameters = parameters, distribution = distribution)
  }))
  finite_rows <- apply(log_lik_matrix, 1, function(r) all(is.finite(r)))
  if (!all(finite_rows)) {
    log_lik_matrix <- log_lik_matrix[finite_rows, , drop = FALSE]
    samples <- samples[finite_rows, , drop = FALSE]
  }

  waic <- .compute_waic(log_lik_matrix)
  post_mean <- colMeans(samples)
  ll_at_mean <- .loglik_pointwise_bayes(post_mean, y = y, predictors = list(t = t),
                                         fn = fn_local, parameters = parameters,
                                         distribution = distribution)
  dic <- .compute_dic(log_lik_matrix, ll_at_mean)

  alpha <- 1 - level
  summary_tbl <- data.frame(
    parameter = colnames(samples),
    mean = colMeans(samples),
    median = apply(samples, 2, stats::median),
    sd = apply(samples, 2, stats::sd),
    ci_lower = apply(samples, 2, stats::quantile, probs = alpha / 2),
    ci_upper = apply(samples, 2, stats::quantile, probs = 1 - alpha / 2),
    rhat = rhat[colnames(samples)],
    ess = ess[colnames(samples)],
    stringsAsFactors = FALSE
  )
  rownames(summary_tbl) <- NULL

  coefficients <- post_mean[parameters]

  out <- list(
    model = model, label = spec$label, spec = spec,
    distribution = distribution,
    t = t, y = y, age = age, size = size, data = data,
    fixed_extra = fixed_extra, fn = fn_local,
    parameters = parameters,
    priors = default_priors,
    n_iter = n_iter, n_burnin = n_burnin, n_chains = n_chains, thin = thin,
    accept_rates = accept_rates,
    chains = chains_post,
    samples = samples,
    log_lik_matrix = log_lik_matrix,
    coefficients = coefficients,
    summary = summary_tbl,
    rhat = rhat, ess = ess,
    waic = waic, dic = dic,
    level = level,
    n = length(y), k = length(parameters)
  )
  class(out) <- "growth_fit_bayes"
  out
}

#' Fit several growth models to the same data by Bayesian inference
#'
#' Bayesian analogue of \code{\link{fit_multimodel}}: fits, in a single
#' call, a collection of growth models to the same data via
#' \code{\link{fit_growth_bayes}}. Models that fail (e.g. the sampler
#' cannot find a finite starting point) are dropped with a warning.
#'
#' @inheritParams fit_growth_bayes
#' @param models Vector of model names to fit (default: all available models).
#' @param ... Additional arguments passed to \code{\link{fit_growth_bayes}}.
#' @return Named list of \code{"growth_fit_bayes"} objects (class
#'   \code{"multimodel_growth_bayes"}).
#' @examples
#' data(growth_data)
#' \donttest{
#' mm_b <- fit_multimodel_bayes(growth_data, "age", "size",
#'                               models = c("von_bertalanffy", "gompertz"),
#'                               n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' bayes_table(mm_b)
#' }
#' @export
fit_multimodel_bayes <- function(data, age, size, models = names(GROWTH_MODELS),
                                  distribution = c("normal", "lognormal"), ...) {
  distribution <- match.arg(distribution)
  results <- list()
  for (m in models) {
    fit <- try(fit_growth_bayes(data, age, size, model = m,
                                 distribution = distribution, ...), silent = TRUE)
    if (inherits(fit, "try-error")) {
      warning(sprintf("Model '%s' could not be fitted by fit_growth_bayes() and was skipped: %s",
                       m, conditionMessage(attr(fit, "condition"))))
    } else {
      results[[m]] <- fit
    }
  }
  class(results) <- "multimodel_growth_bayes"
  results
}

# ---- S3 methods -------------------------------------------------------------

#' Print a growth_fit_bayes object
#'
#' Prints a short summary of a Bayesian growth-model fit: label/formula,
#' error distribution, sample size, number of chains and iterations,
#' WAIC/DIC, and the posterior summary table (mean, median, sd, credible
#' interval, Rhat and effective sample size per parameter).
#'
#' @param x Object of class \code{"growth_fit_bayes"} (from
#'   \code{\link{fit_growth_bayes}}).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' print(fit_b)   # or just `fit_b` at the console
#' }
#' @export
print.growth_fit_bayes <- function(x, ...) {
  cat(sprintf("Bayesian growth model: %s (%s)\n", x$label, x$model))
  cat(sprintf("Formula: %s\n", x$spec$formula))
  cat(sprintf("Error distribution: %s | n = %d | chains = %d | iterations/chain = %d (burn-in = %d)\n",
              x$distribution, x$n, x$n_chains, x$n_iter, x$n_burnin))
  cat(sprintf("WAIC = %.3f | DIC = %.3f\n", x$waic$waic, x$dic$dic))
  cat("Posterior summary:\n")
  y <- x$summary[x$summary$parameter %in% x$parameters, ]
  num_cols <- c("mean", "median", "sd", "ci_lower", "ci_upper", "rhat", "ess")
  for (cc in num_cols) if (cc %in% names(y)) y[[cc]] <- round(y[[cc]], 4)
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

#' Summarise a growth_fit_bayes object
#'
#' Prints everything \code{\link{print.growth_fit_bayes}} prints, plus
#' the per-chain acceptance rate(s), and issues a warning if any
#' parameter's Rhat exceeds 1.1 (a sign the chains may not have converged).
#'
#' @param object Object of class \code{"growth_fit_bayes"}.
#' @param ... Unused.
#' @return \code{object}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' summary(fit_b)
#' }
#' @export
summary.growth_fit_bayes <- function(object, ...) {
  print(object)
  cat(sprintf("Acceptance rate(s): %s\n",
              paste(sprintf("%.2f", object$accept_rates), collapse = ", ")))
  if (any(object$rhat[object$parameters] > 1.1, na.rm = TRUE)) {
    warning("Some Rhat > 1.1: chains may not have converged. Consider increasing ",
            "n_iter/n_burnin, or check plot_trace_bayes().")
  }
  invisible(object)
}

#' Extract coefficients from a growth_fit_bayes object
#'
#' @param object Object of class \code{"growth_fit_bayes"}.
#' @param ... Unused.
#' @return Named numeric vector with the posterior mean of each parameter.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' coef(fit_b)
#' }
#' @export
coef.growth_fit_bayes <- function(object, ...) object$coefficients

#' Predict sizes from a growth_fit_bayes object, with credible bands
#'
#' Unlike \code{\link{predict.growth_fit}}, the credible band here is a
#' genuine joint posterior-predictive summary: the fitted curve is
#' evaluated at every retained posterior sample (not one parameter at a
#' time), and the band is the pointwise envelope of those curves at each
#' requested age - the Bayesian analogue of \code{\link{bootstrap_ci}}'s band.
#'
#' @param object Object of class \code{"growth_fit_bayes"}.
#' @param new_t Vector of ages at which to predict. If \code{NULL}, a
#'   regular grid over the observed range is used.
#' @param n Number of grid points if \code{new_t} is \code{NULL}.
#' @param interval If \code{TRUE}, also returns \code{ci_lower}/\code{ci_upper}
#'   (posterior credible band, using \code{object$level}).
#' @param ... Unused.
#' @return Data frame with columns \code{t}, \code{y_pred} (posterior mean
#'   curve), and (if \code{interval = TRUE}) \code{ci_lower}/\code{ci_upper}.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' predict(fit_b, new_t = c(1, 5, 10))
#' head(predict(fit_b, interval = TRUE))
#' }
#' @export
predict.growth_fit_bayes <- function(object, new_t = NULL, n = 200, interval = FALSE, ...) {
  if (is.null(new_t)) new_t <- seq(min(object$t), max(object$t), length.out = n)
  post_mean <- as.list(object$coefficients)
  y_pred <- do.call(object$fn, c(list(t = new_t), post_mean))
  out <- data.frame(t = new_t, y_pred = as.numeric(y_pred))
  if (isTRUE(interval)) {
    curves <- apply(object$samples, 1, function(p) {
      do.call(object$fn, c(list(t = new_t), as.list(p[object$parameters])))
    })
    alpha <- 1 - object$level
    out$ci_lower <- apply(curves, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE)
    out$ci_upper <- apply(curves, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  }
  out
}

#' Bayesian model comparison table (WAIC / DIC)
#'
#' Bayesian analogue of \code{\link{aic_table}}: compares
#' \code{"growth_fit_bayes"} models by WAIC and/or DIC, with
#' WAIC-weights computed by direct analogy to Akaike weights
#' (\code{exp(-0.5*delta_WAIC)}, renormalised) - a common informal
#' practice, not a formal Bayesian model-probability statement.
#'
#' @param fit_list List of \code{"growth_fit_bayes"} objects (e.g. from
#'   \code{\link{fit_multimodel_bayes}}).
#' @param criterion \code{"waic"} (default) or \code{"dic"}, used for
#'   ranking and weights.
#' @return Data frame (class \code{"bayes_table_growth"}), ordered best to
#'   worst by the chosen criterion.
#' @examples
#' data(growth_data)
#' \donttest{
#' mm_b <- fit_multimodel_bayes(growth_data, "age", "size",
#'                               models = c("von_bertalanffy", "gompertz"),
#'                               n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' bayes_table(mm_b, criterion = "waic")
#' }
#' @export
bayes_table <- function(fit_list, criterion = c("waic", "dic")) {
  criterion <- match.arg(criterion)
  stopifnot(length(fit_list) > 0)

  rows <- lapply(names(fit_list), function(nm) {
    fit <- fit_list[[nm]]
    data.frame(model = nm, label = fit$label, k = fit$k, n = fit$n,
               WAIC = fit$waic$waic, p_waic = fit$waic$p_waic,
               DIC = fit$dic$dic, p_D = fit$dic$p_D,
               max_rhat = suppressWarnings(max(fit$rhat[fit$parameters], na.rm = TRUE)),
               stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)

  criterion_value <- if (criterion == "waic") tbl$WAIC else tbl$DIC
  tbl$delta <- criterion_value - min(criterion_value, na.rm = TRUE)
  tbl$rel_likelihood <- exp(-0.5 * tbl$delta)
  tbl$weight <- tbl$rel_likelihood / sum(tbl$rel_likelihood, na.rm = TRUE)
  tbl <- tbl[order(tbl$delta), ]
  tbl$rank <- seq_len(nrow(tbl))
  rownames(tbl) <- NULL

  class(tbl) <- c("bayes_table_growth", "data.frame")
  attr(tbl, "criterion") <- criterion
  tbl
}

#' Print a bayes_table_growth object
#'
#' Prints the criterion used for ranking, the comparison table with its
#' numeric columns rounded for readability, and (if applicable) a
#' warning that at least one model's max Rhat exceeds 1.1, meaning its
#' WAIC/DIC should be treated with caution.
#'
#' @param x Object of class \code{"bayes_table_growth"} (from
#'   \code{\link{bayes_table}}).
#' @param digits Number of decimal places for the numeric columns
#'   (default 3).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' \donttest{
#' mm_b <- fit_multimodel_bayes(growth_data, "age", "size",
#'                               models = c("von_bertalanffy", "gompertz"),
#'                               n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' print(bayes_table(mm_b))
#' }
#' @export
print.bayes_table_growth <- function(x, digits = 3, ...) {
  cat(sprintf("Bayesian multi-model comparison (criterion: %s)\n", toupper(attr(x, "criterion"))))
  y <- x
  num_cols <- c("WAIC", "p_waic", "DIC", "p_D", "max_rhat", "delta", "rel_likelihood", "weight")
  for (cc in num_cols) if (cc %in% names(y)) y[[cc]] <- round(y[[cc]], digits)
  print.data.frame(y, row.names = FALSE)
  if (any(x$max_rhat > 1.1, na.rm = TRUE)) {
    cat("Warning: at least one model has max Rhat > 1.1 (possible non-convergence) - ",
        "treat its WAIC/DIC with caution.\n", sep = "")
  }
  invisible(x)
}
