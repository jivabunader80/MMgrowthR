#' @importFrom stats coef nls predict logLik AIC pchisq setNames
NULL

# ---- Internal: build a group-dispatching growth function ---------------------

#' Build a growth function that dispatches to group-specific parameters
#'
#' Wraps an already reference-age-fixed model function (\code{fn_local},
#' as returned by \code{.build_fixed_fn()}) into a new function of
#' \code{(t, group, ...)}: for each row, it looks up whether each
#' structural parameter is "free" (one value per group, named
#' \code{"<parameter>.<group>"} in \code{...}) or "common" (one shared
#' value, named \code{"<parameter>"}), builds the right parameter list for
#' that row's group, and evaluates \code{fn_local} on the subset of rows
#' belonging to each group.
#'
#' @param fn_local Reference-age-fixed model function, as returned by
#'   \code{.build_fixed_fn()}.
#' @param parameters Character vector of the model's structural parameter
#'   names (e.g. \code{c("Linf", "K", "t0")}).
#' @param free_params Subset of \code{parameters} that vary by group.
#' @param group_levels Character vector of (sanitised) group levels.
#' @return A function of \code{(t, group, ...)}.
#' @keywords internal
.build_grouped_fn <- function(fn_local, parameters, free_params, group_levels) {
  force(fn_local); force(parameters); force(free_params); force(group_levels)
  function(t, group, ...) {
    args <- list(...)
    out <- rep(NA_real_, length(t))
    group_chr <- as.character(group)
    for (lev in group_levels) {
      idx <- group_chr == lev
      if (!any(idx)) next
      p <- vector("list", length(parameters))
      names(p) <- parameters
      for (pn in parameters) {
        p[[pn]] <- if (pn %in% free_params) args[[paste0(pn, ".", lev)]] else args[[pn]]
      }
      out[idx] <- do.call(fn_local, c(list(t = t[idx]), p))
    }
    out
  }
}

#' @keywords internal
.flat_param_names <- function(parameters, free_params, group_levels) {
  unlist(lapply(parameters, function(pn) {
    if (pn %in% free_params) paste0(pn, ".", group_levels) else pn
  }), use.names = FALSE)
}

# ---- Main interface: fit_growth_grouped ---------------------------------------

#' Fit a growth model with group-specific parameters
#'
#' Fits any of the models in \code{\link{GROWTH_MODELS}} to age-size data
#' from several groups (e.g. sexes, populations, cohorts, years) at once,
#' letting some parameters be estimated \strong{separately per group}
#' (\code{free}) while the rest are held \strong{common} (shared) across
#' groups - a single combined fit rather than one independent
#' \code{\link{fit_growth}} call per group. This is the model-fitting side
#' of the classic growth-comparison framework of Kimura (1980) and Cerrato
#' (1990): \code{free = "all"} reproduces their fully general model
#' (\eqn{\Omega}, every parameter free per group - equivalent to fitting
#' each group separately, but as one combined object with one combined
#' log-likelihood), \code{free = NULL} reproduces their fully common model
#' (\eqn{H_1}, one shared curve, mathematically identical to an ordinary
#' \code{\link{fit_growth}} ignoring group), and any other subset of
#' \code{free} gives an intermediate hypothesis (e.g. only \code{K} shared,
#' \code{Linf}/\code{t0} free). Compare any two nested fits with
#' \code{\link{growth_lrt}}, or the classic full-vs-common test directly
#' with \code{\link{growth_lrt_groups}}.
#'
#' Reference ages for models that need them (\code{t1}/\code{t2} for the
#' Schnute cases, \code{t1}/\code{t2}/\code{t3} for \code{francis_vb}) are
#' fixed once from the \strong{pooled} age range across all groups (as in
#' \code{\link{fit_growth}}), not separately per group, so that groups
#' remain on a directly comparable age scale.
#'
#' @param data Data frame with the data.
#' @param age Name (character) of the age column in \code{data}.
#' @param size Name (character) of the size/weight column in \code{data}.
#' @param group Name (character) of the grouping column in \code{data}
#'   (coerced to character internally; must have at least 2 distinct levels).
#' @param model Name of the model to fit (see \code{names(GROWTH_MODELS)}).
#' @param free Which structural parameters are estimated separately per
#'   group: \code{"all"} (default - every parameter free, the general/
#'   \eqn{\Omega} model), \code{NULL} (no parameter free - a single common
#'   curve, the \eqn{H_1} model), or a character vector naming a subset of
#'   the model's parameters (see \code{GROWTH_MODELS[[model]]$parameters}).
#'   Parameters not named in \code{free} are shared (common) across groups.
#' @param method \code{"ls"} (least squares) or \code{"mle"} (maximum likelihood).
#' @param distribution For \code{method = "mle"}: \code{"normal"} or
#'   \code{"lognormal"}. Ignored if \code{method = "ls"}.
#' @param start Optional named list of starting values, one entry per
#'   flattened parameter name (i.e. \code{"<parameter>"} for common ones,
#'   \code{"<parameter>.<group>"} for free ones - see
#'   \code{names(coef(fit))} after a first attempt with the default
#'   \code{start = NULL} to see the exact names expected). If \code{NULL}
#'   (default), computed automatically: per-group heuristics
#'   (\code{\link{initial_values}} on each group's own data) for free
#'   parameters, and a pooled heuristic for common ones.
#' @param control_ls Control list for \code{minpack.lm::nls.lm.control}.
#' @param control_mle Control list for \code{optim} (via \code{bbmle::mle2}).
#'   With several free-by-group parameters the default \code{maxit} may not
#'   be enough for \code{optim} to fully converge; if you see a
#'   "convergence failure" warning, try e.g.
#'   \code{control_mle = list(maxit = 20000, reltol = 1e-10)}, or use
#'   \code{method = "ls"} instead.
#' @return An object of class \code{"growth_fit_grouped"}, structurally
#'   analogous to \code{"growth_fit"} (same \code{logLik}, \code{k},
#'   \code{aic}/\code{aicc}/\code{bic} fields, so it can be dropped
#'   straight into \code{\link{aic_table}} alongside plain
#'   \code{"growth_fit"} objects to compare hypotheses by information
#'   criterion instead of/alongside a formal LRT).
#' @references
#' Kimura, D.K. (1980). Likelihood methods for the von Bertalanffy growth
#' curve. Fishery Bulletin, 77(4), 765-776.
#'
#' Cerrato, R.M. (1990). Interpretable statistical tests for growth
#' comparisons using parameters in the von Bertalanffy equation. Canadian
#' Journal of Fisheries and Aquatic Sciences, 47(7), 1416-1426.
#' @examples
#' data(growth_data)
#'
#' # Omega: every parameter (Linf, K, t0) free per group
#' fit_full <- fit_growth_grouped(growth_data, "age", "size", "group",
#'                                 model = "von_bertalanffy", free = "all")
#' fit_full
#'
#' # H1: a single common curve (no parameter free)
#' fit_common <- fit_growth_grouped(growth_data, "age", "size", "group",
#'                                   model = "von_bertalanffy", free = NULL)
#'
#' # Intermediate hypothesis: only K held common, Linf/t0 free per group
#' fit_k_common <- fit_growth_grouped(growth_data, "age", "size", "group",
#'                                     model = "von_bertalanffy", free = c("Linf", "t0"))
#'
#' # Compare all three by AIC, exactly like plain growth_fit objects:
#' aic_table(list(full = fit_full, k_common = fit_k_common, common = fit_common))
#' @export
fit_growth_grouped <- function(data, age, size, group, model = "von_bertalanffy",
                                free = "all",
                                method = c("ls", "mle"),
                                distribution = c("normal", "lognormal"),
                                start = NULL, control_ls = NULL, control_mle = NULL) {
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  stopifnot(is.data.frame(data), age %in% names(data), size %in% names(data),
            group %in% names(data))
  stopifnot(model %in% names(GROWTH_MODELS))

  t <- data[[age]]
  y <- data[[size]]
  g_raw <- data[[group]]
  ok <- is.finite(t) & is.finite(y) & !is.na(g_raw)
  if (any(!ok)) warning(sum(!ok), " observations with NA/Inf age, size, or group were excluded.")
  t <- t[ok]; y <- y[ok]; g_raw <- g_raw[ok]

  group_chr <- as.character(g_raw)
  group_levels <- sort(unique(group_chr))
  if (length(group_levels) < 2) {
    stop("fit_growth_grouped(): 'group' must have at least 2 distinct levels ",
         "(found ", length(group_levels), "). Use fit_growth() for a single group.")
  }
  # Sanitise group levels into syntactically valid parameter-name suffixes
  # (e.g. "Linf.Female"); level_lookup keeps the mapping back to the
  # original (possibly non-syntactic) labels for display/predict().
  safe_levels <- make.names(group_levels, unique = TRUE)
  level_lookup <- stats::setNames(group_levels, safe_levels)
  group_safe <- safe_levels[match(group_chr, group_levels)]

  spec <- GROWTH_MODELS[[model]]
  parameters <- spec$parameters

  if (identical(free, "all")) {
    free_params <- parameters
  } else if (is.null(free) || length(free) == 0) {
    free_params <- character(0)
  } else {
    unknown <- setdiff(free, parameters)
    if (length(unknown) > 0) {
      stop("fit_growth_grouped(): 'free' includes name(s) not among the parameters ",
           "of model '", model, "' (", paste(unknown, collapse = ", "), "). Model '",
           model, "' has parameters: ", paste(parameters, collapse = ", "), ".")
    }
    free_params <- free
  }
  common_params <- setdiff(parameters, free_params)

  built <- .build_fixed_fn(spec, t)
  fn_local <- built$fn_local
  fixed_extra <- built$fixed_extra

  grouped_fn <- .build_grouped_fn(fn_local, parameters, free_params, safe_levels)
  flat_params <- .flat_param_names(parameters, free_params, safe_levels)

  if (is.null(start)) {
    pooled_iv <- initial_values(t, y, model = model)
    start <- list()
    for (pn in parameters) {
      if (pn %in% free_params) {
        for (lev in safe_levels) {
          idx <- group_safe == lev
          sub_iv <- tryCatch(initial_values(t[idx], y[idx], model = model),
                              error = function(e) pooled_iv)
          val <- sub_iv[[pn]]
          if (is.null(val) || !is.finite(val)) val <- pooled_iv[[pn]]
          start[[paste0(pn, ".", lev)]] <- val
        }
      } else {
        start[[pn]] <- pooled_iv[[pn]]
      }
    }
  } else {
    missing_start <- setdiff(flat_params, names(start))
    if (length(missing_start) > 0) {
      stop("fit_growth_grouped(): 'start' is missing value(s) for: ",
           paste(missing_start, collapse = ", "), ". With this 'free' setting, ",
           "'start' must be a named list covering exactly: ",
           paste(flat_params, collapse = ", "), ".")
    }
  }

  core <- .fit_core(y = y, predictors = list(t = t, group = group_safe),
                     fn = grouped_fn, parameters = flat_params,
                     method = method, distribution = distribution, start = start,
                     control_ls = control_ls, control_mle = control_mle)

  fit_obj <- core$fit_obj
  n <- length(y)
  if (core$method == "ls") {
    ll <- as.numeric(stats::logLik(fit_obj))
    coefs <- stats::coef(fit_obj)
    sigma <- summary(fit_obj)$sigma
  } else {
    ll <- -fit_obj@min
    coefs <- bbmle::coef(fit_obj)[flat_params]
    sigma <- exp(bbmle::coef(fit_obj)[["log_sigma"]])
  }
  k <- length(flat_params)
  aic <- -2 * ll + 2 * k
  aicc <- if (n - k - 1 > 0) aic + (2 * k * (k + 1)) / (n - k - 1) else NA_real_
  bic <- -2 * ll + k * log(n)

  out <- list(
    model = model, label = spec$label, method = core$method,
    distribution = core$distribution, fit_obj = fit_obj,
    coefficients = coefs, sigma = sigma, logLik = ll, k = k, n = n,
    aic = aic, aicc = aicc, bic = bic,
    t = t, y = y, group = group_chr, group_safe = group_safe,
    group_levels = group_levels, safe_levels = safe_levels, level_lookup = level_lookup,
    free_params = free_params, common_params = common_params,
    age = age, size = size, group_col = group,
    data = data, spec = spec, fixed_extra = fixed_extra,
    fn = grouped_fn, fn_local = fn_local
  )
  class(out) <- "growth_fit_grouped"
  out
}

# ---- S3 methods ----------------------------------------------------------------

#' Print a growth_fit_grouped object
#'
#' Prints label/formula, the grouping column and its levels, which
#' parameters are common vs. free-by-group, fit statistics
#' (logLik/AIC/AICc/BIC), and the flattened coefficient vector.
#'
#' @param x Object of class \code{"growth_fit_grouped"} (from
#'   \code{\link{fit_growth_grouped}}).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' print(fit)   # or just `fit` at the console
#' @export
print.growth_fit_grouped <- function(x, ...) {
  cat(sprintf("Grouped growth model: %s (%s)\n", x$label, x$model))
  cat(sprintf("Formula: %s\n", x$spec$formula))
  cat(sprintf("Group: '%s' (%d levels: %s)\n", x$group_col, length(x$group_levels),
              paste(x$group_levels, collapse = ", ")))
  cat(sprintf("Method: %s | Error distribution: %s\n", toupper(x$method), x$distribution))
  cat(sprintf("n = %d | estimated parameters = %d | logLik = %.3f\n", x$n, x$k, x$logLik))
  cat(sprintf("AIC = %.3f | AICc = %.3f | BIC = %.3f\n", x$aic, x$aicc, x$bic))
  cat("Common (shared across groups):",
      if (length(x$common_params)) paste(x$common_params, collapse = ", ") else "(none)", "\n")
  cat("Free (group-specific):",
      if (length(x$free_params)) paste(x$free_params, collapse = ", ") else "(none)", "\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

#' Summarise a growth_fit_grouped object
#'
#' Prints everything \code{\link{print.growth_fit_grouped}} prints, plus
#' the residual standard deviation (\code{sigma}).
#'
#' @param object Object of class \code{"growth_fit_grouped"}.
#' @param ... Unused.
#' @return \code{object}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' summary(fit)
#' @export
summary.growth_fit_grouped <- function(object, ...) {
  print(object)
  cat(sprintf("Sigma (residual standard deviation): %.4f\n", object$sigma))
  invisible(object)
}

#' Extract coefficients from a growth_fit_grouped object
#'
#' @param object Object of class \code{"growth_fit_grouped"}.
#' @param ... Unused.
#' @return Named numeric vector: common parameters keep their plain name
#'   (e.g. \code{"K"}); free parameters are named \code{"<parameter>.<group>"}
#'   (e.g. \code{"Linf.A"}, \code{"Linf.B"}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' coef(fit)
#' @export
coef.growth_fit_grouped <- function(object, ...) object$coefficients

#' Extract the log-likelihood from a growth_fit_grouped object
#'
#' @param object Object of class \code{"growth_fit_grouped"}.
#' @param ... Unused.
#' @return An object of class \code{"logLik"}, with the total number of
#'   estimated parameters (common + free x groups) as its \code{"df"} attribute.
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' logLik(fit)
#' @export
logLik.growth_fit_grouped <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$k
  class(val) <- "logLik"
  val
}

#' Extract (or recompute) the AIC of a growth_fit_grouped object
#'
#' @param object Object of class \code{"growth_fit_grouped"}.
#' @param ... Unused.
#' @param k Penalty per estimated parameter (default 2, the classic AIC).
#' @return Numeric AIC value.
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' AIC(fit)
#' @export
AIC.growth_fit_grouped <- function(object, ..., k = 2) {
  if (k == 2) object$aic else -2 * object$logLik + k * object$k
}

#' Predict sizes from a growth_fit_grouped object, per group
#'
#' @param object Object of class \code{"growth_fit_grouped"}.
#' @param new_t Vector of ages at which to predict. If \code{NULL}, a
#'   regular grid over the observed range is used.
#' @param group Optional character vector restricting the prediction to a
#'   subset of the fitted group levels. By default (\code{NULL}), predicts
#'   for every group.
#' @param n Number of grid points if \code{new_t} is \code{NULL}.
#' @param ... Unused.
#' @return Data frame with columns \code{t}, \code{group}, and \code{y_pred}
#'   (one row per age x group combination).
#' @examples
#' data(growth_data)
#' fit <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' predict(fit, new_t = c(1, 5, 10))
#' head(predict(fit, group = "Male"))
#' @export
predict.growth_fit_grouped <- function(object, new_t = NULL, group = NULL, n = 200, ...) {
  all_levels <- object$group_levels
  levels_use <- if (is.null(group)) {
    all_levels
  } else {
    unknown <- setdiff(group, all_levels)
    if (length(unknown) > 0) {
      stop("predict.growth_fit_grouped(): 'group' includes level(s) not in the fit (",
           paste(unknown, collapse = ", "), "). Available levels: ",
           paste(all_levels, collapse = ", "), ".")
    }
    group
  }
  if (is.null(new_t)) new_t <- seq(min(object$t), max(object$t), length.out = n)
  parameters <- object$spec$parameters
  coefs <- object$coefficients

  rows <- lapply(levels_use, function(lev) {
    safe_lev <- object$safe_levels[match(lev, object$group_levels)]
    p <- list()
    for (pn in parameters) {
      p[[pn]] <- if (pn %in% object$free_params) {
        unname(coefs[[paste0(pn, ".", safe_lev)]])
      } else {
        unname(coefs[[pn]])
      }
    }
    y_pred <- do.call(object$fn_local, c(list(t = new_t), p))
    data.frame(t = new_t, group = lev, y_pred = as.numeric(y_pred), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- Likelihood-ratio test (Kimura 1980; Cerrato 1990) ------------------------

#' Likelihood-ratio test comparing two nested growth-model fits
#'
#' General-purpose likelihood-ratio test for comparing any two
#' \strong{nested} growth-model fits on the same data - most commonly a
#' \code{\link{fit_growth_grouped}} model against a more constrained one
#' (fewer free-by-group parameters), following Kimura (1980) and Cerrato
#' (1990). \eqn{LR = 2(\ell_{full} - \ell_{reduced})} is compared to a
#' \eqn{\chi^2} distribution with \code{df = k_full - k_reduced} degrees
#' of freedom (Wilks' theorem).
#'
#' For the single most common use - testing whether a fully common curve
#' fits as well as fully separate per-group curves - see the convenience
#' wrapper \code{\link{growth_lrt_groups}}, which fits both extremes for
#' you in one call.
#'
#' @param fit_full The more general/flexible fit (more estimated
#'   parameters) - e.g. a \code{\link{fit_growth_grouped}} result with
#'   more \code{free} parameters.
#' @param fit_reduced The more constrained fit, nested within
#'   \code{fit_full} (fewer estimated parameters) - e.g. the same model
#'   with fewer/no \code{free} parameters, or a plain \code{\link{fit_growth}}.
#' @return An object of class \code{"growth_lrt"}: a list with \code{LR}
#'   (the test statistic), \code{df}, \code{p_value}, and the two fits'
#'   \code{logLik}/\code{k} values and labels.
#' @references
#' Kimura, D.K. (1980). Likelihood methods for the von Bertalanffy growth
#' curve. Fishery Bulletin, 77(4), 765-776.
#'
#' Cerrato, R.M. (1990). Interpretable statistical tests for growth
#' comparisons using parameters in the von Bertalanffy equation. Canadian
#' Journal of Fisheries and Aquatic Sciences, 47(7), 1416-1426.
#' @examples
#' data(growth_data)
#' fit_full <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' fit_common <- fit_growth_grouped(growth_data, "age", "size", "group", free = NULL)
#' growth_lrt(fit_full, fit_common)
#' @export
growth_lrt <- function(fit_full, fit_reduced) {
  ok_classes <- c("growth_fit", "growth_fit_grouped")
  stopifnot(inherits(fit_full, ok_classes), inherits(fit_reduced, ok_classes))
  if (fit_full$n != fit_reduced$n) {
    stop("growth_lrt(): 'fit_full' and 'fit_reduced' were fit to different numbers ",
         "of observations (", fit_full$n, " vs ", fit_reduced$n, "); they must come ",
         "from the same data for the test to be valid.")
  }
  df <- fit_full$k - fit_reduced$k
  if (df <= 0) {
    stop("growth_lrt(): 'fit_full' must have MORE estimated parameters than ",
         "'fit_reduced' (got k_full = ", fit_full$k, ", k_reduced = ", fit_reduced$k,
         "). Check the order of the arguments - the more general/flexible model ",
         "goes first.")
  }
  LR <- 2 * (fit_full$logLik - fit_reduced$logLik)
  if (LR < -1e-6) {
    warning("growth_lrt(): the reduced model has HIGHER log-likelihood than the ",
            "full model (LR = ", round(LR, 4), " < 0), which should not happen if ",
            "'fit_reduced' is genuinely nested in 'fit_full' - check that both fits ",
            "converged (try different starting values). The LR statistic is floored ",
            "at 0 for the test.")
  }
  LR <- max(LR, 0)
  p_value <- stats::pchisq(LR, df = df, lower.tail = FALSE)

  .lrt_label <- function(fit) {
    if (inherits(fit, "growth_fit_grouped")) {
      sprintf("%s [free: %s]", fit$label,
              if (length(fit$free_params)) paste(fit$free_params, collapse = ", ") else "none")
    } else {
      fit$label
    }
  }

  out <- list(LR = LR, df = df, p_value = p_value,
              logLik_full = fit_full$logLik, logLik_reduced = fit_reduced$logLik,
              k_full = fit_full$k, k_reduced = fit_reduced$k,
              label_full = .lrt_label(fit_full), label_reduced = .lrt_label(fit_reduced))
  class(out) <- "growth_lrt"
  out
}

#' Print a growth_lrt object
#'
#' @param x Object of class \code{"growth_lrt"} (from \code{\link{growth_lrt}}
#'   or \code{\link{growth_lrt_groups}}).
#' @param digits Number of decimal places for the numeric output (default 4).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' print(growth_lrt_groups(growth_data, "age", "size", "group"))
#' @export
print.growth_lrt <- function(x, digits = 4, ...) {
  cat("Likelihood-ratio test (Kimura 1980; Cerrato 1990)\n")
  cat(sprintf("Full model:    %s (logLik = %.3f, k = %d)\n", x$label_full, x$logLik_full, x$k_full))
  cat(sprintf("Reduced model: %s (logLik = %.3f, k = %d)\n", x$label_reduced, x$logLik_reduced, x$k_reduced))
  cat(sprintf("LR = %s, df = %d, p-value = %s\n",
              format(round(x$LR, digits), nsmall = digits),
              x$df, format(signif(x$p_value, digits))))
  verdict <- if (x$p_value < 0.05) {
    "Reject H0 (alpha = 0.05): the reduced (more constrained) model fits significantly worse - the full model's extra parameters are supported by the data."
  } else {
    "Do not reject H0 (alpha = 0.05): no significant evidence that the full model's extra parameters improve the fit - the reduced (more parsimonious) model is preferred."
  }
  cat(verdict, "\n")
  invisible(x)
}

#' Classic likelihood-ratio test for a common growth curve across groups
#'
#' Convenience wrapper for the single most common use of
#' \code{\link{growth_lrt}}: testing whether growth can be described by
#' \strong{one common curve} across groups, against the fully general
#' alternative of \strong{completely separate curves per group} - the
#' classic test of Kimura (1980), refined by Cerrato (1990). Internally
#' fits both extremes with \code{\link{fit_growth_grouped}}
#' (\code{free = "all"} and \code{free = NULL}) and compares them with
#' \code{\link{growth_lrt}}.
#'
#' To test an \strong{intermediate} hypothesis (e.g. only \code{K} shared
#' between groups, \code{Linf}/\code{t0} free), fit the two
#' \code{\link{fit_growth_grouped}} models yourself with the appropriate
#' \code{free} arguments and compare them directly with
#' \code{\link{growth_lrt}}.
#'
#' @inheritParams fit_growth_grouped
#' @param ... Additional arguments passed to \code{\link{fit_growth_grouped}}
#'   (e.g. \code{control_ls}).
#' @return An object of class \code{c("growth_lrt_groups", "growth_lrt")}
#'   (prints like \code{\link{growth_lrt}}), with the two underlying
#'   \code{\link{fit_growth_grouped}} fits attached as \code{$fit_full}
#'   and \code{$fit_reduced} for further inspection/plotting.
#' @references
#' Kimura, D.K. (1980). Likelihood methods for the von Bertalanffy growth
#' curve. Fishery Bulletin, 77(4), 765-776.
#'
#' Cerrato, R.M. (1990). Interpretable statistical tests for growth
#' comparisons using parameters in the von Bertalanffy equation. Canadian
#' Journal of Fisheries and Aquatic Sciences, 47(7), 1416-1426.
#' @examples
#' data(growth_data)
#' test <- growth_lrt_groups(growth_data, "age", "size", "group",
#'                            model = "von_bertalanffy")
#' test
#' plot_growth_fit_grouped(test$fit_full)
#' @export
growth_lrt_groups <- function(data, age, size, group, model = "von_bertalanffy",
                               method = c("ls", "mle"),
                               distribution = c("normal", "lognormal"), ...) {
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  fit_full <- fit_growth_grouped(data, age, size, group, model = model, free = "all",
                                  method = method, distribution = distribution, ...)
  fit_reduced <- fit_growth_grouped(data, age, size, group, model = model, free = NULL,
                                     method = method, distribution = distribution, ...)
  result <- growth_lrt(fit_full, fit_reduced)
  result$fit_full <- fit_full
  result$fit_reduced <- fit_reduced
  class(result) <- c("growth_lrt_groups", "growth_lrt")
  result
}
