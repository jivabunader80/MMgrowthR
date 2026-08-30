#' Plot likelihood profile traces (chi-squared probability + log-likelihood)
#'
#' Plots the likelihood profile trace of each parameter of a fitted growth
#' model in the classic dual-axis style used for profile-likelihood
#' confidence intervals (e.g. for maturity-ogive parameters, Roa-Ureta
#' style): a dashed \strong{chi-squared probability} curve,
#' \eqn{P(\theta) = 1 - F_{\chi^2_1}(D(\theta))}, on the left axis (bounded
#' in [0, 1], equal to 1 at the estimate), and a solid \strong{profile
#' log-likelihood} curve, \eqn{\ell(\theta) = \hat{\ell} - D(\theta)/2}, on
#' the right axis, where \eqn{D(\theta) = 2(\hat{\ell} - \ell(\theta))} is
#' the profile deviance (Wilks' theorem: asymptotically \eqn{\chi^2_1}
#' under the null \eqn{\theta = \theta_0}). A horizontal dotted line marks
#' \eqn{\alpha = 1 - \text{level}}: the chi-squared curve crosses it exactly
#' at the profile confidence bounds, which are also drawn as vertical
#' dotted lines (with a solid vertical line at the point estimate) when
#' available from \code{\link{profile_ci}}.
#'
#' Because the two curves live on genuinely different scales, each
#' parameter gets its own panel with its own correctly matched dual axis
#' (via \code{ggplot2::sec_axis}), and the panels are combined into a
#' single figure with \pkg{patchwork} — this avoids the distortion that
#' would result from forcing every parameter to share one axis pair.
#'
#' Internally uses \code{stats::profile.nls} for models fitted by least
#' squares (\code{method = "ls"}), or \code{bbmle::profile} for models
#' fitted by maximum likelihood (\code{method = "mle"}).
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param level Confidence level used for the reference lines (default
#'   0.95). Ignored for the vertical CI bounds if \code{ci} is supplied
#'   (its own bounds are used instead) - but still used for the
#'   chi-squared threshold line itself, so keep the two consistent.
#' @param params Optional character vector to restrict the plot to a
#'   subset of parameters. By default, all model parameters are shown.
#' @param ncol Number of panel columns (passed to \code{patchwork::wrap_plots}).
#' @param line_colour Colour of the chi-squared probability curve (by
#'   default, the fixed colour of the model in \code{\link{MODEL_PALETTE}}).
#'   The log-likelihood curve is always drawn in dark ink to keep it
#'   visually distinct.
#' @param ci Optional: an already-computed \code{\link{profile_ci}} result
#'   for this exact \code{fit} (e.g. \code{ci_profile_mle <- profile_ci(fit, ...)}).
#'   When supplied, the vertical CI reference lines are drawn from this
#'   object exactly as computed - including whichever \code{profile_method}
#'   was used - instead of \code{plot_profile()} silently recomputing its
#'   own (possibly different) CI internally. If \code{NULL} (default),
#'   behaves exactly as before (computed internally). \code{ci} is checked
#'   against \code{fit} both by parameter name and by point estimate, so
#'   passing a \code{ci} computed for a different model errors clearly -
#'   even if that other model happens to share the same parameter names
#'   (e.g. \code{von_bertalanffy} and \code{gompertz} both use \code{Linf},
#'   \code{K}, \code{t0}).
#' @return A \code{patchwork} object (one dual-axis panel per parameter),
#'   or \code{NULL} (with a warning) if the profile could not be computed.
#'   Plot and save it like a regular \code{ggplot} object (e.g. with
#'   \code{ggplot2::ggsave()}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' plot_profile(fit)
#'
#' # Plot exactly the CI the user already computed and inspected:
#' ci_profile <- profile_ci(fit, level = 0.95)
#' plot_profile(fit, ci = ci_profile)
#' @export
plot_profile <- function(fit, level = 0.95, params = NULL, ncol = NULL,
                          line_colour = NULL, ci = NULL) {
  stopifnot(inherits(fit, "growth_fit"))
  if (!is.null(params)) {
    fit_params <- names(fit$coefficients)
    unknown <- setdiff(params, fit_params)
    if (length(unknown) > 0) {
      stop("plot_profile(): 'params' includes name(s) not among the parameters of 'fit' (",
           paste(unknown, collapse = ", "), "). Model '", fit$model, "' (", fit$label,
           ") has parameters: ", paste(fit_params, collapse = ", "), ".")
    }
  }
  trace <- profile_trace(fit, which = params, level = level, ci = ci)
  if (is.null(trace)) {
    warning("plot_profile(): no profile trace available for model '", fit$model, "'.")
    return(NULL)
  }

  df <- trace$data
  ci <- trace$ci
  col <- line_colour %||% unname(MODEL_PALETTE[fit$model]) %||% "#2a78d6"
  alpha_line <- 1 - level

  ll_max <- fit$logLik
  df$loglik <- ll_max - df$deviance / 2
  df$chisq_prob <- 1 - stats::pchisq(df$deviance, df = 1)

  parameters <- unique(df$parameter)
  panels <- lapply(parameters, function(pn) {
    d <- df[df$parameter == pn, , drop = FALSE]
    d <- d[order(d$focal), ]

    rng <- range(d$loglik)
    span <- diff(rng)
    if (!is.finite(span) || span <= 0) span <- 1
    d$loglik_scaled <- (d$loglik - rng[1]) / span

    ci_row <- if (!is.null(ci)) ci[ci$parameter == pn, , drop = FALSE] else NULL

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$focal)) +
      ggplot2::geom_hline(yintercept = alpha_line, colour = "#e34948",
                           linetype = "dotted", linewidth = 0.55)

    if (!is.null(ci_row) && nrow(ci_row) == 1) {
      p <- p +
        ggplot2::geom_vline(xintercept = ci_row$estimate, colour = "#c3c2b7", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = ci_row$ci_lower, colour = "#898781",
                             linetype = "dotted", linewidth = 0.45) +
        ggplot2::geom_vline(xintercept = ci_row$ci_upper, colour = "#898781",
                             linetype = "dotted", linewidth = 0.45)
    }

    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$chisq_prob, linetype = "Chi-squared probability"),
                          colour = col, linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(y = .data$loglik_scaled, linetype = "Log-likelihood"),
                          colour = "#0b0b0b", linewidth = 0.8) +
      ggplot2::scale_y_continuous(
        name = "Chi-squared probability",
        limits = c(0, 1),
        sec.axis = ggplot2::sec_axis(~ . * span + rng[1], name = "Log-likelihood")
      ) +
      ggplot2::scale_linetype_manual(
        name = NULL,
        values = c("Chi-squared probability" = "dashed", "Log-likelihood" = "solid")
      ) +
      ggplot2::labs(title = pn, x = "Parameter value") +
      theme_MMgrowthR()
    p
  })

  combined <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  combined <- combined + patchwork::plot_annotation(
    title = sprintf("Likelihood profile - %s model", fit$label),
    subtitle = sprintf(
      "Method: %s | dotted horizontal line: alpha = %.2f (1 - confidence level) | dotted vertical lines: profile CI",
      toupper(fit$method), alpha_line
    ),
    theme = theme_MMgrowthR()
  )
  combined
}

#' Plot a joint (bivariate) likelihood-profile confidence region
#'
#' Draws the 2D profile-deviance surface of a correlated parameter pair
#' (e.g. \code{Linf} and \code{K} of the von Bertalanffy model) as a
#' heatmap, with the \eqn{\chi^2_{2,\,1-\text{level}}} confidence-region
#' boundary overlaid as a contour line - the two-dimensional analogue of
#' \code{\link{plot_profile}}'s single-parameter chi-squared curve. See
#' \code{\link{profile_ci_bivariate}} for why the joint region (typically
#' an elongated, tilted ellipse-like shape) is not the same as the
#' rectangle implied by the two parameters' independent
#' \code{\link{profile_ci}} intervals, when they are correlated.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param params Character vector of exactly two parameter names to
#'   profile jointly. Required unless \code{region} is supplied (in which
#'   case \code{region$params} is used).
#' @param level Confidence level for the region boundary (default 0.95).
#'   Ignored if \code{region} is supplied (its own level is used instead).
#' @param grid_n Number of grid points per parameter, passed to
#'   \code{\link{profile_ci_bivariate}}. Ignored if \code{region} is supplied.
#' @param range_mult Grid half-width in standard errors, passed to
#'   \code{\link{profile_ci_bivariate}}. Ignored if \code{region} is supplied.
#' @param region Optional: an already-computed \code{\link{profile_ci_bivariate}}
#'   result for this exact \code{fit}. When supplied, it is plotted as-is
#'   instead of being recomputed internally - useful to guarantee that
#'   what gets plotted is exactly the object the user already computed and
#'   inspected (e.g. with a non-default \code{grid_n}), and to avoid paying
#'   for the \code{grid_n^2} refits twice. If \code{NULL} (default),
#'   computed internally via \code{profile_ci_bivariate(fit, params, level,
#'   grid_n, range_mult)}.
#' @param low_colour,high_colour Colours spanning the profile-deviance
#'   heatmap gradient (low deviance = best-supported, high deviance =
#'   least-supported).
#' @return A \code{ggplot} object, or \code{NULL} (with a warning) if the
#'   region could not be computed.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy",
#'                    method = "ls")
#' plot_profile_bivariate(fit = fit, params = c("Linf", "K"), grid_n = 15)
#'
#' # Plot exactly the region the user already computed and inspected:
#' region <- profile_ci_bivariate(fit = fit, params = c("Linf", "K"), grid_n = 15)
#' plot_profile_bivariate(fit = fit, region = region)
#' @export
plot_profile_bivariate <- function(fit, params = NULL, level = 0.95, grid_n = 25,
                                    range_mult = 3.5, region = NULL,
                                    low_colour = "#fff6e0", high_colour = "#1c3f6e") {
  stopifnot(inherits(fit, "growth_fit"))
  if (is.null(region)) {
    if (is.null(params)) {
      stop("plot_profile_bivariate(): supply either 'params' (two parameter names) or ",
           "an already-computed 'region' (from profile_ci_bivariate()).")
    }
    region <- try(profile_ci_bivariate(fit, params = params, level = level,
                                        grid_n = grid_n, range_mult = range_mult),
                  silent = TRUE)
    if (inherits(region, "try-error")) {
      warning("plot_profile_bivariate(): ", conditionMessage(attr(region, "condition")))
      return(NULL)
    }
  } else {
    stopifnot(inherits(region, "profile_ci_bivariate"))
    if (!identical(region$model, fit$model)) {
      warning("plot_profile_bivariate(): 'region' was computed for model '", region$model,
              "' but 'fit' is model '", fit$model, "' - plotting 'region' as-is anyway, ",
              "but double-check this is the intended pair.")
    }
  }

  p1 <- region$params[1]; p2 <- region$params[2]
  d <- region$data

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[p1]], y = .data[[p2]])) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$deviance), interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = .data$deviance), breaks = region$critical_chisq,
                           colour = "#e34948", linewidth = 0.9) +
    ggplot2::geom_point(data = data.frame(x = region$estimate[[1]], y = region$estimate[[2]]),
                         ggplot2::aes(x = .data$x, y = .data$y),
                         inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.1,
                         colour = "#0b0b0b") +
    ggplot2::scale_fill_gradient(
      name = "Profile\ndeviance",
      low = low_colour, high = high_colour,
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      title = sprintf("Joint confidence region - %s model", region$label),
      subtitle = sprintf(
        "%s%% joint region (chi-sq, 2 df, threshold = %.3f) | x = point estimate%s",
        format(region$level * 100, trim = TRUE), region$critical_chisq,
        if (region$n_failed > 0) sprintf(" | %d/%d refits failed", region$n_failed, nrow(d)) else ""
      ),
      x = p1, y = p2
    ) +
    theme_MMgrowthR() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  p
}
