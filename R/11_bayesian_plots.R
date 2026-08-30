#' Plot the fit of a Bayesian growth model, with a posterior credible band
#'
#' Bayesian analogue of \code{\link{plot_growth_fit}}: scatter of the
#' observed data, the posterior-mean curve, and a genuine joint
#' posterior-predictive credible band (via \code{\link{predict.growth_fit_bayes}}).
#'
#' @param fit Object of class \code{"growth_fit_bayes"}.
#' @param n_pred Number of prediction-grid points.
#' @param line_colour Colour of the fitted curve (default: the model's
#'   fixed colour in \code{\link{MODEL_PALETTE}}).
#' @param point_shape,point_colour Aesthetics of the observed-data points.
#' @param ribbon_colour,ribbon_alpha Aesthetics of the credible band
#'   (defaults to \code{line_colour} and \code{0.15}).
#' @param show_formula If \code{TRUE} (default), adds the model's formula
#'   and posterior-mean estimates as a subtitle.
#' @param title Optional plot title.
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' plot_growth_fit_bayes(fit_b)
#' }
#' @export
plot_growth_fit_bayes <- function(fit, n_pred = 200, line_colour = NULL,
                                   point_shape = 16, point_colour = "#52514e",
                                   ribbon_colour = NULL, ribbon_alpha = 0.15,
                                   show_formula = TRUE, title = NULL) {
  stopifnot(inherits(fit, "growth_fit_bayes"))
  pred <- predict(fit, n = n_pred, interval = TRUE)
  col <- line_colour %||% unname(MODEL_PALETTE[fit$model]) %||% "#2a78d6"
  band_col <- ribbon_colour %||% col
  obs_data <- data.frame(t = fit$t, y = fit$y)

  coef_txt <- paste(sprintf("%s = %s", names(fit$coefficients),
                             format(round(fit$coefficients, 4), trim = TRUE, scientific = FALSE)),
                     collapse = ", ")

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pred, ggplot2::aes(x = .data$t, ymin = .data$ci_lower, ymax = .data$ci_upper),
                          fill = band_col, alpha = ribbon_alpha, colour = NA) +
    ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$t, y = .data$y),
                         colour = point_colour, alpha = 0.55, size = 1.8, shape = point_shape) +
    ggplot2::geom_line(data = pred, ggplot2::aes(x = .data$t, y = .data$y_pred),
                        colour = col, linewidth = 1) +
    ggplot2::labs(
      title = title %||% sprintf("%s model (Bayesian)", fit$label),
      subtitle = if (isTRUE(show_formula)) paste(fit$spec$formula, coef_txt, sep = "\n") else NULL,
      x = "Age", y = "Size",
      caption = sprintf("n = %d | %d chains x %d it. (burn-in %d) | WAIC = %.1f | %.0f%% credible band (posterior predictive)",
                         fit$n, fit$n_chains, fit$n_iter, fit$n_burnin, fit$waic$waic, fit$level * 100)
    ) +
    theme_MMgrowthR()
}

#' Plot a multi-model Bayesian growth comparison
#'
#' Bayesian analogue of \code{\link{plot_multimodel}}: overlays the
#' posterior-mean curves of several \code{"growth_fit_bayes"} models on
#' the same observed data, using the fixed categorical palette from
#' \code{\link{MODEL_PALETTE}}. Accepts the same \code{facet}/
#' \code{legend_position} options, with \code{WAIC} shown instead of AIC.
#'
#' @param fit_list List of \code{"growth_fit_bayes"} objects (e.g. the
#'   result of \code{\link{fit_multimodel_bayes}}).
#' @param show_waic If \code{TRUE} (default), includes each model's WAIC
#'   in the legend/facet strip.
#' @param n_pred Number of prediction-grid points per curve.
#' @param title Plot title.
#' @param legend_position Where to place the legend: \code{"right"}
#'   (default), \code{"bottom"}, \code{"inside"}, or a numeric
#'   \code{c(x, y)}. Ignored when \code{facet = TRUE}.
#' @param point_shape,point_colour Aesthetics of the observed-data points.
#' @param facet If \code{TRUE}, draws one panel per model instead of
#'   overlaying every curve on a single panel; every curve is then drawn
#'   in a single flat colour (\code{facet_colour}).
#' @param facet_colour Line colour used for every curve when \code{facet = TRUE}.
#' @param facet_ncol Number of facet columns when \code{facet = TRUE}.
#' @param show_rhat_warning If \code{TRUE} (default), appends a "(Rhat > 1.1)"
#'   marker to the label of any model whose \code{max(fit$rhat) > 1.1},
#'   flagging possible non-convergence directly on the plot.
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' \donttest{
#' mm_b <- fit_multimodel_bayes(growth_data, "age", "size",
#'                               models = c("von_bertalanffy", "gompertz"),
#'                               n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' plot_multimodel_bayes(mm_b)
#' plot_multimodel_bayes(mm_b, facet = TRUE)
#' }
#' @export
plot_multimodel_bayes <- function(fit_list, show_waic = TRUE, n_pred = 200,
                                   title = "Multi-model growth comparison (Bayesian)",
                                   legend_position = "right",
                                   point_shape = 19, point_colour = "#898781",
                                   facet = FALSE, facet_colour = "#2a78d6", facet_ncol = NULL,
                                   show_rhat_warning = TRUE) {
  stopifnot(length(fit_list) > 0)
  first <- fit_list[[1]]
  obs_data <- data.frame(t = first$t, y = first$y)

  make_label <- function(fit) {
    lbl <- fit$label
    if (isTRUE(show_waic)) lbl <- sprintf("%s (WAIC=%.1f)", lbl, fit$waic$waic)
    if (isTRUE(show_rhat_warning) && suppressWarnings(max(fit$rhat[fit$parameters], na.rm = TRUE)) > 1.1) {
      lbl <- paste(lbl, "[Rhat>1.1]")
    }
    lbl
  }

  curves <- do.call(rbind, lapply(names(fit_list), function(nm) {
    fit <- fit_list[[nm]]
    pr <- stats::predict(fit, n = n_pred)
    data.frame(t = pr$t, y_pred = pr$y_pred, label = make_label(fit), model_id = nm)
  }))

  label_order <- vapply(names(fit_list), function(nm) make_label(fit_list[[nm]]), character(1))
  curves$label <- factor(curves$label, levels = label_order)

  if (isTRUE(facet)) {
    g <- ggplot2::ggplot() +
      ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$t, y = .data$y),
                           colour = point_colour, alpha = 0.45, size = 1.6, shape = point_shape) +
      ggplot2::geom_line(data = curves,
                          ggplot2::aes(x = .data$t, y = .data$y_pred, group = .data$label),
                          colour = facet_colour, linewidth = 1) +
      ggplot2::labs(title = title, x = "Age", y = "Size") +
      theme_MMgrowthR() +
      ggplot2::facet_wrap(~label, ncol = facet_ncol) +
      ggplot2::theme(legend.position = "none",
                      strip.text = ggplot2::element_text(size = ggplot2::rel(0.75)))
    return(g)
  }

  colours <- unname(MODEL_PALETTE[names(fit_list)])
  names(colours) <- label_order

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$t, y = .data$y),
                         colour = point_colour, alpha = 0.45, size = 1.6, shape = point_shape) +
    ggplot2::geom_line(data = curves,
                        ggplot2::aes(x = .data$t, y = .data$y_pred, colour = .data$label),
                        linewidth = 1) +
    ggplot2::scale_colour_manual(values = colours, name = "Model") +
    ggplot2::labs(title = title, x = "Age", y = "Size") +
    theme_MMgrowthR()

  n_models <- length(fit_list)
  inside_theme <- function(xy, justification = xy) {
    ggplot2::theme(
      legend.position = xy,
      legend.justification = justification,
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("#fcfcfb", alpha.f = 0.85),
                                                  colour = "#e1e0d9"),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.size = ggplot2::unit(0.7, "lines"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.7)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.spacing.y = ggplot2::unit(0.1, "lines")
    )
  }

  if (identical(legend_position, "inside")) {
    g <- g + inside_theme(c(0.99, 0.01), justification = c(1, 0)) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = if (n_models > 5) 2 else 1))
  } else if (is.numeric(legend_position) && length(legend_position) == 2) {
    g <- g + inside_theme(legend_position) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = if (n_models > 5) 2 else 1))
  } else {
    g <- g + ggplot2::theme(legend.position = legend_position)
  }
  g
}

#' Plot posterior distributions of a Bayesian growth fit's parameters
#'
#' Histograms of the (pooled, post-burn-in) posterior samples for each
#' parameter, with the posterior mean and credible interval bounds marked
#' - the Bayesian analogue of \code{\link{plot_bootstrap}}.
#'
#' @param fit Object of class \code{"growth_fit_bayes"}.
#' @param bins Number of histogram bins.
#' @param ncol Number of facet columns.
#' @param bar_colour,bar_alpha Histogram bar aesthetics.
#' @param estimate_colour,ci_colour Colours of the mean and credible-interval lines.
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' plot_posterior(fit_b)
#' }
#' @export
plot_posterior <- function(fit, bins = 30, ncol = NULL,
                            bar_colour = "#2a78d6", bar_alpha = 0.55,
                            estimate_colour = "#0b0b0b", ci_colour = "#e34948") {
  stopifnot(inherits(fit, "growth_fit_bayes"))
  samp <- fit$samples[, fit$parameters, drop = FALSE]
  df <- do.call(rbind, lapply(fit$parameters, function(pn) {
    data.frame(parameter = pn, value = samp[, pn])
  }))
  df$parameter <- factor(df$parameter, levels = fit$parameters)
  summ <- fit$summary[fit$summary$parameter %in% fit$parameters, ]
  summ$parameter <- factor(summ$parameter, levels = fit$parameters)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = bar_colour, alpha = bar_alpha, colour = NA) +
    ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = .data$mean),
                         colour = estimate_colour, linewidth = 0.7) +
    ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = .data$ci_lower),
                         colour = ci_colour, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = .data$ci_upper),
                         colour = ci_colour, linetype = "dashed", linewidth = 0.5) +
    ggplot2::facet_wrap(~parameter, scales = "free", ncol = ncol) +
    ggplot2::labs(title = sprintf("Posterior distributions - %s (Bayesian)", fit$label),
                  subtitle = sprintf("Solid line: posterior mean | dashed lines: %.0f%% credible interval",
                                      fit$level * 100),
                  x = "Parameter value", y = "Count") +
    theme_MMgrowthR()
}

#' Trace plots for MCMC convergence diagnostics
#'
#' One panel per parameter, one line per chain, showing the sampled value
#' across iterations (post-burn-in only by default). Well-mixed,
#' converged chains look like flat, overlapping "hairy caterpillars" with
#' no trend or separation between chains; visible trends, drift, or
#' chains occupying different regions are signs of non-convergence (also
#' flagged numerically by \eqn{\hat R} in \code{fit$rhat}).
#'
#' @param fit Object of class \code{"growth_fit_bayes"}.
#' @param include_burnin If \code{TRUE}, also plots the (discarded)
#'   burn-in iterations, shaded, so the adaptation/burn-in period can be
#'   inspected too. Requires the original chains to still include it
#'   (always true, since \code{fit$chains} stores only post-burn-in
#'   samples - set \code{include_burnin = FALSE}, the default, unless you
#'   refit keeping the raw chains).
#' @param ncol Number of facet columns.
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' \donttest{
#' fit_b <- fit_growth_bayes(growth_data, "age", "size", model = "von_bertalanffy",
#'                            n_iter = 3000, n_burnin = 800, n_chains = 3, seed = 1)
#' plot_trace_bayes(fit_b)
#' }
#' @export
plot_trace_bayes <- function(fit, include_burnin = FALSE, ncol = NULL) {
  stopifnot(inherits(fit, "growth_fit_bayes"))
  if (isTRUE(include_burnin)) {
    warning("plot_trace_bayes(): fit$chains only stores post-burn-in samples; ",
            "plotting those (include_burnin has no effect).")
  }
  rows <- list()
  for (ch_i in seq_along(fit$chains)) {
    ch <- fit$chains[[ch_i]][, fit$parameters, drop = FALSE]
    for (pn in fit$parameters) {
      rows[[length(rows) + 1]] <- data.frame(
        iteration = seq_len(nrow(ch)), chain = factor(ch_i),
        parameter = pn, value = ch[, pn]
      )
    }
  }
  df <- do.call(rbind, rows)
  df$parameter <- factor(df$parameter, levels = fit$parameters)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, colour = .data$chain)) +
    ggplot2::geom_line(alpha = 0.75, linewidth = 0.35) +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = ncol) +
    ggplot2::labs(title = sprintf("MCMC trace plots (post-burn-in) - %s", fit$label),
                  subtitle = "Well-mixed chains: flat, overlapping, no trend or separation",
                  x = "Iteration (post-burn-in, thinned)", y = "Parameter value", colour = "Chain") +
    theme_MMgrowthR()
}
