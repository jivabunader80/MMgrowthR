#' @importFrom grDevices adjustcolor
NULL

#' MMgrowthR's own ggplot2 theme
#'
#' Minimalist theme for \pkg{MMgrowthR} plots: light background, muted
#' gridlines, neutral text tones, and emphasis on the data. Built on top
#' of \code{ggplot2::theme_minimal()}.
#'
#' @param base_size Base font size.
#' @return A \code{ggplot2} theme object, addable to any \code{ggplot} with \code{+}.
#' @examples
#' library(ggplot2)
#' data(growth_data)
#' ggplot(growth_data, aes(age, size)) +
#'   geom_point() +
#'   theme_MMgrowthR()
#' @export
theme_MMgrowthR <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "") +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "#fcfcfb", colour = NA),
      panel.background  = ggplot2::element_rect(fill = "#fcfcfb", colour = NA),
      panel.grid.major  = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.35),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(colour = "#c3c2b7", linewidth = 0.4),
      axis.ticks        = ggplot2::element_line(colour = "#c3c2b7", linewidth = 0.4),
      axis.text         = ggplot2::element_text(colour = "#52514e"),
      axis.title        = ggplot2::element_text(colour = "#0b0b0b", face = "plain"),
      plot.title        = ggplot2::element_text(colour = "#0b0b0b", face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle     = ggplot2::element_text(colour = "#52514e", size = ggplot2::rel(0.85)),
      plot.caption      = ggplot2::element_text(colour = "#898781", size = ggplot2::rel(0.7)),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_text(colour = "#0b0b0b"),
      legend.text       = ggplot2::element_text(colour = "#52514e"),
      strip.text        = ggplot2::element_text(colour = "#0b0b0b", face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "#f0efec", colour = NA)
    )
}

#' Build a curve confidence band from the likelihood profile
#'
#' For each parameter, sweeps the fitted growth curve across that
#' parameter's likelihood-profile values that fall within the chi-squared
#' confidence threshold (Wilks' theorem), holding every other parameter at
#' its point estimate, and returns the pointwise envelope (min/max) across
#' every one of those curves plus the fitted curve itself. This is not a
#' true joint/simultaneous confidence region for the curve (the profile
#' only varies one parameter at a time), but it is a standard, easy to
#' compute visual summary of how the curve's uncertainty - as quantified by
#' the likelihood profile of each parameter - propagates to predictions.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param level Confidence level (default 0.95).
#' @param grid_t Vector of ages/times at which to evaluate the band.
#' @return A data frame with columns \code{t}, \code{lower}, \code{upper},
#'   or \code{NULL} (with a warning) if the profile could not be computed.
#' @keywords internal
.profile_band <- function(fit, level, grid_t) {
  trace <- profile_trace(fit, level = level)
  if (is.null(trace)) return(NULL)

  df <- trace$data
  df <- df[is.finite(df$deviance) & df$deviance <= trace$critical_chisq, , drop = FALSE]

  base_params <- as.list(fit$coefficients)
  curves <- list(do.call(fit$fn, c(list(t = grid_t), base_params)))
  for (pname in unique(df$parameter)) {
    sub <- df[df$parameter == pname, , drop = FALSE]
    for (val in sub$focal) {
      p <- base_params
      p[[pname]] <- val
      curves[[length(curves) + 1]] <- do.call(fit$fn, c(list(t = grid_t), p))
    }
  }
  mat <- do.call(cbind, curves)
  data.frame(t = grid_t,
             lower = apply(mat, 1, min, na.rm = TRUE),
             upper = apply(mat, 1, max, na.rm = TRUE))
}

#' Plot the fit of a growth model
#'
#' Scatter plot of the observed data with the fitted curve overlaid and,
#' optionally, a confidence band for the curve: either from a
#' \strong{likelihood profile} (\code{profile_band = TRUE}, via
#' \code{\link{profile_trace}} - see \code{\link{.profile_band}}) or from
#' \strong{bootstrap} replicate predictions (\code{boot}, the result of
#' \code{\link{bootstrap_ci}}).
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param boot Optional object of class \code{"bootstrap_ci_growth"}
#'   (result of \code{\link{bootstrap_ci}}) to draw a bootstrap confidence
#'   band of the curve. Ignored if \code{profile_band = TRUE}.
#' @param profile_band If \code{TRUE}, draws the confidence band from the
#'   model's likelihood profile instead of from \code{boot} (see
#'   \code{\link{.profile_band}} for how the band is built). Takes priority
#'   over \code{boot} when both are supplied. If the profile cannot be
#'   computed (weakly identified parameters, or non-convergence while
#'   profiling - see \code{\link{profile_trace}}), a warning is issued and
#'   no band is drawn; \code{boot}/\code{\link{bootstrap_ci}} is the
#'   recommended fallback in that case.
#' @param level Confidence level used when \code{profile_band = TRUE}
#'   (default 0.95).
#' @param n_pred Number of points in the prediction grid.
#' @param line_colour Colour of the fitted curve (by default, the fixed
#'   colour of the model in \code{\link{MODEL_PALETTE}}).
#' @param point_shape Shape of the observed-data points (ggplot2 shape
#'   code; default \code{16}, a filled circle).
#' @param point_colour Colour of the observed-data points (default
#'   \code{"#52514e"}), independent of \code{line_colour}.
#' @param ribbon_colour Fill colour of the confidence band. By default
#'   (\code{NULL}) it matches \code{line_colour}, but it can be set
#'   independently (e.g. a neutral grey band under a brightly coloured
#'   line).
#' @param ribbon_alpha Opacity of the confidence band fill (default
#'   \code{0.15}).
#' @param show_formula If \code{TRUE}, includes the model's formula as a
#'   subtitle, together with a second line showing the fitted parameter
#'   estimates (via \code{\link{.format_coefficients}}).
#' @param title Plot title (by default, the model's name).
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' plot_growth_fit(fit, profile_band = TRUE)   # likelihood-profile band
#'
#' \donttest{
#' boot <- bootstrap_ci(fit, R = 200, seed = 1)
#' plot_growth_fit(fit, boot = boot,
#'                 point_shape = 17, point_colour = "steelblue",
#'                 line_colour = "black", ribbon_colour = "grey60", ribbon_alpha = 0.3)
#' }
#' @export
plot_growth_fit <- function(fit, boot = NULL, profile_band = FALSE, level = 0.95,
                             n_pred = 200, line_colour = NULL,
                             point_shape = 16, point_colour = "#52514e",
                             ribbon_colour = NULL, ribbon_alpha = 0.15,
                             show_formula = TRUE, title = NULL) {
  stopifnot(inherits(fit, "growth_fit"))
  pred <- stats::predict(fit, n = n_pred)
  col <- line_colour %||% unname(MODEL_PALETTE[fit$model]) %||% "#2a78d6"
  band_col <- ribbon_colour %||% col
  obs_data <- data.frame(t = fit$t, y = fit$y)

  g <- ggplot2::ggplot()
  band_caption <- NULL

  if (isTRUE(profile_band)) {
    if (!is.null(boot)) {
      warning("plot_growth_fit(): both 'boot' and 'profile_band = TRUE' were supplied; ",
              "using the likelihood-profile band and ignoring 'boot'.")
    }
    band_df <- .profile_band(fit, level = level, grid_t = pred$t)
    if (is.null(band_df)) {
      warning("plot_growth_fit(): the likelihood profile could not be computed for model '",
              fit$model, "'; no confidence band drawn. Consider boot = bootstrap_ci(fit) instead.")
    } else {
      g <- g + ggplot2::geom_ribbon(data = band_df,
                                     ggplot2::aes(x = .data$t, ymin = .data$lower, ymax = .data$upper),
                                     fill = band_col, alpha = ribbon_alpha, colour = NA)
      band_caption <- sprintf(" | %.0f%% CI band: likelihood profile", level * 100)
    }
  } else if (!is.null(boot)) {
    stopifnot(inherits(boot, "bootstrap_ci_growth"))
    fit_params <- names(fit$coefficients)
    boot_params <- colnames(boot$replicates)
    if (!identical(sort(fit_params), sort(boot_params))) {
      stop("plot_growth_fit(): 'boot' was computed for a different model (parameters: ",
           paste(boot_params, collapse = ", "), ") than 'fit' (parameters: ",
           paste(fit_params, collapse = ", "), "). Pass a bootstrap_ci() object ",
           "computed for this exact 'fit' (e.g. bootstrap_ci(fit)).")
    }
    grid_t <- pred$t
    band <- apply(boot$replicates, 1, function(p) {
      do.call(fit$fn, c(list(t = grid_t), as.list(p)))
    })
    alpha_b <- 1 - boot$level
    lower <- apply(band, 1, stats::quantile, probs = alpha_b / 2, na.rm = TRUE)
    upper <- apply(band, 1, stats::quantile, probs = 1 - alpha_b / 2, na.rm = TRUE)
    band_df <- data.frame(t = grid_t, lower = lower, upper = upper)
    g <- g + ggplot2::geom_ribbon(data = band_df,
                                   ggplot2::aes(x = .data$t, ymin = .data$lower, ymax = .data$upper),
                                   fill = band_col, alpha = ribbon_alpha, colour = NA)
    band_caption <- sprintf(" | %.0f%% CI band: bootstrap (%s)", boot$level * 100, boot$type)
  }

  g <- g +
    ggplot2::geom_point(data = obs_data, ggplot2::aes(x = .data$t, y = .data$y),
                         colour = point_colour, alpha = 0.55, size = 1.8, shape = point_shape) +
    ggplot2::geom_line(data = pred, ggplot2::aes(x = .data$t, y = .data$y_pred),
                        colour = col, linewidth = 1) +
    ggplot2::labs(
      title = title %||% sprintf("%s model", fit$label),
      subtitle = if (isTRUE(show_formula)) {
        paste(fit$spec$formula, .format_coefficients(fit), sep = "\n")
      } else {
        NULL
      },
      x = "Age", y = "Size",
      caption = paste0(sprintf("Method: %s | n = %d | AIC = %.1f", toupper(fit$method), fit$n, fit$aic),
                        band_caption %||% "")
    ) +
    theme_MMgrowthR()
  g
}

#' Plot a multi-model growth comparison
#'
#' Overlays the fitted curves of several models on the same observed
#' data, using the fixed categorical palette from \code{\link{MODEL_PALETTE}}.
#'
#' @param fit_list List of \code{"growth_fit"} objects (e.g. the result of
#'   \code{\link{fit_multimodel}}).
#' @param show_aic If \code{TRUE} (default), includes each model's AIC in
#'   the legend.
#' @param n_pred Number of prediction-grid points per curve.
#' @param title Plot title.
#' @param legend_position Where to place the legend: \code{"right"}
#'   (default, outside the panel), \code{"bottom"}, \code{"inside"} (placed
#'   inside the plot panel, saving the lateral space that an outside
#'   legend uses - handy with many models), or a numeric vector
#'   \code{c(x, y)} (each in [0, 1], panel-relative coordinates) for exact
#'   custom placement inside the panel. Ignored when \code{facet = TRUE}
#'   (each panel's strip already identifies its model, so no legend is
#'   drawn).
#' @param point_shape Shape of the observed-data points (ggplot2 shape
#'   code; default \code{19}, a solid circle).
#' @param point_colour Colour of the observed-data points (default
#'   \code{"#898781"}).
#' @param facet If \code{TRUE}, draws one panel per model (via
#'   \code{ggplot2::facet_wrap}) instead of overlaying every curve on a
#'   single panel (the default, \code{FALSE}). Handy when there are many
#'   models and the overlaid curves get hard to tell apart. In this mode
#'   every curve is drawn in a single flat colour (\code{facet_colour})
#'   instead of the per-model \code{\link{MODEL_PALETTE}} colours, since
#'   each panel's strip already identifies its model and a legend would be
#'   redundant.
#' @param facet_colour Line colour used for every curve when
#'   \code{facet = TRUE} (default \code{"#2a78d6"}). Ignored when
#'   \code{facet = FALSE} (there, each model keeps its own
#'   \code{\link{MODEL_PALETTE}} colour so the overlaid curves stay
#'   distinguishable).
#' @param facet_ncol Number of facet columns when \code{facet = TRUE}
#'   (passed to \code{ggplot2::facet_wrap}'s \code{ncol}; \code{NULL},
#'   the default, lets ggplot2 choose automatically).
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' plot_multimodel(mm)
#' plot_multimodel(mm, facet = TRUE)
#' @export
plot_multimodel <- function(fit_list, show_aic = TRUE, n_pred = 200,
                             title = "Multi-model growth comparison",
                             legend_position = "right",
                             point_shape = 19, point_colour = "#898781",
                             facet = FALSE, facet_colour = "#2a78d6", facet_ncol = NULL) {
  stopifnot(length(fit_list) > 0)
  first <- fit_list[[1]]
  obs_data <- data.frame(t = first$t, y = first$y)

  curves <- do.call(rbind, lapply(names(fit_list), function(nm) {
    fit <- fit_list[[nm]]
    pr <- stats::predict(fit, n = n_pred)
    lbl <- if (isTRUE(show_aic)) {
      sprintf("%s (AIC=%.1f)", fit$label, fit$aic)
    } else {
      fit$label
    }
    data.frame(t = pr$t, y_pred = pr$y_pred, label = lbl, model_id = nm)
  }))

  label_order <- vapply(names(fit_list), function(nm) {
    fit <- fit_list[[nm]]
    if (isTRUE(show_aic)) sprintf("%s (AIC=%.1f)", fit$label, fit$aic) else fit$label
  }, character(1))
  curves$label <- factor(curves$label, levels = label_order)

  if (isTRUE(facet)) {
    # A single flat colour for every curve: with one model per panel, the
    # per-model palette colouring of the overlaid view serves no purpose
    # and only adds visual noise. obs_data has no 'label' column, so
    # ggplot2 automatically repeats it in every facet panel (the classic
    # "shared background layer" trick). Strip text is shrunk a bit, since
    # panel strips are narrower than the full-width legend labels
    # (especially with show_aic = TRUE, whose "(AIC=...)" suffix is the
    # main source of long labels - set show_aic = FALSE for shorter strip
    # titles if needed).
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

  # Legend placement (only relevant when facet = FALSE). The installed
  # ggplot2 (< 3.5.0) does not have the modern legend.position = "inside" /
  # legend.position.inside API, so an "inside" placement uses the legacy
  # legend.position = c(x, y) numeric form together with legend.justification,
  # with a translucent background so the legend stays legible over the
  # curves. With several models the legend is also shrunk (smaller
  # text/keys) and split into two columns so it hugs a corner instead of
  # covering the middle of the panel.
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

#' Plot the bootstrap distribution of the parameters
#'
#' Histograms (one per parameter) of the bootstrap replicates, with
#' vertical lines for the point estimate and the confidence-interval
#' bounds.
#'
#' @param boot Object of class \code{"bootstrap_ci_growth"} (result of
#'   \code{\link{bootstrap_ci}}).
#' @param bins Number of histogram bins.
#' @param ncol Number of facet columns (passed to
#'   \code{ggplot2::facet_wrap}'s \code{ncol}). \code{NULL} (default) lets
#'   ggplot2 choose automatically.
#' @param bar_colour Fill colour of the histogram bars (default
#'   \code{"#2a78d6"}).
#' @param bar_alpha Opacity of the histogram bars (default \code{0.55}).
#' @param estimate_colour Colour of the vertical line marking the point
#'   estimate (default \code{"#0b0b0b"}).
#' @param ci_colour Colour of the two vertical dashed lines marking the
#'   confidence-interval bounds (default \code{"#e34948"}).
#' @return A \code{ggplot} object (faceted by parameter).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' \donttest{
#' boot <- bootstrap_ci(fit, R = 200, seed = 1)
#' plot_bootstrap(boot, bins = 20, ncol = 2,
#'                 bar_colour = "steelblue", estimate_colour = "black",
#'                 ci_colour = "darkorange")
#' }
#' @export
plot_bootstrap <- function(boot, bins = 30, ncol = NULL,
                            bar_colour = "#2a78d6", bar_alpha = 0.55,
                            estimate_colour = "#0b0b0b", ci_colour = "#e34948") {
  stopifnot(inherits(boot, "bootstrap_ci_growth"))
  long <- as.data.frame(boot$replicates)
  long <- utils::stack(long)
  names(long) <- c("value", "parameter")

  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = bar_colour, alpha = bar_alpha, colour = NA) +
    ggplot2::geom_vline(data = boot$table,
                         ggplot2::aes(xintercept = .data$estimate),
                         colour = estimate_colour, linewidth = 0.7) +
    ggplot2::geom_vline(data = boot$table,
                         ggplot2::aes(xintercept = .data$ci_lower),
                         colour = ci_colour, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(data = boot$table,
                         ggplot2::aes(xintercept = .data$ci_upper),
                         colour = ci_colour, linetype = "dashed", linewidth = 0.5) +
    ggplot2::facet_wrap(~parameter, scales = "free", ncol = ncol) +
    ggplot2::labs(
      title = sprintf("Bootstrap parameter distribution (%s)", boot$label),
      subtitle = sprintf("%s resampling | %d/%d replicates converged | %.0f%% CI",
                          boot$type, boot$R_converged, boot$R_requested, boot$level * 100),
      x = "Parameter value", y = "Frequency"
    ) +
    theme_MMgrowthR()
  g
}
