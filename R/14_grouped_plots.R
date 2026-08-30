NULL

#' Plot a grouped growth-model fit
#'
#' Scatter plot of the observed data (coloured by group) with one fitted
#' curve per group overlaid, from a \code{\link{fit_growth_grouped}} result.
#' A quick visual companion to \code{\link{growth_lrt}}/
#' \code{\link{growth_lrt_groups}}: common parameters produce curves that
#' largely differ only where the free parameters let them, while
#' \code{free = NULL} collapses every group onto a single shared curve.
#'
#' @param fit Object of class \code{"growth_fit_grouped"} (from
#'   \code{\link{fit_growth_grouped}}).
#' @param n_pred Number of prediction-grid points per curve.
#' @param colours Optional named character vector of colours, one per group
#'   level (names must match \code{fit$group_levels}). By default
#'   (\code{NULL}), a built-in categorical palette is used, recycled if
#'   there are more groups than palette colours.
#' @param point_alpha Opacity of the observed-data points (default 0.55).
#' @param point_size Size of the observed-data points (default 1.8).
#' @param line_width Width of the fitted curves (default 1).
#' @param show_formula If \code{TRUE} (default), includes the model's
#'   formula and which parameters are common vs. free-by-group as a
#'   subtitle.
#' @param title Plot title (by default, the model's name).
#' @param legend_position Where to place the group legend (default
#'   \code{"bottom"}, matching \code{\link{theme_MMgrowthR}}).
#' @return A \code{ggplot} object.
#' @examples
#' data(growth_data)
#' fit_full <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' plot_growth_fit_grouped(fit_full)
#'
#' fit_k_common <- fit_growth_grouped(growth_data, "age", "size", "group",
#'                                     free = c("Linf", "t0"))
#' plot_growth_fit_grouped(fit_k_common)
#' @export
plot_growth_fit_grouped <- function(fit, n_pred = 200, colours = NULL,
                                     point_alpha = 0.55, point_size = 1.8,
                                     line_width = 1, show_formula = TRUE,
                                     title = NULL, legend_position = "bottom") {
  stopifnot(inherits(fit, "growth_fit_grouped"))

  obs_data <- data.frame(t = fit$t, y = fit$y, group = fit$group)
  pred <- predict(fit, n = n_pred)

  palette_base <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                     "#008300", "#4a3aa7", "#e34948", "#a0522d", "#17becf")
  levels_use <- fit$group_levels
  if (is.null(colours)) {
    cols <- rep_len(palette_base, length(levels_use))
    names(cols) <- levels_use
  } else {
    missing_cols <- setdiff(levels_use, names(colours))
    if (length(missing_cols) > 0) {
      stop("plot_growth_fit_grouped(): 'colours' is missing entries for group level(s): ",
           paste(missing_cols, collapse = ", "), ".")
    }
    cols <- colours[levels_use]
  }

  free_txt <- if (length(fit$free_params)) paste(fit$free_params, collapse = ", ") else "(none)"
  common_txt <- if (length(fit$common_params)) paste(fit$common_params, collapse = ", ") else "(none)"

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = obs_data,
                         ggplot2::aes(x = .data$t, y = .data$y, colour = .data$group),
                         alpha = point_alpha, size = point_size) +
    ggplot2::geom_line(data = pred,
                        ggplot2::aes(x = .data$t, y = .data$y_pred, colour = .data$group),
                        linewidth = line_width) +
    ggplot2::scale_colour_manual(values = cols, name = fit$group_col) +
    ggplot2::labs(
      title = title %||% sprintf("%s model, by %s", fit$label, fit$group_col),
      subtitle = if (isTRUE(show_formula)) {
        sprintf("%s\nCommon: %s | Free by group: %s", fit$spec$formula, common_txt, free_txt)
      } else {
        NULL
      },
      x = "Age", y = "Size",
      caption = sprintf("Method: %s | n = %d | k = %d | AIC = %.1f",
                         toupper(fit$method), fit$n, fit$k, fit$aic)
    ) +
    theme_MMgrowthR() +
    ggplot2::theme(legend.position = legend_position)
  g
}
