# ---- Growth rate: computing and plotting dL/dt from a fitted model ---------

#' Compute growth rate at a vector of ages from a fitted model
#'
#' The growth curve \eqn{L(t)} itself describes size-at-age; the
#' \strong{growth rate} describes how fast size is changing at a given
#' age, \eqn{dL/dt}. This is computed \strong{numerically} (central finite
#' differences of the fitted curve, with automatic one-sided fallback near
#' a model's domain boundary - e.g. \code{power_growth}/\code{gamma_growth}
#' being undefined at \code{t = 0}), so it works identically for
#' \strong{every} model in \code{\link{GROWTH_MODELS}}, including any
#' registered at runtime with \code{\link{add_growth_model}}: no
#' model-specific derivative formula needs to be implemented.
#'
#' Two kinds of rate are available via \code{type}:
#' \itemize{
#'   \item \code{"absolute"} (default): the absolute growth rate,
#'     \eqn{AGR(t) = dL/dt}, in size units per age unit (e.g. cm/year) -
#'     how much size is added per unit of time at age \code{t}.
#'   \item \code{"relative"}: the relative (specific) growth rate,
#'     \eqn{RGR(t) = (1/L(t)) \, dL/dt}, in proportion per age unit - the
#'     instantaneous per-unit-size growth rate, useful for comparing
#'     growth speed across models/individuals of very different absolute
#'     size (multiply by 100 for \%/age-unit).
#'   \item \code{"both"}: computes \emph{both} in a single call (at no
#'     extra cost - \eqn{RGR} is just \eqn{AGR} divided by size, and
#'     \eqn{AGR} is always computed internally regardless of \code{type}),
#'     returned as two separate columns instead of one \code{rate} column
#'     - see Value.
#' }
#'
#' @param fit Object of class \code{"growth_fit"} (from \code{\link{fit_growth}}),
#'   for any model in \code{\link{GROWTH_MODELS}} - built-in or added via
#'   \code{\link{add_growth_model}}.
#' @param t Vector of ages at which to compute the rate. If \code{NULL}
#'   (default), a regular grid of \code{n} points over the observed age
#'   range is used (same default grid as \code{\link{predict.growth_fit}}).
#' @param type \code{"absolute"} (\eqn{dL/dt}, default), \code{"relative"}
#'   (\eqn{(1/L)\,dL/dt}), or \code{"both"} - see Details.
#' @param h Step size for the central finite difference. If \code{NULL}
#'   (default), it is set automatically to 1/2000th of the observed age
#'   range (or of \code{max(abs(fit$t), 1)} if that range is degenerate).
#'   Rarely needs changing; pass a larger value only if \code{fn} is very
#'   noisy/flat and the default step gives numerically unstable rates.
#' @param n Number of grid points used when \code{t} is \code{NULL}.
#' @return A data frame (class \code{c("growth_rate_df", "data.frame")})
#'   with columns \code{t} (age) and \code{size} (the fitted \eqn{L(t)}),
#'   plus: a single \code{rate} column (\eqn{AGR} or \eqn{RGR}, per
#'   \code{type}) when \code{type} is \code{"absolute"} or
#'   \code{"relative"}; or, when \code{type = "both"}, two columns,
#'   \code{rate_absolute} and \code{rate_relative}, instead. The
#'   \code{"type"}, \code{"model"}, and \code{"label"} are attached as
#'   attributes (used by \code{\link{plot_growth_rate}}). If a queried age
#'   falls exactly where the model itself is undefined (e.g. \code{t = 0}
#'   for \code{power_growth}/\code{gamma_growth}/\code{logarithmic_growth}/
#'   \code{extended_power_growth}, which return \code{NA} there by
#'   convention - see \code{\link{power_growth}}), \code{size} and every
#'   rate column are \code{NA} at that age too, exactly as
#'   \code{\link{predict.growth_fit}} already behaves; ages arbitrarily
#'   close to (but not exactly at) such a boundary still get a valid rate
#'   from a one-sided difference.
#' @seealso \code{\link{plot_growth_rate}} to plot the result.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' growth_rate(fit, t = c(0, 1, 2, 5, 10))
#' head(growth_rate(fit, type = "relative"))
#' growth_rate(fit, t = c(1, 5, 10), type = "both")   # both rates, one call
#'
#' # Works the same way for a sigmoid model with an interior inflection
#' # point (peak growth rate at an intermediate age, not at t = 0):
#' fit_g <- fit_growth(growth_data, "age", "size", model = "gompertz")
#' r <- growth_rate(fit_g)
#' r$t[which.max(r$rate)]   # age of fastest growth
#' @export
growth_rate <- function(fit, t = NULL, type = c("absolute", "relative", "both"),
                         h = NULL, n = 200) {
  stopifnot(inherits(fit, "growth_fit"))
  type <- match.arg(type)

  if (is.null(t)) {
    t <- seq(min(fit$t), max(fit$t), length.out = n)
  }
  stopifnot(is.numeric(t), length(t) > 0)

  if (is.null(h)) {
    span <- diff(range(fit$t))
    if (!is.finite(span) || span <= 0) span <- max(abs(fit$t), 1)
    h <- span / 2000
  }
  stopifnot(is.numeric(h), length(h) == 1, is.finite(h), h > 0)

  params <- as.list(fit$coefficients)
  size_t  <- as.numeric(do.call(fit$fn, c(list(t = t),     params)))
  size_hi <- as.numeric(do.call(fit$fn, c(list(t = t + h), params)))
  size_lo <- as.numeric(do.call(fit$fn, c(list(t = t - h), params)))

  # Central difference where possible; fall back to a one-sided difference
  # right at a domain boundary (e.g. t - h < 0 for a model undefined there),
  # so a curve like power_growth still gets a rate near t = 0 instead of NA.
  agr <- (size_hi - size_lo) / (2 * h)
  need_fwd <- !is.finite(agr) & is.finite(size_t) & is.finite(size_hi)
  agr[need_fwd] <- (size_hi[need_fwd] - size_t[need_fwd]) / h
  need_bwd <- !is.finite(agr) & is.finite(size_t) & is.finite(size_lo)
  agr[need_bwd] <- (size_t[need_bwd] - size_lo[need_bwd]) / h

  rgr <- agr / size_t

  out <- if (type == "both") {
    data.frame(t = t, size = size_t, rate_absolute = agr, rate_relative = rgr)
  } else {
    data.frame(t = t, size = size_t, rate = if (type == "absolute") agr else rgr)
  }
  attr(out, "type") <- type
  attr(out, "model") <- fit$model
  attr(out, "label") <- fit$label
  class(out) <- c("growth_rate_df", "data.frame")
  out
}

#' @keywords internal
.growth_rate_panel <- function(t, rate, col, y_lab, panel_title, show_peak,
                                linewidth, linetype, zero_line, zero_line_colour,
                                peak_colour) {
  d <- data.frame(t = t, rate = rate)
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$t, y = .data$rate))

  if (isTRUE(zero_line)) {
    g <- g + ggplot2::geom_hline(yintercept = 0, colour = zero_line_colour, linewidth = 0.4)
  }
  g <- g + ggplot2::geom_line(colour = col, linewidth = linewidth, linetype = linetype)

  if (isTRUE(show_peak)) {
    finite_rows <- is.finite(d$rate)
    if (any(finite_rows)) {
      peak <- d[finite_rows, , drop = FALSE]
      peak <- peak[which.max(peak$rate), , drop = FALSE]
      y_range <- range(d$rate, na.rm = TRUE)
      label_y <- y_range[2] + 0.08 * diff(y_range)
      g <- g +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
        ggplot2::geom_vline(xintercept = peak$t, colour = peak_colour,
                             linetype = "dashed", linewidth = 0.5) +
        ggplot2::annotate("text", x = peak$t, y = label_y,
                           label = sprintf(" peak at t = %.2f", peak$t),
                           hjust = 0, vjust = 0, size = 3.2, colour = peak_colour)
    }
  }

  g + ggplot2::labs(title = panel_title, x = "Age", y = y_lab) + theme_MMgrowthR()
}

#' Plot growth rate against age
#'
#' Draws the growth-rate curve (\eqn{dL/dt}, the relative rate, or both -
#' see \code{\link{growth_rate}}) computed from a fitted growth model, in
#' the same visual style as \code{\link{plot_growth_fit}}.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param t Vector of ages at which to evaluate the rate. If \code{NULL}
#'   (default), a regular grid of \code{n} points over the observed age
#'   range is used.
#' @param type \code{"absolute"} (\eqn{dL/dt}, default), \code{"relative"}
#'   (\eqn{(1/L)\,dL/dt}), or \code{"both"} - see \code{\link{growth_rate}}.
#'   \code{"absolute"}/\code{"relative"} draw a single panel; \code{"both"}
#'   draws \strong{both} rates as two stacked panels (via \pkg{patchwork}),
#'   each on its own correctly scaled y-axis - the two rates live on
#'   different scales/units, so overlaying them on one shared y-axis would
#'   be misleading (the same reasoning \code{\link{plot_profile}} uses for
#'   its own dual-scale curves).
#' @param h Step size for the underlying finite difference; see
#'   \code{\link{growth_rate}}.
#' @param n Number of grid points when \code{t} is \code{NULL}.
#' @param line_colour Colour of the rate curve (by default, the model's
#'   fixed colour in \code{\link{MODEL_PALETTE}}). Applied to both panels
#'   when \code{type = "both"} (the two panels are the same model's two
#'   rates, so they share one colour by default; pass \code{colours}
#'   instead for two different colours).
#' @param colours Only used when \code{type = "both"}: a length-2 character
#'   vector, \code{c(colour_for_absolute, colour_for_relative)}, to give
#'   the two panels different colours instead of the single
#'   \code{line_colour}. \code{NULL} (default) uses \code{line_colour} for
#'   both.
#' @param linewidth Thickness of the rate curve (default \code{1}).
#' @param linetype Line type of the rate curve (default \code{"solid"};
#'   any \code{ggplot2} linetype, e.g. \code{"dashed"}, \code{"dotted"}).
#' @param zero_line If \code{TRUE} (default), draws a thin horizontal
#'   reference line at \code{rate = 0}. Set \code{FALSE} to omit it (e.g.
#'   for a relative rate that never goes negative, where it adds little).
#' @param zero_line_colour Colour of the zero reference line (default a
#'   light neutral grey, \code{"#c3c2b7"}). Ignored if \code{zero_line = FALSE}.
#' @param show_peak If \code{TRUE}, marks the age of maximum rate within
#'   the plotted range with a dashed vertical line and a label. Most
#'   informative for a sigmoid model with a genuine interior inflection
#'   point (e.g. Gompertz, logistic, Richards); for a model whose rate is
#'   monotonic over the plotted range (e.g. von Bertalanffy, whose rate
#'   peaks at the youngest age shown), the mark simply falls at the edge
#'   of the range and is less informative. Default \code{FALSE}. Applied
#'   to each panel independently when \code{type = "both"}.
#' @param peak_colour Colour of the peak marker line/label (default
#'   \code{"#52514e"}, dark neutral grey). Ignored if \code{show_peak = FALSE}.
#' @param title Plot title (by default, the model's label).
#' @param ncol Number of panel columns when \code{type = "both"}
#'   (\code{1}, the default, stacks the two panels vertically so they
#'   share a directly comparable x-axis; use \code{2} to place them
#'   side by side instead). Ignored otherwise.
#' @return A \code{ggplot} object (single panel), or - when
#'   \code{type = "both"} - a \code{patchwork} object (two panels);
#'   both can be printed, further modified with \code{+}/\code{&}, and
#'   saved with \code{ggplot2::ggsave()} in the usual way.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' plot_growth_rate(fit)
#' plot_growth_rate(fit, type = "relative")
#' plot_growth_rate(fit, type = "both")   # both rates, one figure
#'
#' # Customising the look:
#' plot_growth_rate(fit, line_colour = "firebrick", linewidth = 1.5, linetype = "dashed")
#' plot_growth_rate(fit, type = "both", colours = c("steelblue", "darkorange"))
#'
#' fit_g <- fit_growth(growth_data, "age", "size", model = "gompertz")
#' plot_growth_rate(fit_g, show_peak = TRUE, peak_colour = "firebrick")
#' @export
plot_growth_rate <- function(fit, t = NULL, type = c("absolute", "relative", "both"),
                              h = NULL, n = 200, line_colour = NULL, colours = NULL,
                              linewidth = 1, linetype = "solid",
                              zero_line = TRUE, zero_line_colour = "#c3c2b7",
                              show_peak = FALSE, peak_colour = "#52514e",
                              title = NULL, ncol = 1) {
  stopifnot(inherits(fit, "growth_fit"))
  type <- match.arg(type)

  rate_df <- growth_rate(fit, t = t, type = type, h = h, n = n)
  col <- line_colour %||% unname(MODEL_PALETTE[fit$model]) %||% "#2a78d6"
  main_title <- title %||% sprintf("%s growth rate", fit$label)

  if (type != "both") {
    y_lab <- if (type == "absolute") {
      "Growth rate (size / age unit)"
    } else {
      "Relative growth rate (proportion / age unit)"
    }
    g <- .growth_rate_panel(rate_df$t, rate_df$rate, col, y_lab, NULL, show_peak,
                             linewidth, linetype, zero_line, zero_line_colour, peak_colour)
    return(
      g + ggplot2::labs(
        title = main_title,
        caption = sprintf("Type: %s | Method: %s | n = %d",
                           if (type == "absolute") "Absolute (dL/dt)" else "Relative ((1/L) dL/dt)",
                           toupper(fit$method), fit$n)
      )
    )
  }

  if (!is.null(colours)) {
    stopifnot(is.character(colours), length(colours) == 2)
    col_abs <- colours[1]; col_rel <- colours[2]
  } else {
    col_abs <- col; col_rel <- col
  }

  panel_abs <- .growth_rate_panel(rate_df$t, rate_df$rate_absolute, col_abs,
                                   "Absolute (size / age unit)",
                                   "Absolute growth rate (dL/dt)", show_peak,
                                   linewidth, linetype, zero_line, zero_line_colour, peak_colour)
  panel_rel <- .growth_rate_panel(rate_df$t, rate_df$rate_relative, col_rel,
                                   "Relative (proportion / age unit)",
                                   "Relative growth rate ((1/L) dL/dt)", show_peak,
                                   linewidth, linetype, zero_line, zero_line_colour, peak_colour)

  combined <- patchwork::wrap_plots(list(panel_abs, panel_rel), ncol = ncol)
  combined + patchwork::plot_annotation(
    title = main_title,
    subtitle = sprintf("Method: %s | n = %d", toupper(fit$method), fit$n),
    theme = theme_MMgrowthR()
  )
}
