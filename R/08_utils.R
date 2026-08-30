#' @keywords internal
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

#' Fixed categorical palette for the growth models
#'
#' A categorical palette (one colour per model in \code{\link{GROWTH_MODELS}})
#' chosen for colour-blind differentiation and contrast: assigned in a
#' fixed order to each model, never recycled or dynamically reordered, so
#' that a given model always has the same colour across plots.
#'
#' @examples
#' MODEL_PALETTE["von_bertalanffy"]
#' MODEL_PALETTE
#' @export
MODEL_PALETTE <- c(
  von_bertalanffy   = "#2a78d6",
  gompertz          = "#eb6834",
  gompertz_laird    = "#1baf7a",
  logistic          = "#eda100",
  richards          = "#e87ba4",
  schnute_richards  = "#008300",
  vb_seasonal       = "#4a3aa7",
  schnute           = "#e34948",
  schnute_case2     = "#a0522d",
  schnute_case3     = "#17becf",
  schnute_case4     = "#bcbd22",
  gallucci_quinn    = "#6c5ce7",
  francis_vb        = "#d63384",
  biphasic_growth   = "#20c997",
  persistence       = "#b5651d",
  tanaka            = "#495057",
  linear            = "#5c6bc0",
  power             = "#26a69a",
  exponential       = "#ef5350",
  logarithmic       = "#8d6e63",
  hyperbolic        = "#00838f",
  ricker            = "#c77c02",
  beverton_holt     = "#7cb342",
  gamma             = "#ad1457",
  weibull           = "#3949ab",
  extended_power    = "#e64a19",
  johnson           = "#8e24aa"
)

#' Format a growth_fit object's estimated coefficients as text
#'
#' Builds a compact, comma-separated \code{"name = value"} string with the
#' fitted parameter estimates, used to annotate plots (e.g.
#' \code{\link{plot_growth_fit}}) alongside the model's symbolic formula.
#'
#' @param fit Object of class \code{"growth_fit"}.
#' @param digits Number of decimal places to round each estimate to.
#' @return A single character string.
#' @keywords internal
.format_coefficients <- function(fit, digits = 4) {
  coefs <- round(fit$coefficients, digits)
  paste(sprintf("%s = %s", names(coefs), format(coefs, trim = TRUE, scientific = FALSE)),
        collapse = ", ")
}
