#' MMgrowthR: Multi-Model Fitting of Growth Curves
#'
#' \pkg{MMgrowthR} provides a complete workflow to fit, compare, and plot
#' individual growth models (von Bertalanffy, Gompertz, Gompertz-Laird,
#' logistic, Richards, Schnute, Schnute-Richards, and a seasonal von
#' Bertalanffy model, among others - see \code{\link{GROWTH_MODELS}}) under
#' a multi-model approach: fitting by nonlinear least squares, maximum
#' likelihood, or Bayesian inference, confidence/credible intervals via
#' likelihood profiles (with profile plots), bootstrap, or posterior
#' samples, model comparison via AIC/AICc/BIC (or WAIC/DIC) with Akaike
#' weights, model-averaged predictions, growth-curve comparison between
#' groups, results tables, and elegant plots with \pkg{ggplot2}.
#'
#' A bundled practice dataset, \code{\link{growth_data}}, lets every
#' function below be tried immediately via \code{data(growth_data)}.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Load data: \code{data(growth_data)}, or simulate your own:
#'     \code{\link{simulate_growth_data}}.
#'   \item Fit a model: \code{\link{fit_growth}}, or several at once:
#'     \code{\link{fit_multimodel}}.
#'   \item Confidence intervals: \code{\link{profile_ci}} or \code{\link{bootstrap_ci}}.
#'   \item Compare models: \code{\link{aic_table}}.
#'   \item Plot: \code{\link{plot_growth_fit}}, \code{\link{plot_multimodel}},
#'     \code{\link{plot_bootstrap}}, \code{\link{plot_profile}}.
#'   \item Compare growth between groups (sexes, populations, cohorts...):
#'     \code{\link{fit_growth_grouped}}, \code{\link{growth_lrt}},
#'     \code{\link{growth_lrt_groups}}, \code{\link{plot_growth_fit_grouped}}.
#'   \item Report: \code{\link{parameter_table}}, \code{\link{summary_table}},
#'     \code{\link{export_table}}.
#' }
#'
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"
