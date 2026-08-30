#' Simulate age-size data from a growth model
#'
#' Generates a synthetic age-size dataset by adding observation error
#' (additive or multiplicative) to the predictions of a growth model. Very
#' useful for testing the package's full workflow and for teaching
#' purposes (e.g. to check that the fit recovers the "true" parameters
#' used in the simulation).
#'
#' @param model Model name (see \code{names(GROWTH_MODELS)}).
#' @param parameters Named list with the "true" model parameters.
#' @param n Number of individuals to simulate (ignored if \code{ages} is supplied).
#' @param ages Optional vector of ages to use instead of generating them randomly.
#' @param age_min,age_max Range of simulated ages (if \code{ages} is \code{NULL}).
#' @param cv Coefficient of variation of the error (used if \code{error_sd} is \code{NULL}).
#' @param error_sd Standard deviation of the additive error (if specified,
#'   takes priority over \code{cv} in the additive case).
#' @param error_type \code{"multiplicative"} (default, more realistic for
#'   sizes that grow in magnitude) or \code{"additive"}.
#' @param seed Random seed for reproducibility (\code{NULL} to not set one).
#' @return Data frame with columns \code{age}, \code{size} (observed), and
#'   \code{mean_size} (expected value under the model, without noise).
#' @examples
#' data <- simulate_growth_data(
#'   model = "von_bertalanffy",
#'   parameters = list(Linf = 80, K = 0.35, t0 = -0.3),
#'   n = 150
#' )
#' head(data)
#' @export
simulate_growth_data <- function(model = "von_bertalanffy", parameters,
                                  n = 150, ages = NULL,
                                  age_min = 0.5, age_max = 12,
                                  cv = 0.06, error_sd = NULL,
                                  error_type = c("multiplicative", "additive"),
                                  seed = 123) {
  error_type <- match.arg(error_type)
  stopifnot(model %in% names(GROWTH_MODELS))
  if (!is.null(seed)) set.seed(seed)

  if (is.null(ages)) ages <- sort(stats::runif(n, age_min, age_max))
  spec <- GROWTH_MODELS[[model]]

  if (isTRUE(spec$needs_t1_t2)) {
    t1 <- min(ages); t2 <- max(ages)
    mean_size <- do.call(spec$fn, c(list(t = ages), parameters[spec$parameters],
                                     list(t1 = t1, t2 = t2)))
  } else if (isTRUE(spec$needs_t1_t3)) {
    t1 <- min(ages); t3 <- max(ages); t2 <- (t1 + t3) / 2
    mean_size <- do.call(spec$fn, c(list(t = ages), parameters[spec$parameters],
                                     list(t1 = t1, t2 = t2, t3 = t3)))
  } else {
    mean_size <- do.call(spec$fn, c(list(t = ages), parameters))
  }

  if (error_type == "multiplicative") {
    obs_size <- mean_size * (1 + stats::rnorm(length(ages), mean = 0, sd = cv))
  } else {
    sd_use <- if (is.null(error_sd)) cv * mean(mean_size, na.rm = TRUE) else error_sd
    obs_size <- mean_size + stats::rnorm(length(ages), mean = 0, sd = sd_use)
  }
  obs_size <- pmax(obs_size, 0.01)

  data.frame(age = ages, size = obs_size, mean_size = mean_size)
}
