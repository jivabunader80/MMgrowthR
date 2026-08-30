#' @title Individual growth model functions
#' @name growth_model_family
#' @description
#' Family of functions implementing the size-at-age equations for the most
#' commonly used age-size growth models in fisheries biology and ecology:
#' von Bertalanffy, Gompertz, Gompertz-Laird, logistic, Richards, Schnute,
#' Schnute-Richards, and a seasonally oscillating von Bertalanffy model
#' (Somers 1988) - plus four classical, non-asymptotic single-predictor
#' curves (linear, power, exponential, logarithmic) useful as simple
#' baselines or for growth phases where an asymptotic model is not
#' appropriate (see \code{\link{linear_growth}} and neighbouring functions),
#' and seven further named curves drawn from the growth-model literature:
#' two more saturating/asymptotic two-parameter forms
#' (\code{\link{hyperbolic_growth}}, Gulland & Holt 1959;
#' \code{\link{beverton_holt_growth}}, Beverton & Holt 1959), two
#' dome-shaped forms mainly suited to data that genuinely decline at old
#' age (\code{\link{ricker_growth}}, Ricker 1954; \code{\link{gamma_growth}},
#' Troynikov & Gorfine 1998), the four-parameter \code{\link{weibull_growth}}
#' (Seber & Wild 1989), the three-parameter \code{\link{extended_power_growth}}
#' (an extension of the plain power model, Mercier et al. 2011), and
#' \code{\link{johnson_growth}} (Ricker 1975), an asymptotic curve without
#' an evident inflection point, alongside von Bertalanffy.
#'
#' All functions are vectorised over \code{t} and return the predicted
#' size. They are used both to simulate data and inside the fitting
#' engines (\code{\link{fit_growth}}).
NULL

#' Von Bertalanffy growth model
#'
#' L(t) = Linf * (1 - exp(-K * (t - t0)))
#'
#' @param t Vector of ages.
#' @param Linf Asymptotic size.
#' @param K Growth coefficient (rate at which Linf is approached).
#' @param t0 Hypothetical age at zero size.
#' @return Vector of predicted sizes.
#' @examples
#' von_bertalanffy(t = 1:10, Linf = 85, K = 0.28, t0 = -0.35)
#'
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' fit
#' @export
von_bertalanffy <- function(t, Linf, K, t0) {
  Linf * (1 - exp(-K * (t - t0)))
}

#' Gompertz growth model
#'
#' L(t) = Linf * exp(-exp(-K * (t - t0)))
#'
#' @param t Vector of ages.
#' @param Linf Asymptotic size.
#' @param K Relative instantaneous growth rate at the inflection point.
#' @param t0 Age at the inflection point of the curve (size = Linf/e).
#' @return Vector of predicted sizes.
#' @examples
#' gompertz(t = 1:10, Linf = 85, K = 0.6, t0 = 1)
#' @export
gompertz <- function(t, Linf, K, t0) {
  Linf * exp(-exp(-K * (t - t0)))
}

#' Gompertz-Laird growth model
#'
#' Alternative parameterisation of Gompertz (Laird 1965), widely used for
#' animal growth curves, expressed in terms of the initial size \code{L0}.
#'
#' L(t) = L0 * exp( (A / k) * (1 - exp(-k * t)) )
#'
#' @param t Vector of ages (or time since birth/hatching).
#' @param L0 Initial size (t = 0).
#' @param A Initial specific growth rate.
#' @param k Exponential decay rate of the growth rate.
#' @return Vector of predicted sizes.
#' @examples
#' gompertz_laird(t = 0:10, L0 = 5, A = 0.9, k = 0.3)
#' @export
gompertz_laird <- function(t, L0, A, k) {
  L0 * exp((A / k) * (1 - exp(-k * t)))
}

#' Logistic growth model
#'
#' L(t) = Linf / (1 + exp(-K * (t - t0)))
#'
#' @param t Vector of ages.
#' @param Linf Asymptotic size.
#' @param K Growth rate.
#' @param t0 Age at the inflection point (size = Linf/2).
#' @return Vector of predicted sizes.
#' @examples
#' logistic_growth(t = 1:10, Linf = 85, K = 0.7, t0 = 2)
#' @export
logistic_growth <- function(t, Linf, K, t0) {
  Linf / (1 + exp(-K * (t - t0)))
}

#' Richards growth model
#'
#' Flexible generalisation that includes von Bertalanffy (b = 1) as a
#' special case, and approaches Gompertz as b -> 0 and the logistic model
#' as b -> -1.
#'
#' L(t) = Linf * (1 - exp(-K * (t - t0)))^(1 / b)
#'
#' @param t Vector of ages.
#' @param Linf Asymptotic size.
#' @param K Growth rate.
#' @param t0 Hypothetical age at zero size.
#' @param b Shape parameter (b > 0).
#' @return Vector of predicted sizes.
#' @examples
#' richards(t = 1:10, Linf = 85, K = 0.3, t0 = -0.5, b = 1.5)
#' @export
richards <- function(t, Linf, K, t0, b) {
  base <- 1 - exp(-K * (t - t0))
  base[base < 0] <- NA_real_
  Linf * base^(1 / b)
}

#' Schnute growth model (general case, 4 parameters)
#'
#' General formulation of Schnute (1981), which avoids indeterminacies
#' when \code{a} or \code{b} approach zero by using the corresponding
#' analytical limits. Requires fixing two reference ages \code{t1} and
#' \code{t2} (typically the minimum and maximum of the observed ages),
#' which correspond to sizes \code{y1} and \code{y2}.
#'
#' @param t Vector of ages.
#' @param a Parameter related to the growth rate (may be 0).
#' @param b Shape parameter (may be 0).
#' @param y1 Estimated size at the reference age \code{t1}.
#' @param y2 Estimated size at the reference age \code{t2}.
#' @param t1 Smaller (fixed, not estimated) reference age.
#' @param t2 Larger (fixed, not estimated) reference age.
#' @return Vector of predicted sizes.
#' @examples
#' schnute(t = 1:10, a = 0.3, b = 1, y1 = 15, y2 = 80, t1 = 1, t2 = 10)
#' @export
schnute <- function(t, a, b, y1, y2, t1, t2) {
  eps <- 1e-6
  if (abs(a) > eps && abs(b) > eps) {
    out <- (y1^b + (y2^b - y1^b) * (1 - exp(-a * (t - t1))) /
              (1 - exp(-a * (t2 - t1))))^(1 / b)
  } else if (abs(a) > eps && abs(b) <= eps) {
    out <- y1 * exp(log(y2 / y1) * (1 - exp(-a * (t - t1))) /
                       (1 - exp(-a * (t2 - t1))))
  } else if (abs(a) <= eps && abs(b) > eps) {
    out <- (y1^b + (y2^b - y1^b) * (t - t1) / (t2 - t1))^(1 / b)
  } else {
    out <- y1 * exp(log(y2 / y1) * (t - t1) / (t2 - t1))
  }
  out
}

#' Schnute growth model (case 2: a != 0, b = 0)
#'
#' Special case of the general Schnute (1981) model when the shape
#' parameter \code{b} approaches zero (exponential-type growth in size),
#' keeping the rate parameter \code{a} different from zero. As with the
#' general case, \code{t1} and \code{t2} are fixed (not estimated)
#' reference ages.
#'
#' L(t) = y1 * exp( log(y2/y1) * (1 - exp(-a(t-t1))) / (1 - exp(-a(t2-t1))) )
#'
#' @param t Vector of ages.
#' @param a Parameter related to the growth rate (a != 0).
#' @param y1 Estimated size at the reference age \code{t1}.
#' @param y2 Estimated size at the reference age \code{t2}.
#' @param t1 Smaller (fixed, not estimated) reference age.
#' @param t2 Larger (fixed, not estimated) reference age.
#' @return Vector of predicted sizes.
#' @examples
#' schnute_case2(t = 1:10, a = 0.3, y1 = 15, y2 = 80, t1 = 1, t2 = 10)
#' @export
schnute_case2 <- function(t, a, y1, y2, t1, t2) {
  y1 * exp(log(y2 / y1) * (1 - exp(-a * (t - t1))) /
             (1 - exp(-a * (t2 - t1))))
}

#' Schnute growth model (case 3: a = 0, b != 0)
#'
#' Special case of the general Schnute (1981) model when the rate
#' parameter \code{a} approaches zero (linear-type growth in size^b),
#' keeping the shape parameter \code{b} different from zero. \code{t1}
#' and \code{t2} are fixed (not estimated) reference ages.
#'
#' L(t) = [y1^b + (y2^b - y1^b) * (t-t1)/(t2-t1)]^(1/b)
#'
#' @param t Vector of ages.
#' @param b Shape parameter (b != 0).
#' @param y1 Estimated size at the reference age \code{t1}.
#' @param y2 Estimated size at the reference age \code{t2}.
#' @param t1 Smaller (fixed, not estimated) reference age.
#' @param t2 Larger (fixed, not estimated) reference age.
#' @return Vector of predicted sizes.
#' @examples
#' schnute_case3(t = 1:10, b = 1, y1 = 15, y2 = 80, t1 = 1, t2 = 10)
#' @export
schnute_case3 <- function(t, b, y1, y2, t1, t2) {
  (y1^b + (y2^b - y1^b) * (t - t1) / (t2 - t1))^(1 / b)
}

#' Schnute growth model (case 4: a = 0, b = 0)
#'
#' Special case of the general Schnute (1981) model when both the rate
#' parameter \code{a} and the shape parameter \code{b} approach zero:
#' simple exponential growth between the two reference sizes. \code{t1}
#' and \code{t2} are fixed (not estimated) reference ages.
#'
#' L(t) = y1 * exp( log(y2/y1) * (t-t1)/(t2-t1) )
#'
#' @param t Vector of ages.
#' @param y1 Estimated size at the reference age \code{t1}.
#' @param y2 Estimated size at the reference age \code{t2}.
#' @param t1 Smaller (fixed, not estimated) reference age.
#' @param t2 Larger (fixed, not estimated) reference age.
#' @return Vector of predicted sizes.
#' @examples
#' schnute_case4(t = 1:10, y1 = 15, y2 = 80, t1 = 1, t2 = 10)
#' @export
schnute_case4 <- function(t, y1, y2, t1, t2) {
  y1 * exp(log(y2 / y1) * (t - t1) / (t2 - t1))
}

#' Schnute-Richards growth model (unified curve, 5 parameters)
#'
#' Unified formulation of Schnute and Richards (1990) that generalises
#' most classic sigmoid curves through the shape parameters \code{d} and
#' \code{g}. It is the most flexible model in the family included in the
#' package and usually needs good starting values.
#'
#' L(t) = Linf * (1 - b * exp(-c * t^d))^(1 / g)
#'
#' @param t Vector of ages.
#' @param Linf Asymptotic size.
#' @param b Scale parameter (typically 0 < b <= 1).
#' @param c Growth rate.
#' @param d Early-life shape parameter.
#' @param g Late-life shape parameter (g > 0).
#' @return Vector of predicted sizes.
#' @examples
#' schnute_richards(t = 1:10, Linf = 85, b = 0.9, c = 0.3, d = 1, g = 1)
#' @export
schnute_richards <- function(t, Linf, b, c, d, g) {
  base <- 1 - b * exp(-c * t^d)
  base[base < 0] <- NA_real_
  Linf * base^(1 / g)
}

#' Seasonal von Bertalanffy growth model (Somers 1988)
#'
#' Extension of the von Bertalanffy model that incorporates a periodic
#' annual component, useful when growth shows marked seasonality (common
#' in temperate-zone fish).
#'
#' @param t Vector of ages (in years, or fraction of a year).
#' @param Linf Asymptotic size.
#' @param K Average annual growth coefficient.
#' @param t0 Hypothetical age at zero size.
#' @param C Amplitude of the seasonal oscillation (0 = no seasonality,
#'   1 = maximum oscillation, typically zero growth at the coldest point).
#' @param ts Phase: point in the year (fraction, 0-1) at which the first
#'   "winter" (minimum growth rate) starts, relative to the start of the year.
#' @return Vector of predicted sizes.
#' @examples
#' vb_seasonal(t = seq(0, 5, by = 0.25), Linf = 85, K = 0.28, t0 = -0.35,
#'              C = 0.4, ts = 0.2)
#' @export
vb_seasonal <- function(t, Linf, K, t0, C, ts) {
  two_pi <- 2 * pi
  S_t  <- (C * K / two_pi) * sin(two_pi * (t - ts))
  S_t0 <- (C * K / two_pi) * sin(two_pi * (t0 - ts))
  Linf * (1 - exp(-(K * (t - t0) + S_t - S_t0)))
}

#' Gallucci & Quinn (1979) omega-reparametrisation of von Bertalanffy
#'
#' Reparametrises the von Bertalanffy model in terms of
#' \code{omega = K * Linf} (the growth rate near size zero) instead of
#' \code{Linf} directly. Algebraically identical curve to
#' \code{\link{von_bertalanffy}} (\code{Linf = omega / K}), but \code{omega}
#' and \code{K} are typically much less correlated than \code{Linf} and
#' \code{K}, which can make the fit more numerically stable - the original
#' motivation in Gallucci & Quinn (1979).
#'
#' L(t) = (omega / K) * (1 - exp(-K * (t - t0)))
#'
#' @param t Vector of ages.
#' @param omega Growth rate parameter, \code{omega = K * Linf}.
#' @param K Growth coefficient.
#' @param t0 Hypothetical age at zero size.
#' @return Vector of predicted sizes.
#' @references Gallucci, V.F. and Quinn, T.J. (1979). Reparameterizing,
#'   fitting, and testing a simple growth model. Transactions of the
#'   American Fisheries Society, 108(1), 14-25.
#' @examples
#' gallucci_quinn(t = 1:10, omega = 0.28 * 85, K = 0.28, t0 = -0.35)
#' @export
gallucci_quinn <- function(t, omega, K, t0) {
  (omega / K) * (1 - exp(-K * (t - t0)))
}

#' Francis (1988) reparametrisation of von Bertalanffy
#'
#' Reparametrises the von Bertalanffy model in terms of the expected sizes
#' \code{L1}, \code{L2}, \code{L3} at three fixed (not estimated) reference
#' ages \code{t1 < t2 < t3}, with \code{t2 = (t1 + t3) / 2}. Chosen so that
#' \code{t1}/\code{t3} span the observed age range, this parametrisation is
#' well known to have much better statistical properties (lower parameter
#' correlation, more stable convergence) than fitting \code{Linf}/\code{K}/\code{t0}
#' directly, and to be directly comparable across studies using different
#' age ranges.
#'
#' L(t) = L1 + (L3 - L1) * (1 - r^(2*(t - t1)/(t3 - t1))) / (1 - r^2),
#' with r = (L3 - L2) / (L2 - L1)
#'
#' @param t Vector of ages.
#' @param L1 Estimated expected size at the reference age \code{t1}.
#' @param L2 Estimated expected size at the reference age \code{t2}.
#' @param L3 Estimated expected size at the reference age \code{t3}.
#' @param t1 Smallest (fixed, not estimated) reference age.
#' @param t2 Middle (fixed, not estimated) reference age, \code{= (t1+t3)/2}.
#' @param t3 Largest (fixed, not estimated) reference age.
#' @return Vector of predicted sizes.
#' @references Francis, R.I.C.C. (1988). Are growth parameters estimated
#'   from tagging and age-length data comparable? Canadian Journal of
#'   Fisheries and Aquatic Sciences, 45(6), 936-942.
#' @examples
#' francis_vb(t = 1:10, L1 = 15, L2 = 75, L3 = 84, t1 = 1, t2 = 5.5, t3 = 10)
#' @export
francis_vb <- function(t, L1, L2, L3, t1, t2, t3) {
  r <- (L3 - L2) / (L2 - L1)
  L1 + (L3 - L1) * (1 - r^(2 * (t - t1) / (t3 - t1))) / (1 - r^2)
}

#' Biphasic growth model (Lester-type, practical parametrisation)
#'
#' Piecewise model combining a linear immature phase with a von
#' Bertalanffy mature phase, joined continuously at the age at maturity
#' \code{T} - the general structure proposed by Lester et al. (2004) for
#' species whose growth rate decelerates markedly at maturity, when energy
#' starts being diverted from growth to reproduction.
#'
#' L(t) = h * (t - t0)  for t <= T
#' L(t) = Linf - (Linf - L_T) * exp(-K * (t - T))  for t > T,
#' with L_T = h * (T - t0) (continuity at T)
#'
#' \strong{Note on parametrisation:} in Lester et al.'s (2004) original
#' formulation, \code{Linf} and \code{K} of the mature phase are not
#' estimated directly - they are \emph{derived} from the immature growth
#' rate \code{h} and a reproductive-investment parameter (an
#' energy-weighted gonadosomatic index), via a relationship that involves
#' an intermediate, model-specific transformed size variable. That exact
#' reduction could not be confirmed here from a primary source with
#' enough confidence to implement without risk of a mathematical error, so
#' it is intentionally not included. This function instead estimates
#' \code{h}, \code{t0}, \code{T}, \code{Linf}, and \code{K} directly,
#' enforcing only the continuity constraint at \code{T} - the same
#' immature-linear / mature-von-Bertalanffy biphasic structure, without
#' the reproductive-investment derivation of \code{Linf}/\code{K}.
#'
#' @param t Vector of ages.
#' @param h Immature-phase linear growth rate.
#' @param t0 Immature-phase hypothetical age at zero size.
#' @param T Age at maturity (phase transition).
#' @param Linf Asymptotic size of the mature phase.
#' @param K Mature-phase growth coefficient.
#' @return Vector of predicted sizes.
#' @references Lester, N.P., Shuter, B.J. and Abrams, P.A. (2004).
#'   Interpreting the von Bertalanffy model of somatic growth in fishes:
#'   the cost of reproduction. Proceedings of the Royal Society B,
#'   271(1548), 1625-1631.
#' @examples
#' biphasic_growth(t = 1:10, h = 16, t0 = -0.6, T = 2, Linf = 85, K = 0.3)
#' @export
biphasic_growth <- function(t, h, t0, T, Linf, K) {
  L_T <- h * (T - t0)
  ifelse(t <= T,
         h * (t - t0),
         Linf - (Linf - L_T) * exp(-K * (t - T)))
}

#' Persistence growth model (Tjorve 2009)
#'
#' A derivative of the power function in which the exponent itself decays
#' exponentially with age, giving a quasi-sigmoidal (but non-asymptotic)
#' shape that approaches the power function \code{a * t^b} as age
#' increases. Formula and parametrisation taken directly from Table 1 of
#' Mercier et al. (2011), who found it among the best-supported models
#' (by AIC weight) for describing gilthead seabream (\emph{Sparus aurata})
#' growth in their whole-population dataset.
#'
#' L(t) = a * t^(b * exp(-c/t))
#'
#' @param t Vector of ages (t > 0; t = 0 is defined by continuity and
#'   evaluates to \code{a}).
#' @param a Scale parameter.
#' @param b Power-function exponent approached as t -> Inf.
#' @param c Controls how quickly the exponent decays away from \code{b} at
#'   small ages.
#' @return Vector of predicted sizes.
#' @references Mercier, L., Panfili, J., Paillon, C., N'diaye, A.,
#'   Mouillot, D. and Darnaude, A.M. (2011). Otolith reading and
#'   multi-model inference for improved estimation of age and growth in
#'   the gilthead seabream Sparus aurata (L.). Estuarine, Coastal and
#'   Shelf Science, 92(4), 534-545 (Table 1). Originally proposed as a
#'   derivative of the power function by Tjorve, E. (2009). Shapes and
#'   functions of species-area curves (II): a review of new models and
#'   parameterizations. Journal of Biogeography, 36(8), 1435-1445.
#' @examples
#' persistence_growth(t = 1:10, a = 33, b = 0.37, c = -0.4)
#' @export
persistence_growth <- function(t, a, b, c) {
  a * t^(b * exp(-c / t))
}

#' Tanaka growth model (1982)
#'
#' A quasi-sigmoidal, non-asymptotic growth model. Algebraically, it is
#' equivalent to an inverse-hyperbolic-sine curve,
#' \code{L(t) = (1/sqrt(f)) * asinh(sqrt(f/a) * (t - c)) + D} (with
#' \code{D} absorbing \code{d} and the remaining constants) - unbounded
#' but decelerating growth, similar in spirit to a logarithmic curve.
#' Formula and parametrisation taken directly from Table 1 of Mercier et
#' al. (2011), who found it the best-supported model (by AIC) for
#' describing juvenile gilthead seabream (\emph{Sparus aurata}) growth.
#'
#' L(t) = (1/sqrt(f)) * ln(2f(t-c) + 2*sqrt(fa + f^2(t-c)^2)) + d
#'
#' @param t Vector of ages.
#' @param a Scale parameter (a > 0).
#' @param c Age offset.
#' @param d Vertical offset.
#' @param f Rate/curvature parameter (f > 0).
#' @return Vector of predicted sizes.
#' @references Tanaka, M. (1982). A new growth curve which expresses
#'   infinite increase. Publication of the Amakusa Marine Biology
#'   Laboratory, 6, 167-177 (as tabulated in Mercier et al. 2011, Table 1).
#' @examples
#' tanaka_growth(t = 1:10, a = 0.003, c = 1.5, d = 111, f = 0.004)
#' @export
tanaka_growth <- function(t, a, c, d, f) {
  # Domain guards (f > 0; the sqrt/log arguments non-negative/positive):
  # during optimisation, trial parameter values that fall outside the
  # valid domain are mapped to NA rather than emitting base R's sqrt()/
  # log() "NaNs produced" warnings, matching how richards()/
  # schnute_richards() already guard their own domain restrictions.
  if (!is.finite(f) || f <= 0) return(rep(NA_real_, length(t)))
  inner <- f * a + f^2 * (t - c)^2
  inner[inner < 0] <- NA_real_
  arg <- 2 * f * (t - c) + 2 * sqrt(inner)
  arg[arg <= 0] <- NA_real_
  (1 / sqrt(f)) * log(arg) + d
}

#' Linear growth model
#'
#' The simplest possible age-size relationship: a constant absolute
#' growth rate, with no deceleration and no asymptote. Useful mainly as a
#' naive baseline to compare the more realistic (decelerating/asymptotic)
#' models above against via \code{\link{aic_table}} - if none of them beat
#' this one, the data may not show enough curvature over the observed age
#' range to justify a more complex model.
#'
#' L(t) = a + b * t
#'
#' @param t Vector of ages.
#' @param a Size at age zero (intercept).
#' @param b Constant growth rate per unit of age (slope).
#' @return Vector of predicted sizes.
#' @examples
#' linear_growth(t = 1:10, a = 10, b = 6)
#' @export
linear_growth <- function(t, a, b) {
  a + b * t
}

#' Power growth model
#'
#' Classic allometric/power-law relationship between size and age. Like
#' \code{\link{linear_growth}}, it is non-asymptotic (no size plateau),
#' making it a common, simple descriptor of early-life/juvenile growth or
#' of species that keep growing (slowly) throughout life. \code{b > 1}
#' gives accelerating growth, \code{b = 1} is linear growth through the
#' origin, and \code{0 < b < 1} gives decelerating but still unbounded
#' growth.
#'
#' L(t) = a * t^b
#'
#' @param t Vector of ages. Undefined at/below \code{t = 0} for
#'   non-integer \code{b}, so values with \code{t <= 0} return \code{NA}
#'   rather than a warning or a complex number.
#' @param a Scale parameter (predicted size at age 1).
#' @param b Power-law exponent (growth curvature).
#' @return Vector of predicted sizes.
#' @examples
#' power_growth(t = 1:10, a = 8, b = 0.7)
#' @export
power_growth <- function(t, a, b) {
  out <- a * t^b
  out[t <= 0] <- NA_real_
  out
}

#' Exponential growth model
#'
#' Constant \emph{relative} (per-unit-age) growth rate \code{b}: size
#' increases (\code{b > 0}) or decreases (\code{b < 0}) by a fixed
#' proportion per unit of age. Non-asymptotic, like
#' \code{\link{linear_growth}} and \code{\link{power_growth}} - useful as
#' a simple baseline, or to describe an early/juvenile growth phase before
#' deceleration sets in (several of the asymptotic models above are, in
#' fact, locally well approximated by an exponential at small \code{t}).
#'
#' L(t) = a * exp(b * t)
#'
#' @param t Vector of ages.
#' @param a Size at age zero.
#' @param b Instantaneous relative growth rate.
#' @return Vector of predicted sizes.
#' @examples
#' exponential_growth(t = 1:10, a = 5, b = 0.25)
#' @export
exponential_growth <- function(t, a, b) {
  a * exp(b * t)
}

#' Logarithmic growth model
#'
#' Size increases in proportion to the \emph{logarithm} of age: a simple,
#' strongly decelerating but still non-asymptotic curve - a lightweight
#' alternative to the asymptotic models above when growth clearly slows
#' down but the data give no real evidence of an eventual size plateau
#' (see also \code{\link{tanaka_growth}}, which is algebraically related
#' to this same shape at large \code{t}).
#'
#' L(t) = a + b * log(t)
#'
#' @param t Vector of ages. Undefined at \code{t = 0}, so values with
#'   \code{t <= 0} return \code{NA} rather than a warning.
#' @param a Size at age 1 (where \code{log(1) = 0}).
#' @param b Growth-rate parameter (change in size per log-unit of age).
#' @return Vector of predicted sizes.
#' @examples
#' logarithmic_growth(t = 1:10, a = 10, b = 25)
#' @export
logarithmic_growth <- function(t, a, b) {
  out <- a + b * log(t)
  out[t <= 0] <- NA_real_
  out
}

#' Hyperbolic growth model (Gulland & Holt)
#'
#' A saturating (Michaelis-Menten-type) hyperbolic curve: size rises
#' quickly at first and approaches the asymptote \code{a} as age grows,
#' with \code{b} controlling how quickly (the age at which size reaches
#' half of \code{a}). Two parameters only, so it is a lightweight
#' asymptotic alternative to von Bertalanffy/Gompertz/logistic when the
#' data support fewer degrees of freedom.
#'
#' L(t) = a * t / (b + t)
#'
#' @param t Vector of ages.
#' @param a Asymptotic size (as \code{t -> Inf}).
#' @param b Half-saturation age: the age at which size equals \code{a/2}.
#'   Undefined where \code{b + t = 0}, so those values return \code{NA}
#'   rather than an infinite or division-by-zero result.
#' @return Vector of predicted sizes.
#' @references Gulland, J. A. & S. J. Holt. 1959. Estimation of growth
#'   parameters for data at unequal time intervals. Journal du Conseil
#'   CIEM 25:47-49.
#' @examples
#' hyperbolic_growth(t = 1:10, a = 40, b = 2)
#' @export
hyperbolic_growth <- function(t, a, b) {
  denom <- b + t
  out <- a * t / denom
  out[denom == 0] <- NA_real_
  out
}

#' Ricker growth model
#'
#' A dome-shaped curve (size rises to a maximum, then declines): unlike
#' the other models here, it is not monotonic in \code{t}, so it is
#' mainly useful for weight-at-age (or other) data that genuinely decline
#' at old age (e.g. post-spawning condition loss), rather than as a
#' general-purpose size-at-age descriptor.
#'
#' L(t) = a * t * exp(-b * t)
#'
#' @param t Vector of ages.
#' @param a Scale parameter.
#' @param b Rate at which the curve declines past its peak (at
#'   \code{t = 1/b}). \code{b > 0} gives the classic rise-then-fall shape.
#' @return Vector of predicted sizes.
#' @references Ricker, W. E. 1954. Stock and recruitment. Journal of the
#'   Fisheries Research Board of Canada 11:559-623.
#' @examples
#' ricker_growth(t = 1:10, a = 20, b = 0.2)
#' @export
ricker_growth <- function(t, a, b) {
  a * t * exp(-b * t)
}

#' Beverton-Holt growth model
#'
#' Another saturating hyperbolic curve (compare \code{\link{hyperbolic_growth}}),
#' in the classic Beverton-Holt parametrisation: size approaches the
#' asymptote \code{a/b} as age grows, with initial slope \code{a}.
#'
#' L(t) = a * t / (1 + b * t)
#'
#' @param t Vector of ages.
#' @param a Initial (per-unit-age) growth rate, at \code{t} near 0.
#' @param b Controls the asymptote, \code{a/b} (as \code{t -> Inf}).
#'   Undefined where \code{1 + b*t = 0}, so those values return \code{NA}.
#' @return Vector of predicted sizes.
#' @references Beverton, R. J. H. & S. J. Holt. 1959. On the dynamics of
#'   exploited fish populations. Fisheries Investigation, London.
#' @examples
#' beverton_holt_growth(t = 1:10, a = 20, b = 0.5)
#' @export
beverton_holt_growth <- function(t, a, b) {
  denom <- 1 + b * t
  out <- a * t / denom
  out[denom == 0] <- NA_real_
  out
}

#' Gamma growth model
#'
#' A flexible, dome-shaped or decelerating curve (depending on \code{g}):
#' a power-law rise (\code{t^g}) combined with an exponential decline
#' (\code{exp(-b*t)}), the same functional shape as a Gamma probability
#' density (up to a normalising constant). Like \code{\link{ricker_growth}},
#' it is not necessarily monotonic, so it is best suited to data that show
#' a genuine peak-then-decline pattern.
#'
#' L(t) = a * t^g * exp(-b * t)
#'
#' @param t Vector of ages. Undefined at/below \code{t = 0} for
#'   non-integer \code{g}, so values with \code{t <= 0} return \code{NA}.
#' @param a Scale parameter.
#' @param b Exponential decay rate.
#' @param g Power-law (shape) exponent.
#' @return Vector of predicted sizes.
#' @references Troynikov, V. S. & H. K. Gorfine. 1998. Alternative
#'   approach for establishing legal minimum lengths for abalone based on
#'   stochastic growth models for length increment data. Journal of
#'   Shellfish Research 17:827-831.
#' @examples
#' gamma_growth(t = 1:10, a = 5, b = 0.15, g = 1.5)
#' @export
gamma_growth <- function(t, a, b, g) {
  out <- a * t^g * exp(-b * t)
  out[t <= 0] <- NA_real_
  out
}

#' Weibull growth model
#'
#' A flexible asymptotic sigmoid: \code{alpha} is the asymptotic size, and
#' the extra shape parameter \code{delta} (on top of the age-offset
#' \code{gamma} and rate \code{kappa}) lets the curve's inflection sit at
#' different relative heights than the fixed-shape models above (e.g. von
#' Bertalanffy or the logistic).
#'
#' L(t) = alpha * (1 - exp(-kappa * (t - gamma)))^delta
#'
#' @param t Vector of ages.
#' @param alpha Asymptotic size (as \code{t -> Inf}).
#' @param kappa Growth-rate parameter.
#' @param gamma Age offset (theoretical age at which the base term is zero).
#' @param delta Shape exponent. Where the base term
#'   \code{1 - exp(-kappa*(t-gamma))} is negative (i.e. \code{t < gamma}
#'   when \code{kappa > 0}), raising it to a non-integer power \code{delta}
#'   is undefined, so those values return \code{NA} rather than \code{NaN}.
#' @return Vector of predicted sizes.
#' @references Seber, G. A. & C. J. Wild. 1989. Nonlinear regression.
#'   John Wiley and Sons, Inc.
#' @examples
#' weibull_growth(t = 1:10, alpha = 40, kappa = 0.3, gamma = 0, delta = 1.5)
#' @export
weibull_growth <- function(t, alpha, kappa, gamma, delta) {
  base <- 1 - exp(-kappa * (t - gamma))
  out <- alpha * base^delta
  out[base < 0] <- NA_real_
  out
}

#' Extended power growth model
#'
#' Extends \code{\link{power_growth}} (\code{L(t) = a*t^b}) with an extra
#' term \code{c/t} added to the exponent, which fades out at large
#' \code{t} (so the curve approaches simple power-law growth with exponent
#' \code{b} for old individuals) but lets the curvature at small \code{t}
#' depart from that same exponent - a one-extra-parameter refinement over
#' the plain power model, useful when a single fixed exponent does not fit
#' both young and old individuals well.
#'
#' L(t) = a * t^(b + c/t)
#'
#' @param t Vector of ages. Undefined at/below \code{t = 0}, so values
#'   with \code{t <= 0} return \code{NA}.
#' @param a Scale parameter.
#' @param b Power-law exponent approached at large \code{t}.
#' @param c Controls the departure from exponent \code{b} at small \code{t}.
#' @return Vector of predicted sizes.
#' @references Mercier, L., J. Panfili, C. Paillon, AN'diaye, D. Mouillot
#'   & A. M. Darnaudea. 2011. Otolith reading and multi-model inference
#'   for improved estimation of age and growth in the gilthead seabream
#'   (Sparus aurata, L.). Estuarine, Coastal and Shelf Science 92:534-545.
#' @examples
#' extended_power_growth(t = 1:10, a = 8, b = 0.7, c = 1)
#' @export
extended_power_growth <- function(t, a, b, c) {
  out <- a * t^(b + c / t)
  out[t <= 0] <- NA_real_
  out
}

#' Johnson growth model
#'
#' An asymmetric asymptotic curve whose inflection point sits very low
#' (close to zero size) or is not reached at all over realistic ages - in
#' practice, this makes it look similar to von Bertalanffy (also without
#' an evident inflection point) but with a different curvature at young
#' ages. \code{Linf} is the asymptotic size and \code{t0} plays the same
#' role as in von Bertalanffy (theoretical age at size zero, typically
#' negative).
#'
#' L(t) = Linf * exp(-1 / (k * (t - t0)))
#'
#' @param t Vector of ages. Undefined where \code{k*(t-t0) = 0} (i.e. at
#'   \code{t = t0}), so those values return \code{NA}.
#' @param Linf Asymptotic size (as \code{t -> Inf}).
#' @param k Growth-rate parameter.
#' @param t0 Theoretical age at which size is (in the limit) zero.
#' @return Vector of predicted sizes.
#' @references Ricker, W. E. 1975. Computation and interpretation of
#'   biological statistics of fish populations. Bulletin of the Fisheries
#'   Research Board of Canada 191. Model attributed to Johnson; as
#'   tabulated in Oribe-Perez, I. A., I. Velazquez-Abunader & G. R.
#'   Poot-Lopez. 2020. Age and multi-model growth estimation of white
#'   grunt, Haemulon plumieri. Regional Studies in Marine Science 34:101069.
#' @examples
#' johnson_growth(t = 1:10, Linf = 40, k = 0.2, t0 = -3)
#' @export
johnson_growth <- function(t, Linf, k, t0) {
  denom <- k * (t - t0)
  out <- Linf * exp(-1 / denom)
  out[denom == 0] <- NA_real_
  out
}

#' Central registry of (age-size) growth models
#'
#' A list with the specification of each classic age-size growth model:
#' the function implementing it, its parameter names, a readable label,
#' and its formula as text. Used internally by \code{\link{fit_growth}},
#' \code{\link{aic_table}}, and the plotting functions to operate
#' generically on any model without repeating model-specific code.
#'
#' @examples
#' names(GROWTH_MODELS)
#' GROWTH_MODELS$von_bertalanffy$formula
#' GROWTH_MODELS$von_bertalanffy$parameters
#' @export
GROWTH_MODELS <- list(
  von_bertalanffy = list(
    fn = von_bertalanffy,
    parameters = c("Linf", "K", "t0"),
    label = "von Bertalanffy",
    formula = "L(t) = Linf * (1 - exp(-K * (t - t0)))"
  ),
  gompertz = list(
    fn = gompertz,
    parameters = c("Linf", "K", "t0"),
    label = "Gompertz",
    formula = "L(t) = Linf * exp(-exp(-K * (t - t0)))"
  ),
  gompertz_laird = list(
    fn = gompertz_laird,
    parameters = c("L0", "A", "k"),
    label = "Gompertz-Laird",
    formula = "L(t) = L0 * exp((A/k) * (1 - exp(-k * t)))"
  ),
  logistic = list(
    fn = logistic_growth,
    parameters = c("Linf", "K", "t0"),
    label = "Logistic",
    formula = "L(t) = Linf / (1 + exp(-K * (t - t0)))"
  ),
  richards = list(
    fn = richards,
    parameters = c("Linf", "K", "t0", "b"),
    label = "Richards",
    formula = "L(t) = Linf * (1 - exp(-K * (t - t0)))^(1/b)"
  ),
  schnute_richards = list(
    fn = schnute_richards,
    parameters = c("Linf", "b", "c", "d", "g"),
    label = "Schnute-Richards",
    formula = "L(t) = Linf * (1 - b * exp(-c * t^d))^(1/g)"
  ),
  vb_seasonal = list(
    fn = vb_seasonal,
    parameters = c("Linf", "K", "t0", "C", "ts"),
    label = "Seasonal von Bertalanffy",
    formula = "L(t) = Linf * (1 - exp(-(K(t-t0) + S(t) - S(t0))))"
  ),
  schnute = list(
    fn = schnute,
    parameters = c("a", "b", "y1", "y2"),
    label = "Schnute (case 1: general)",
    formula = "L(t) = [y1^b + (y2^b - y1^b)*(1-exp(-a(t-t1)))/(1-exp(-a(t2-t1)))]^(1/b)",
    needs_t1_t2 = TRUE
  ),
  schnute_case2 = list(
    fn = schnute_case2,
    parameters = c("a", "y1", "y2"),
    label = "Schnute (case 2: a != 0, b = 0)",
    formula = "L(t) = y1 * exp(log(y2/y1) * (1-exp(-a(t-t1)))/(1-exp(-a(t2-t1))))",
    needs_t1_t2 = TRUE
  ),
  schnute_case3 = list(
    fn = schnute_case3,
    parameters = c("b", "y1", "y2"),
    label = "Schnute (case 3: a = 0, b != 0)",
    formula = "L(t) = [y1^b + (y2^b - y1^b)*(t-t1)/(t2-t1)]^(1/b)",
    needs_t1_t2 = TRUE
  ),
  schnute_case4 = list(
    fn = schnute_case4,
    parameters = c("y1", "y2"),
    label = "Schnute (case 4: a = 0, b = 0)",
    formula = "L(t) = y1 * exp(log(y2/y1) * (t-t1)/(t2-t1))",
    needs_t1_t2 = TRUE
  ),
  gallucci_quinn = list(
    fn = gallucci_quinn,
    parameters = c("omega", "K", "t0"),
    label = "Gallucci & Quinn (1979) omega-VB",
    formula = "L(t) = (omega/K) * (1 - exp(-K * (t - t0)))"
  ),
  francis_vb = list(
    fn = francis_vb,
    parameters = c("L1", "L2", "L3"),
    label = "Francis (1988) reparametrised VB",
    formula = "L(t) = L1 + (L3-L1)*(1-r^(2(t-t1)/(t3-t1)))/(1-r^2), r=(L3-L2)/(L2-L1)",
    needs_t1_t3 = TRUE
  ),
  biphasic_growth = list(
    fn = biphasic_growth,
    parameters = c("h", "t0", "T", "Linf", "K"),
    label = "Biphasic growth (Lester-type, practical)",
    formula = "L(t) = h(t-t0) if t<=T; Linf-(Linf-L_T)exp(-K(t-T)) if t>T"
  ),
  persistence = list(
    fn = persistence_growth,
    parameters = c("a", "b", "c"),
    label = "Persistence model (Tjorve 2009)",
    formula = "L(t) = a * t^(b * exp(-c/t))"
  ),
  tanaka = list(
    fn = tanaka_growth,
    parameters = c("a", "c", "d", "f"),
    label = "Tanaka (1982)",
    formula = "L(t) = (1/sqrt(f)) * ln(2f(t-c) + 2*sqrt(fa+f^2(t-c)^2)) + d"
  ),
  linear = list(
    fn = linear_growth,
    parameters = c("a", "b"),
    label = "Linear",
    formula = "L(t) = a + b * t"
  ),
  power = list(
    fn = power_growth,
    parameters = c("a", "b"),
    label = "Power",
    formula = "L(t) = a * t^b"
  ),
  exponential = list(
    fn = exponential_growth,
    parameters = c("a", "b"),
    label = "Exponential",
    formula = "L(t) = a * exp(b * t)"
  ),
  logarithmic = list(
    fn = logarithmic_growth,
    parameters = c("a", "b"),
    label = "Logarithmic",
    formula = "L(t) = a + b * log(t)"
  ),
  hyperbolic = list(
    fn = hyperbolic_growth,
    parameters = c("a", "b"),
    label = "Hyperbolic (Gulland & Holt, 1959)",
    formula = "L(t) = a * t / (b + t)"
  ),
  ricker = list(
    fn = ricker_growth,
    parameters = c("a", "b"),
    label = "Ricker (1954)",
    formula = "L(t) = a * t * exp(-b * t)"
  ),
  beverton_holt = list(
    fn = beverton_holt_growth,
    parameters = c("a", "b"),
    label = "Beverton-Holt (1959)",
    formula = "L(t) = a * t / (1 + b * t)"
  ),
  gamma = list(
    fn = gamma_growth,
    parameters = c("a", "b", "g"),
    label = "Gamma (Troynikov & Gorfine, 1998)",
    formula = "L(t) = a * t^g * exp(-b * t)"
  ),
  weibull = list(
    fn = weibull_growth,
    parameters = c("alpha", "kappa", "gamma", "delta"),
    label = "Weibull (Seber & Wild, 1989)",
    formula = "L(t) = alpha * (1 - exp(-kappa * (t - gamma)))^delta"
  ),
  extended_power = list(
    fn = extended_power_growth,
    parameters = c("a", "b", "c"),
    label = "Extended power (Mercier et al., 2011)",
    formula = "L(t) = a * t^(b + c/t)"
  ),
  johnson = list(
    fn = johnson_growth,
    parameters = c("Linf", "k", "t0"),
    label = "Johnson (Ricker, 1975)",
    formula = "L(t) = Linf * exp(-1 / (k * (t - t0)))"
  )
)

#' Names of the models shipped with the package (before any
#' \code{\link{add_growth_model}} calls)
#'
#' Captured once, right after \code{\link{GROWTH_MODELS}} is first built,
#' so \code{add_growth_model()}/\code{\link{remove_growth_model}} can tell
#' a built-in model apart from one added at runtime - regardless of how
#' many custom models have been added or removed since.
#'
#' @keywords internal
.BUILTIN_GROWTH_MODELS <- names(GROWTH_MODELS)
