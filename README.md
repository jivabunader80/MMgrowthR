# MMgrowthR

<!-- badges: start -->
[![R-CMD-check](https://github.com/jivabunader/MMgrowthR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jivabunader/MMgrowthR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package to fit, compare, and plot **individual growth models** under a
**multi-model approach**: von Bertalanffy, Gompertz, Gompertz-Laird,
logistic, Richards, the four classic Schnute (1981) parametrisation cases
(`schnute`, `schnute_case2`, `schnute_case3`, `schnute_case4`),
Schnute-Richards, a seasonal von Bertalanffy model, two alternative von
Bertalanffy reparametrisations with better statistical properties
(Gallucci & Quinn 1979 `gallucci_quinn`, Francis 1988 `francis_vb`), a
practical biphasic (immature-linear / mature-VB) growth model
(`biphasic_growth`, Lester-type), two non-asymptotic, quasi-sigmoidal
models (`persistence`, Tjorve 2009; `tanaka`, 1982) taken directly from
Mercier et al. (2011), who found them the best-supported models for
juvenile gilthead seabream growth, four classical, non-asymptotic
single-predictor baseline curves (`linear`, `power`, `exponential`,
`logarithmic`) - useful as naive baselines to compare the more realistic
models against via `aic_table()`, or to describe growth phases where an
asymptote is not appropriate - and seven further named curves from the
growth-model literature: two more two-parameter saturating/asymptotic
forms (`hyperbolic`, Gulland & Holt 1959; `beverton_holt`, Beverton &
Holt 1959), two dome-shaped forms for data that genuinely decline at old
age (`ricker`, Ricker 1954; `gamma`, Troynikov & Gorfine 1998), the
four-parameter `weibull` (Seber & Wild 1989), the three-parameter
`extended_power` (Mercier et al. 2011), and `johnson` (Ricker 1975), an
asymptotic curve without an evident inflection point (alongside von
Bertalanffy) reported as the best-supported model for white grunt,
*Haemulon plumieri* (Oribe-Perez et al. 2020).

## Installation

From GitHub (recommended - always the latest version):

```r
install.packages("remotes")  # if not already installed
remotes::install_github("jivabunader/MMgrowthR")
```

From a local source tarball:

```r
install.packages(c("ggplot2", "patchwork", "minpack.lm", "bbmle"))  # dependencies
install.packages("MMgrowthR_0.7.4.tar.gz", repos = NULL, type = "source")
```

## Citation

If MMgrowthR is useful in your work, please cite it. Once the package is
installed, run:

```r
citation("MMgrowthR")
```

which prints the entry to use (and its BibTeX form, via
`toBibtex(citation("MMgrowthR"))`). As plain text:

> Velázquez Abunader, J.I. (2026). MMgrowthR: Multi-Model Fitting of
> Growth Curves. R package version 0.7.4.
> https://github.com/jivabunader/MMgrowthR

The repository also includes a `CITATION.cff` file, so GitHub shows a
**"Cite this repository"** button (top right of the repo page) that
generates the same citation in APA or BibTeX format directly - useful for
anyone who wants to cite the package without installing it first.

## Documentation

Every exported function - and every `print`/`summary`/`coef` method for
its result objects - has full help available with `?` once the package
is loaded, e.g. `?fit_growth`, `?plot_multimodel_bayes`,
`?print.growth_fit_bayes`. Every one of those help pages includes a
runnable **Examples** section (`example("fit_growth")` runs it directly).
A complete PDF reference manual (title, usage, arguments, value, and
examples for every documented function) is also included:
`MMgrowthR-manual.pdf`.

## Practice dataset

The package bundles a ready-to-use simulated dataset, `growth_data`
(180 simulated individuals, columns `age`/`size`/`mean_size`/`group`), so
you can try every function immediately without simulating or importing
your own data first. Ages are whole-number age classes 1-11 with several
individuals per class (a stylised declining-frequency pattern, more young
individuals than old) rather than one unique continuous age per
individual - matching how real age-reading data (otoliths, scales,
growth rings) actually looks, unlike a plain `runif()` draw. See
`?growth_data` for the exact reproducible recipe:

```r
library(MMgrowthR)
data(growth_data)
head(growth_data)
fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
fit

# group ("Male"/"Female") is ready to use with the group-comparison
# functions too - see the "Comparing growth between groups" example below.
# It was assigned by an independent random split (see ?growth_data), so it
# does NOT reflect a real sex-based growth difference in this data.
fit_full <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
fit_full
```

See `?growth_data` for details (including the exact calls used to
generate it, if you want to regenerate a fresh copy or vary its parameters).

## Typical workflow

```r
library(MMgrowthR)

# 1. Data (simulated here; with real data a data.frame with age and size
#    columns is all you need)
data <- simulate_growth_data(
  model = "von_bertalanffy",
  parameters = list(Linf = 85, K = 0.28, t0 = -0.4),
  n = 180
)

# 2. Fit ONE model (least squares or maximum likelihood)
fit <- fit_growth(data, "age", "size", model = "von_bertalanffy", method = "ls")
summary(fit)

# 3. Fit ALL models (multi-model approach)
models <- fit_multimodel(data, "age", "size", method = "ls")

# 4. Compare by AIC/AICc/BIC (with Akaike weights)
comp <- aic_table(models)
print(comp)
best <- models[[best_model(comp)]]

# 5. Confidence intervals for the parameters
profile_ci(best)                      # likelihood profile (one parameter at a time)
plot_profile(best)                    # likelihood-profile PLOTS (per parameter)
boot <- bootstrap_ci(best, R = 500)   # bootstrap (case resampling)

# 5b. JOINT confidence region for a correlated pair (e.g. Linf and K are
#     typically strongly correlated in von Bertalanffy): a bivariate
#     likelihood profile referenced to chi-squared with 2 df, rather than
#     two independent 1-df intervals read together as a (misleading)
#     rectangle.
region <- profile_ci_bivariate(best, params = c("Linf", "K"))
region
plot_profile_bivariate(best, region = region)   # 2D heatmap + confidence-region contour

# 6. Elegant plots (ggplot2)
plot_growth_fit(best, boot = boot)               # curve + bootstrap CI band
plot_growth_fit(best, profile_band = TRUE)        # curve + likelihood-profile CI band
plot_growth_fit(best, boot = boot,                 # fully customised aesthetics
                point_shape = 17, point_colour = "steelblue",
                line_colour = "black", ribbon_colour = "orange", ribbon_alpha = 0.3)
plot_multimodel(models)                            # overlaid curves (default)
plot_multimodel(models, facet = TRUE)              # one panel per model, single flat colour
plot_multimodel(models, facet = TRUE, facet_colour = "firebrick")  # ...customisable
plot_bootstrap(boot, bar_colour = "steelblue", estimate_colour = "black", ci_colour = "darkorange")

# 6b. Growth rate (dL/dt) at any vector of ages - computed numerically, so
#     it works the same way for every model, built-in or custom.
growth_rate(best, t = c(1, 5, 10))
plot_growth_rate(best, show_peak = TRUE)   # marks the age of fastest growth

# 7. Results tables (exportable to CSV/HTML)
parameter_table(best, ci = boot)
summary_table(models)
export_table(summary_table(models), "summary.csv")

# 8. Model-averaged prediction (model averaging, Burnham & Anderson)
averaged_prediction(models, comp)                       # every model
averaged_prediction(models, comp, top_n = 3)             # only the top 3 by weight
averaged_prediction(models, comp, delta_max = 4)         # only models with delta <= 4

# 9. Bayesian fitting (self-contained adaptive MCMC, no external engine
#    required), credible intervals, and Bayesian model comparison (WAIC/DIC)
fit_b <- fit_growth_bayes(data, "age", "size", model = "von_bertalanffy",
                           n_iter = 8000, n_burnin = 2000, n_chains = 3)
summary(fit_b)                                # posterior summary + Rhat + acceptance rate
plot_growth_fit_bayes(fit_b)                  # curve + joint posterior-predictive credible band
plot_posterior(fit_b)                         # posterior histograms per parameter
plot_trace_bayes(fit_b)                       # MCMC trace plots (convergence diagnostic)

models_b <- fit_multimodel_bayes(data, "age", "size",
                                  models = c("von_bertalanffy", "gompertz", "richards"))
bayes_table(models_b, criterion = "waic")     # WAIC-based comparison (recommended)
bayes_table(models_b, criterion = "dic")      # DIC-based comparison
plot_multimodel_bayes(models_b)               # overlaid curves, WAIC in the legend
plot_multimodel_bayes(models_b, facet = TRUE) # one panel per model

# 10. Comparing growth between groups (sexes, populations, cohorts...):
#     group-covariate fitting + the classic Kimura (1980)/Cerrato (1990)
#     likelihood-ratio test
data$group <- sample(c("Male", "Female"), nrow(data), replace = TRUE)  # your real grouping column

fit_full   <- fit_growth_grouped(data, "age", "size", "group", free = "all")  # Omega: all free
fit_common <- fit_growth_grouped(data, "age", "size", "group", free = NULL)   # H1: one shared curve
growth_lrt(fit_full, fit_common)          # generic LRT between any two nested fits

test <- growth_lrt_groups(data, "age", "size", "group")  # fits both extremes and tests, in one call
test
plot_growth_fit_grouped(test$fit_full)    # one curve per group, coloured by group

# intermediate hypothesis: only K shared, Linf/t0 free per group
fit_k_common <- fit_growth_grouped(data, "age", "size", "group", free = c("Linf", "t0"))
growth_lrt(fit_full, fit_k_common)
aic_table(list(full = fit_full, k_common = fit_k_common, common = fit_common))  # or compare by AIC
parameter_table(fit_full)
```

## Full demonstration script

An end-to-end script (simulation -> multi-model fit -> AIC -> CI -> profile
plots -> tables -> plots) ships with the package:

```r
source(system.file("scripts/full_workflow.R", package = "MMgrowthR"))
```

## Extending the catalogue

The ~27 built-in curves are a starting point, not a ceiling. If the shape you
need is not in the catalogue, write it as a plain R function and register it
with `add_growth_model()` - it then works exactly like a built-in model with
`fit_growth()`, `fit_multimodel()`, `fit_growth_bayes()`,
`fit_multimodel_bayes()`, `simulate_growth_data()`, every plotting function,
`aic_table()`, and so on, for the rest of the current R session:

```r
add_growth_model(
  fn = function(t, a, b) a * sqrt(t) + b,     # 't' must be the first argument
  name = "sqrt_growth",                       # what you'll pass as model = "sqrt_growth"
  label = "Square-root growth",               # display label (defaults to 'name')
  formula = "L(t) = a * sqrt(t) + b",         # cosmetic only, shown in plots/tables
  colour = "#009e73"                          # used consistently across every plot
)

names(GROWTH_MODELS)   # "sqrt_growth" is now listed, right alongside the built-ins
fit_growth(data = growth_data, age = "age", size = "size", model = "sqrt_growth")

remove_growth_model(name = "sqrt_growth")   # undo it when no longer needed
```

Note that `start` (a starting-value heuristic) was **not** supplied above, and
the fit still converges. A short smoke test runs automatically inside
`add_growth_model()` (calling `fn` once with test values), so a common
mistake - wrong argument names, or a non-vectorised `if()` where `ifelse()`
is needed - is reported immediately with a clear message rather than failing
obscurely later inside the optimiser; but writing a good *starting value*
for an arbitrary model requires understanding that model's own mathematics,
which defeats the point of a quick escape hatch for a model that just isn't
in the catalogue. So when `start` is left out, `fit_growth()` runs an
**automatic multi-start search** instead of just guessing 1 for every
parameter: many candidate parameter combinations - spanning several orders
of magnitude, scaled to the data's own age/size ranges - are scored directly
against the data (no optimiser run needed for the scoring itself), and the
best-scoring combination seeds `minpack.lm`/`bbmle`. This is what makes
`start` genuinely optional rather than merely optional-in-name. It can still
be registered - as a fixed list or as a `function(t, y) list(...)` heuristic,
exactly as before - which remains useful for speed or for an unusually
difficult/multi-modal model, and always takes priority over the automatic
search when supplied. Built-in models can never be overwritten or removed
this way. Registration is **not persisted** - it lasts for the current R
session only, so `add_growth_model()` needs to be re-run after a fresh
`library(MMgrowthR)`.

## Growth rate

The growth curve `L(t)` describes size-at-age; `growth_rate()` and
`plot_growth_rate()` describe how *fast* size is changing at a given age -
the derivative `dL/dt` - computed **numerically** (finite differences of the
fitted curve, with an automatic one-sided fallback right at a model's
domain boundary, e.g. `power_growth` being undefined at `t = 0`). Because
this is numerical rather than a per-model analytical formula, it works
identically for every model in `GROWTH_MODELS` - built-in or added via
`add_growth_model()` - with no extra code needed:

```r
growth_rate(best, t = c(0, 1, 2, 5, 10))          # absolute rate (dL/dt) at chosen ages
growth_rate(best, type = "relative")               # relative rate, (1/L) dL/dt, default 200-age grid
growth_rate(best, type = "both")                   # BOTH rates in one call (two columns, one data frame)

plot_growth_rate(best)                              # the rate curve, ggplot2
plot_growth_rate(best, type = "relative")
plot_growth_rate(best, type = "both")               # BOTH curves, one figure (two stacked panels)
plot_growth_rate(best, show_peak = TRUE)            # marks the age of fastest growth
```

`type = "absolute"` (the default) gives the absolute growth rate in size
units per age unit (e.g. cm/year); `type = "relative"` gives the relative
(specific) growth rate, `(1/L(t)) * dL/dt`, in proportion per age unit -
useful for comparing growth speed across models or individuals of very
different absolute size. `type = "both"` computes/plots both at once (at
no extra cost in `growth_rate()` - the relative rate is just the absolute
rate divided by size, and the absolute rate is always computed internally
regardless of `type`): `growth_rate(..., type = "both")` returns
`rate_absolute` and `rate_relative` as two separate columns instead of one
`rate` column, and `plot_growth_rate(..., type = "both")` returns a
`patchwork` object with the two rates as separate, independently-scaled
panels (`ncol = 1` stacks them vertically, the default; `ncol = 2` places
them side by side) - since the two rates live on different scales/units,
overlaying them on one shared y-axis would be misleading, the same
reasoning `plot_profile()` uses for its own dual-scale curves. `show_peak
= TRUE` marks the age of maximum rate within the plotted range (applied
to each panel independently when `type = "both"`): most informative for a
sigmoid model with a genuine interior inflection point (Gompertz,
logistic, Richards), less so for a model like von Bertalanffy whose rate
simply decreases from the youngest age shown (there, the mark falls at
the edge of the range).

`plot_growth_rate()` also has full visual customization, mirroring
`plot_growth_fit()`'s styling options:

```r
plot_growth_rate(best, line_colour = "firebrick", linewidth = 1.5, linetype = "dashed")
plot_growth_rate(best, type = "both", colours = c("steelblue", "darkorange"))  # one colour per panel
plot_growth_rate(best, zero_line = FALSE)                                     # hide the y = 0 reference line
plot_growth_rate(best, show_peak = TRUE, peak_colour = "firebrick")           # recolour the peak marker
```

`line_colour` sets the curve colour (defaults to the model's fixed
`MODEL_PALETTE` colour); `colours` is a length-2 vector `c(colour_abs,
colour_rel)` used only with `type = "both"`, to give the two panels
different colours instead of sharing `line_colour`; `linewidth` and
`linetype` control the curve's line style; `zero_line`/`zero_line_colour`
toggle and recolour the horizontal reference line at rate = 0; and
`peak_colour` recolours the vertical marker/label added by `show_peak =
TRUE`.

## Methodological notes

- **Least-squares fitting** (`method = "ls"`): `minpack.lm::nlsLM`
  (Levenberg-Marquardt), robust to imperfect starting values.
- **Maximum-likelihood fitting** (`method = "mle"`): `bbmle::mle2`, with
  additive normal or multiplicative lognormal error (`distribution`).
- **Confidence intervals**:
  - `profile_ci()`: numeric likelihood-profile CI (via `stats::confint.nls`
    or `bbmle::confint`). When a parameter's profile is poorly behaved
    (flat and/or non-monotonic - typically a weakly identified parameter,
    e.g. `Linf` when no individuals are old enough to approach the
    asymptote), `bbmle` normally reports this as several cryptic low-level
    warnings; `profile_ci()`/`profile_trace()`/`plot_profile()` catch
    those and replace them with a single clear warning, and flag the
    affected row(s) of the result via a `reliable` column, rather than
    silently returning a possibly-bad CI.
  - `plot_growth_fit(..., profile_band = TRUE)`: draws the curve's
    confidence band from the likelihood profile rather than from
    bootstrap: for each parameter, the curve is swept across that
    parameter's profile values within the chi-squared confidence
    threshold (holding the other parameters at their point estimate), and
    the band is the pointwise envelope (min/max) of all of those curves.
    Falls back to a warning (no band drawn) whenever the profile itself
    cannot be computed - `boot = bootstrap_ci(fit)` is the recommended
    alternative then.
  - `plot_profile()`: the corresponding likelihood-profile **plots**, in
    the classic dual-axis style (chi-squared probability, dashed, left
    axis; profile log-likelihood, solid, right axis; one panel per
    parameter, combined with `patchwork`). The chi-squared curve crosses
    the horizontal `alpha = 1 - level` reference line exactly at the
    profile CI bounds, which are also drawn as vertical lines. Pass a
    `ci` you already computed (e.g. `ci <- profile_ci(fit, profile_method =
    "spline"); plot_profile(fit, ci = ci)`) to guarantee the plot shows
    exactly that result instead of `plot_profile()` silently recomputing
    its own CI internally with default arguments. `ci` is validated
    against `fit` both by parameter name and by point estimate, so a `ci`
    computed for a different model errors clearly - even when the two
    models happen to share parameter names (e.g. `von_bertalanffy` and
    `gompertz` both use `Linf`, `K`, `t0`).
  - `profile_ci_bivariate()` / `plot_profile_bivariate()`: for models
    where two parameters are strongly correlated (the classic case is
    `Linf` and `K` in von Bertalanffy: many slightly-smaller-`Linf`/
    slightly-larger-`K` combinations fit almost as well as the joint
    optimum), the two independent `profile_ci()` intervals, read together
    as a rectangle, generally do **not** match the shape of the actual
    joint region of statistically-as-good-as-the-optimum parameter pairs.
    `profile_ci_bivariate()` profiles the pair jointly over a 2D grid
    (refitting every other parameter at each grid point) and compares the
    resulting deviance to `chi-squared(df = 2)` rather than `df = 1`
    (Wilks' theorem, now with one degree of freedom per jointly-profiled
    parameter); `plot_profile_bivariate()` draws the deviance surface as a
    heatmap with the confidence-region boundary overlaid as a contour -
    typically a tilted, elongated ellipse-like shape along the ridge of
    correlation, rather than the axis-aligned rectangle a reader might
    otherwise (incorrectly) infer from the two separate intervals.
  - `bootstrap_ci()`: non-parametric bootstrap (case or residual
    resampling). Recommended fallback whenever profiling fails to
    converge (this happens for some of the more flexible models —
    Richards, Schnute-Richards, seasonal von Bertalanffy — when the extra
    shape parameters are weakly identified by the data; `profile_ci()`
    and `plot_profile()` then return `NULL` with an explanatory warning).
- **Bayesian fitting** (`fit_growth_bayes()`): a self-contained adaptive
  Metropolis-Hastings sampler (Haario et al. 2001) - no external MCMC
  engine (Stan/JAGS/NIMBLE) is required, so the package stays installable
  as a single source tarball. Runs several independent chains
  (`n_chains`, default 3) so convergence can be checked via the
  Gelman-Rubin statistic (`fit$rhat`, target ~1.00-1.05) and effective
  sample size (`fit$ess`); the proposal covariance adapts during
  `n_burnin` and is then frozen (diminishing adaptation, required for a
  valid chain). Default priors are weakly informative Normal
  distributions centred on the automatic starting values, overridable
  per-parameter via `priors`. `predict(fit_b, interval = TRUE)` and
  `plot_growth_fit_bayes()` give a genuine **joint** posterior-predictive
  credible band (the curve evaluated at every posterior draw), unlike
  `plot_growth_fit(..., profile_band = TRUE)`'s one-parameter-at-a-time
  heuristic. `plot_posterior()` shows the marginal posterior of each
  parameter and `plot_trace_bayes()` the MCMC trace plots - always check
  these (and Rhat) before trusting the results, especially for models
  already known to be weakly identified by profile-likelihood/bootstrap
  (Richards, Schnute-Richards, seasonal von Bertalanffy - the same
  difficulty shows up here as poor MCMC mixing, not a bug).
- **`plot_multimodel_bayes()`**: Bayesian analogue of `plot_multimodel()`
  (same overlay/`facet`/`legend_position` interface), showing WAIC
  instead of AIC, and flagging any model with `max(rhat) > 1.1` directly
  in its legend/facet label (`[Rhat>1.1]`) as a visual convergence warning.
- **Bayesian model comparison** (`bayes_table()`, analogous to
  `aic_table()`): compares `fit_growth_bayes()` models by **WAIC**
  (recommended - the direct Bayesian analogue of AIC, computed from the
  posterior pointwise log-likelihood; more robust to the kind of
  non-convergence flagged by Rhat) or **DIC** (older, simpler, but known
  to behave erratically - e.g. a nonsensical negative effective number of
  parameters - when the underlying chains have not converged; always
  check `max_rhat` in the table, printed with a warning when > 1.1,
  before trusting DIC in particular).
- **Model selection**: AIC, AICc (small-sample corrected), and BIC, with
  Akaike weights and delta relative to the best model (`aic_table()`). The
  reported number of parameters `k` is the growth model's own structural
  parameter count (e.g. `k = 3` for von Bertalanffy: `Linf`, `K`, `t0`) -
  it does not add an extra parameter for the residual error variance
  (sigma), following the convention most common in the fish-growth
  model-selection literature (e.g. Cerrato 1990; Katsanevakis 2006). Since
  every model here estimates exactly one such variance parameter, this
  choice never changes AIC/BIC-based rankings or Akaike weights (it only
  shifts every model's AIC/BIC by the same constant); it can make a small
  difference to AICc's non-linear small-sample correction term.
- **Comparing growth between groups** (`fit_growth_grouped()`,
  `growth_lrt()`, `growth_lrt_groups()`): the classic framework of Kimura
  (1980) and Cerrato (1990) for testing whether growth differs between
  sexes, populations, cohorts, or any other grouping. `fit_growth_grouped()`
  fits a single combined model in which each structural parameter is
  either `free` (estimated **separately per group**, e.g. `Linf.A`,
  `Linf.B`) or common (**shared** across groups) - `free = "all"`
  reproduces Kimura/Cerrato's general model (Omega, equivalent to fitting
  each group independently but as one combined log-likelihood), `free =
  NULL` reproduces their fully common model (H1, mathematically identical
  to an ordinary `fit_growth()` ignoring group), and any other subset of
  `free` gives an intermediate hypothesis (e.g. only `K` shared). Compare
  any two nested fits with `growth_lrt()` (Wilks' theorem: `LR =
  2*(logLik_full - logLik_reduced) ~ chi-sq(df = k_full - k_reduced)`), or
  run the classic full-vs-common test directly with `growth_lrt_groups()`,
  which fits both extremes for you in one call. `growth_fit_grouped`
  objects are structurally analogous to plain `growth_fit` objects (same
  `logLik`/`k`/`aic`/`aicc`/`bic` fields), so `aic_table()` and
  `parameter_table()` already work with them - no separate comparison
  machinery needed to also compare hypotheses by information criterion
  instead of/alongside a formal LRT. `plot_growth_fit_grouped()` draws one
  fitted curve per group, coloured by group, on top of the observed data.
- Starting values are computed automatically (`initial_values()`) via a
  classic linearisation of the von Bertalanffy model; they can always be
  overridden with the `start` argument.
- **Schnute (1981) model**: the four classic parametrisation cases are
  included as independent, separately fittable models, so they can be
  compared against each other (and against every other model) via
  `aic_table()`: `schnute` (case 1, general: a != 0, b != 0),
  `schnute_case2` (a != 0, b = 0), `schnute_case3` (a = 0, b != 0), and
  `schnute_case4` (a = 0, b = 0, simple exponential growth). All four fix
  the reference ages `t1`/`t2` at the minimum/maximum observed age and
  estimate only the reference sizes `y1`/`y2` (plus `a` and/or `b` as
  applicable).
- **`gallucci_quinn`** (Gallucci & Quinn 1979): reparametrises von
  Bertalanffy in terms of `omega = K * Linf` instead of `Linf`. Produces
  the *exact same curve* as `von_bertalanffy` (`Linf = omega / K`) - it
  exists because `omega` and `K` are typically much less correlated than
  `Linf` and `K`, which can make optimisation more stable, especially
  when the data do not clearly reach the asymptote.
- **`francis_vb`** (Francis 1988): reparametrises von Bertalanffy in terms
  of the expected sizes `L1`, `L2`, `L3` at three fixed reference ages
  `t1 < t2 < t3` (`t1`/`t3` set at the minimum/maximum observed age,
  `t2` at their midpoint, following the same `needs_t1_t2`-style
  convention as the Schnute cases). Also produces essentially the same
  curve as `von_bertalanffy`, with the same statistical-stability
  motivation as `gallucci_quinn`, and is the parametrisation most directly
  comparable across studies with different observed age ranges.
- **`biphasic_growth`** (Lester-type, practical parametrisation): a
  piecewise model with a linear immature phase (`L(t) = h*(t - t0)`) and a
  von Bertalanffy mature phase, joined continuously at the estimated age
  at maturity `T`. This captures the same immature/mature biphasic
  structure proposed by Lester et al. (2004) for species whose growth
  decelerates sharply at maturity, but note: in Lester et al.'s original
  formulation `Linf` and `K` of the mature phase are *derived* from the
  immature growth rate `h` and a reproductive-investment parameter (an
  energy-weighted gonadosomatic index) via a relationship involving an
  intermediate transformed size variable - that exact derivation could
  not be confirmed here from a primary source with enough confidence to
  implement without risk of a mathematical error, so `biphasic_growth`
  instead estimates `h`, `t0`, `T`, `Linf`, and `K` directly (enforcing
  only continuity at `T`).
- **`persistence`** (Tjorve 2009) and **`tanaka`** (Tanaka 1982): two
  non-asymptotic, quasi-sigmoidal models, both taken verbatim from Table 1
  of Mercier et al. (2011), who compared six growth models for gilthead
  seabream (*Sparus aurata*) and found the Tanaka model best-supported
  (by AIC) for juveniles, with the Persistence model contributing
  substantial weight to their multi-model-averaged curve. Neither has a
  horizontal asymptote (`Linf`): `persistence` behaves like the power
  function `a*t^b` at large ages, and `tanaka` is algebraically an
  inverse-hyperbolic-sine curve, `L(t) = (1/sqrt(f))*asinh(sqrt(f/a)*(t-c)) + D`
  - unbounded but decelerating growth. Useful alternatives when a species'
  growth does not clearly plateau over the observed age range (the same
  situation where VBGF is known to perform poorly).
