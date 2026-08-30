## ---------------------------------------------------------------------------
## MMgrowthR - Full workflow (demonstration)
##
## Simulates age-size data, fits a multi-model growth approach, compares the
## models by AIC/AICc, computes confidence intervals for the best model's
## parameters (likelihood profile, likelihood-profile PLOTS, and bootstrap),
## generates results tables, and produces elegant ggplot2 plots.
##
## Run with:  source(system.file("scripts/full_workflow.R", package = "MMgrowthR"))
## ---------------------------------------------------------------------------

library(MMgrowthR)
library(ggplot2)

out_dir <- "MMgrowthR_results"
if (!dir.exists(out_dir)) dir.create(out_dir)

## 1. Simulate age-size data --------------------------------------------------
set.seed(123)
data <- simulate_growth_data(
  model = "von_bertalanffy",
  parameters = list(Linf = 85, K = 0.28, t0 = -0.4),
  n = 180, age_min = 0.3, age_max = 14, cv = 0.07
)
cat("First rows of the simulated data:\n")
print(head(data))

## 2. Fit all available age-size models (least squares) -----------------------
## Includes the 4 classic Schnute (1981) parametrisation cases (schnute,
## schnute_case2, schnute_case3, schnute_case4) as independent models.
models <- fit_multimodel(data, age = "age", size = "size", method = "ls")
cat("\nFitted models:", paste(names(models), collapse = ", "), "\n")

## 3. Multi-model comparison by AICc ------------------------------------------
comp_table <- aic_table(models, criterion = "aicc")
cat("\n--- Model comparison table (AICc) ---\n")
print(comp_table)
export_table(comp_table, file.path(out_dir, "aic_table.csv"))

best_name <- best_model(comp_table)
best <- models[[best_name]]
cat(sprintf("\nBest model according to AICc: %s\n", best$label))
print(best)

## 4. Confidence intervals for the best model ---------------------------------
prof <- profile_ci(best)
cat("\n--- Likelihood-profile CI ---\n")
print(prof)

set.seed(1)
boot <- bootstrap_ci(best, R = 500, type = "cases", seed = 1)
cat("\n--- Bootstrap CI (case resampling) ---\n")
print(boot)

## 5. Results tables -----------------------------------------------------------
param_table <- parameter_table(best, ci = boot)
cat("\n--- Parameter table (best model) ---\n")
print(param_table)
export_table(param_table, file.path(out_dir, "parameter_table_best_model.csv"))

results_summary <- summary_table(models)
export_table(results_summary, file.path(out_dir, "summary_table_multimodel.html"))

## 6. Plots ----------------------------------------------------------------------
g_fit <- plot_growth_fit(best, boot = boot, title = sprintf("Best model: %s", best$label))
ggsave(file.path(out_dir, "plot_best_model.png"), g_fit, width = 8, height = 5.5, dpi = 200)

## legend_position = "inside" places the legend inside the panel instead of
## the default outside/right placement, saving lateral space (handy with
## many models); it splits into 2 columns automatically when needed.
g_multi <- plot_multimodel(models, legend_position = "inside")
ggsave(file.path(out_dir, "plot_multimodel.png"), g_multi, width = 9, height = 6, dpi = 200)

## ncol controls the number of facet columns in the bootstrap histograms
g_boot <- plot_bootstrap(boot, ncol = 2)
ggsave(file.path(out_dir, "plot_bootstrap.png"), g_boot, width = 9, height = 5.5, dpi = 200)

## 6b. Likelihood-profile plots (one panel per parameter) ---------------------
g_profile <- plot_profile(best)
if (!is.null(g_profile)) {
  ggsave(file.path(out_dir, "plot_profile_best_model.png"), g_profile,
         width = 13, height = 4.6, dpi = 200)
} else {
  cat("\nNote: the likelihood profile could not be computed for this model/data;",
      "see the warning above (bootstrap_ci() is the recommended fallback).\n")
}

## 7. Model-averaged prediction (model averaging) ------------------------------
avg_pred <- averaged_prediction(models, comp_table)
g_avg <- ggplot(avg_pred, aes(x = t, y = y_pred_avg)) +
  geom_line(colour = "#2a78d6", linewidth = 1.1) +
  labs(title = "Model-averaged growth curve (weighted by Akaike weights)",
       x = "Age", y = "Predicted size") +
  theme_MMgrowthR()
ggsave(file.path(out_dir, "plot_averaged_prediction.png"), g_avg,
       width = 8, height = 5.5, dpi = 200)

cat(sprintf("\nDone. Results saved in the '%s' folder.\n", out_dir))
