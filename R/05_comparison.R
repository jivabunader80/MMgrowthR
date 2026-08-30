#' Model comparison table (AIC / AICc / BIC)
#'
#' Builds a multi-model selection table with the number of parameters,
#' log-likelihood, AIC, AICc, and BIC of each model, together with the
#' delta relative to the best model and the Akaike weights (interpretable
#' as the relative probability that each model is the best one, within
#' the set of models considered).
#'
#' @param fit_list List of \code{"growth_fit"} objects (e.g. the result of
#'   \code{\link{fit_multimodel}}).
#' @param criterion Criterion used for ranking and weights: \code{"aicc"}
#'   (default, recommended with small samples), \code{"aic"}, or \code{"bic"}.
#' @return Data frame (class \code{"aic_table_growth"}) ordered from best
#'   to worst model according to the chosen criterion.
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' aic_table(mm)
#' @export
aic_table <- function(fit_list, criterion = c("aicc", "aic", "bic")) {
  criterion <- match.arg(criterion)
  stopifnot(length(fit_list) > 0)

  rows <- lapply(names(fit_list), function(nm) {
    fit <- fit_list[[nm]]
    data.frame(model = nm, label = fit$label, k = fit$k, n = fit$n,
               logLik = fit$logLik, AIC = fit$aic, AICc = fit$aicc, BIC = fit$bic,
               stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)

  criterion_value <- switch(criterion,
    aicc = tbl$AICc,
    aic  = tbl$AIC,
    bic  = tbl$BIC
  )
  tbl$delta <- criterion_value - min(criterion_value, na.rm = TRUE)
  tbl$rel_likelihood <- exp(-0.5 * tbl$delta)
  tbl$akaike_weight <- tbl$rel_likelihood / sum(tbl$rel_likelihood, na.rm = TRUE)
  tbl <- tbl[order(tbl$delta), ]
  tbl$rank <- seq_len(nrow(tbl))
  rownames(tbl) <- NULL

  class(tbl) <- c("aic_table_growth", "data.frame")
  attr(tbl, "criterion") <- criterion
  tbl
}

#' Print an aic_table_growth object
#'
#' Prints the criterion used for ranking, then the comparison table with
#' the numeric columns rounded for readability.
#'
#' @param x Object of class \code{"aic_table_growth"} (from
#'   \code{\link{aic_table}}).
#' @param digits Number of decimal places for the numeric columns
#'   (default 3).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' print(aic_table(mm))
#' @export
print.aic_table_growth <- function(x, digits = 3, ...) {
  cat(sprintf("Multi-model comparison (criterion: %s)\n", toupper(attr(x, "criterion"))))
  y <- x
  num_cols <- c("logLik", "AIC", "AICc", "BIC", "delta", "rel_likelihood", "akaike_weight")
  for (cc in num_cols) if (cc %in% names(y)) y[[cc]] <- round(y[[cc]], digits)
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

#' Name of the best model according to a comparison table
#'
#' @param table Object returned by \code{\link{aic_table}}.
#' @return Character with the name (key) of the model with the lowest criterion.
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' comp <- aic_table(mm)
#' best_model(comp)
#' @export
best_model <- function(table) {
  stopifnot(inherits(table, "aic_table_growth"))
  table$model[table$rank == 1]
}

#' Model-averaged prediction
#'
#' Computes the growth curve averaged across several models, weighting
#' each model's prediction by its Akaike weight (Burnham & Anderson
#' multi-model inference approach). Useful when no single model clearly
#' dominates the others (weights spread across several models).
#'
#' By default every model in \code{fit_list} is included, each weighted by
#' its Akaike weight from \code{table} (renormalised to sum to 1 over the
#' models actually used). Two optional filters let you restrict the
#' averaging to a smaller, better-supported subset of models before that
#' renormalisation - the usual way multi-model averaging is done in
#' practice, rather than including every model regardless of how poorly
#' supported it is:
#' \itemize{
#'   \item \code{top_n}: keep only the \code{top_n} models with the
#'     highest Akaike weight.
#'   \item \code{delta_max}: keep only models whose delta (relative to the
#'     best model, as reported in \code{table} for whichever criterion it
#'     was built with) is at most \code{delta_max}. Common rules of thumb
#'     (Burnham & Anderson 2002): \code{delta_max = 2} ("substantial
#'     support"), \code{delta_max = 4} or \code{7} ("considerably less
#'     support" but still sometimes retained), \code{delta_max = 10}
#'     (essentially no support, usually excluded).
#' }
#' If both are supplied, models must satisfy both (their intersection).
#' If neither is supplied, every model in \code{fit_list} is used, exactly
#' as before.
#'
#' @param fit_list List of \code{"growth_fit"} objects.
#' @param table Comparison table (result of \code{\link{aic_table}}) with
#'   the same models as \code{fit_list}.
#' @param new_t Optional vector of ages at which to evaluate the averaged
#'   prediction. If \code{NULL}, a regular grid over the observed range is used.
#' @param n Number of grid points if \code{new_t} is \code{NULL}.
#' @param top_n Optional integer: restrict the averaging to the
#'   \code{top_n} models with the highest Akaike weight (see Details).
#' @param delta_max Optional numeric: restrict the averaging to models
#'   with delta \code{<= delta_max} (see Details).
#' @return Data frame with columns \code{t} and \code{y_pred_avg}, with an
#'   attribute \code{"models_used"} listing the model names that were
#'   actually included in the average (useful to check what a filter kept).
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' comp <- aic_table(mm)
#' head(averaged_prediction(mm, comp))                # every model
#' head(averaged_prediction(mm, comp, top_n = 2))      # only the top 2 by weight
#' head(averaged_prediction(mm, comp, delta_max = 4))  # only delta <= 4
#' @export
averaged_prediction <- function(fit_list, table, new_t = NULL, n = 200,
                                 top_n = NULL, delta_max = NULL) {
  stopifnot(length(fit_list) > 0)

  model_names <- names(fit_list)
  idx <- match(model_names, table$model)
  if (any(is.na(idx))) stop("The comparison table does not contain all models in 'fit_list'.")
  weights_all <- table$akaike_weight[idx]
  delta_all <- table$delta[idx]

  keep <- rep(TRUE, length(model_names))
  if (!is.null(delta_max)) {
    keep <- keep & (delta_all <= delta_max)
  }
  if (!is.null(top_n)) {
    top_n <- min(top_n, length(model_names))
    ord <- order(weights_all, decreasing = TRUE)
    keep_top <- rep(FALSE, length(model_names))
    keep_top[ord[seq_len(top_n)]] <- TRUE
    keep <- keep & keep_top
  }
  if (!any(keep)) {
    stop("No models satisfy the given 'top_n'/'delta_max' filters; try relaxing them.")
  }

  fit_list <- fit_list[keep]
  model_names <- model_names[keep]

  if (is.null(new_t)) {
    t_range <- range(unlist(lapply(fit_list, function(a) a$t)))
    new_t <- seq(t_range[1], t_range[2], length.out = n)
  }

  preds <- sapply(model_names, function(nm) {
    stats::predict(fit_list[[nm]], new_t = new_t)$y_pred
  })
  weights <- table$akaike_weight[match(model_names, table$model)]
  weights <- weights / sum(weights)

  y_avg <- if (is.matrix(preds)) as.numeric(preds %*% weights) else preds * weights
  out <- data.frame(t = new_t, y_pred_avg = y_avg)
  attr(out, "models_used") <- model_names
  out
}
