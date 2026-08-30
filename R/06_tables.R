#' @keywords internal
.std_error <- function(fit) {
  parameters <- names(fit$coefficients)
  if (fit$method == "ls") {
    tab <- summary(fit$fit_obj)$coefficients
    se <- tab[parameters, "Std. Error"]
  } else {
    tab <- bbmle::summary(fit$fit_obj)@coef
    se <- tab[parameters, "Std. Error"]
  }
  stats::setNames(as.numeric(se), parameters)
}

#' Parameter table for a growth model
#'
#' A report-ready table with the point estimate, standard error, and
#' (optionally) the confidence interval of each model parameter.
#'
#' @param fit Object of class \code{"growth_fit"} or
#'   \code{"growth_fit_grouped"} (from \code{\link{fit_growth_grouped}} -
#'   the \code{parameter} column then holds the flattened names, e.g.
#'   \code{"Linf.Male"} for free-by-group parameters).
#' @param ci Optional: result of \code{\link{profile_ci}} or
#'   \code{\link{bootstrap_ci}} for the same \code{fit}, added as
#'   \code{ci_lower}/\code{ci_upper} columns.
#' @return Data frame (class \code{"parameter_table_growth"}).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy", method = "mle")
#' parameter_table(fit)
#' parameter_table(fit, ci = profile_ci(fit))
#'
#' fit_g <- fit_growth_grouped(growth_data, "age", "size", "group", free = "all")
#' parameter_table(fit_g)
#' @export
parameter_table <- function(fit, ci = NULL) {
  stopifnot(inherits(fit, c("growth_fit", "growth_fit_grouped")))
  parameters <- names(fit$coefficients)
  se <- .std_error(fit)

  tbl <- data.frame(
    model = fit$label,
    parameter = parameters,
    estimate = as.numeric(fit$coefficients[parameters]),
    std_error = as.numeric(se[parameters]),
    stringsAsFactors = FALSE
  )

  if (!is.null(ci)) {
    ci_df <- if (inherits(ci, "bootstrap_ci_growth")) ci$table else ci
    stopifnot(all(c("parameter", "ci_lower", "ci_upper") %in% names(ci_df)))
    if (!all(parameters %in% ci_df$parameter)) {
      stop("parameter_table(): 'ci' does not cover all parameters of 'fit' (fit: ",
           paste(parameters, collapse = ", "), "; ci: ",
           paste(unique(ci_df$parameter), collapse = ", "), "). It looks like 'ci' ",
           "was computed for a different model - pass a profile_ci()/bootstrap_ci() ",
           "result computed for this exact 'fit'.")
    }
    tbl <- merge(tbl, ci_df[, c("parameter", "ci_lower", "ci_upper")],
                 by = "parameter", sort = FALSE)
    tbl <- tbl[match(parameters, tbl$parameter), ]
  }
  rownames(tbl) <- NULL
  class(tbl) <- c("parameter_table_growth", "data.frame")
  tbl
}

#' Print a parameter_table_growth object
#'
#' Prints the parameter table with its numeric columns (estimate,
#' standard error, and confidence interval bounds if present) rounded
#' for readability.
#'
#' @param x Object of class \code{"parameter_table_growth"} (from
#'   \code{\link{parameter_table}}).
#' @param digits Number of decimal places for the numeric columns
#'   (default 4).
#' @param ... Unused.
#' @return \code{x}, invisibly (called for its printed side effect).
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' print(parameter_table(fit))
#' @export
print.parameter_table_growth <- function(x, digits = 4, ...) {
  y <- x
  num_cols <- c("estimate", "std_error", "ci_lower", "ci_upper")
  for (cc in num_cols) if (cc %in% names(y)) y[[cc]] <- round(y[[cc]], digits)
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

#' Multi-model fit summary table
#'
#' Combines the information-criteria comparison table
#' (\code{\link{aic_table}}) with a summary of each model's estimated
#' parameters, in a single table ready to report or export.
#'
#' @param fit_list List of \code{"growth_fit"} objects.
#' @param criterion Criterion for the ranking (see \code{\link{aic_table}}).
#' @return Data frame (class \code{"summary_table_growth"}).
#' @examples
#' data(growth_data)
#' mm <- fit_multimodel(growth_data, "age", "size",
#'                       models = c("von_bertalanffy", "gompertz", "logistic"))
#' summary_table(mm)
#' @export
summary_table <- function(fit_list, criterion = c("aicc", "aic", "bic")) {
  criterion <- match.arg(criterion)
  comp <- aic_table(fit_list, criterion = criterion)
  coef_txt <- vapply(comp$model, function(nm) {
    fit <- fit_list[[nm]]
    paste(sprintf("%s = %.4g", names(fit$coefficients), fit$coefficients), collapse = "; ")
  }, character(1))
  comp$estimated_parameters <- coef_txt
  class(comp) <- c("summary_table_growth", "aic_table_growth", "data.frame")
  comp
}

#' Export a results table to CSV or HTML
#'
#' @param table Data frame (or \pkg{MMgrowthR} table) to export.
#' @param file Output file path; the extension (\code{.csv}, \code{.html})
#'   determines the format.
#' @param digits Number of decimal places to use in the HTML export.
#' @return Invisibly, the path of the file written.
#' @examples
#' data(growth_data)
#' fit <- fit_growth(growth_data, "age", "size", model = "von_bertalanffy")
#' out_csv <- tempfile(fileext = ".csv")
#' export_table(parameter_table(fit), out_csv)
#' read.csv(out_csv)
#' @export
export_table <- function(table, file, digits = 4) {
  ext <- tolower(tools::file_ext(file))
  df <- as.data.frame(table)
  if (ext == "csv") {
    utils::write.csv(df, file, row.names = FALSE)
  } else if (ext %in% c("html", "htm")) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("The 'knitr' package is required to export to HTML.")
    }
    html_table <- knitr::kable(df, format = "html", digits = digits)
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      html_table <- kableExtra::kable_styling(
        html_table, full_width = FALSE,
        bootstrap_options = c("striped", "hover", "condensed")
      )
    }
    writeLines(as.character(html_table), file)
  } else {
    stop("Unsupported file format: use a .csv or .html extension.")
  }
  invisible(file)
}
