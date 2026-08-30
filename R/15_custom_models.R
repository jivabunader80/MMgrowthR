#' Push an updated value for a top-level exported package object
#'
#' \code{\link{GROWTH_MODELS}} and \code{\link{MODEL_PALETTE}} are ordinary
#' package-level objects, so their bindings are locked once the package is
#' loaded (as for every R package) - a plain \code{<<-} or \code{assign()}
#' would fail with "cannot change value of locked binding". This helper
#' updates \strong{both} places a value can be looked up from: the
#' package's namespace (what \code{MMgrowthR:::GROWTH_MODELS} and every
#' internal function that references \code{GROWTH_MODELS} by lexical
#' scoping see - via \code{utils::assignInNamespace()}) and, if the
#' package is currently attached, the \code{package:MMgrowthR} entry on
#' the search path (what a bare \code{GROWTH_MODELS} typed at the console
#' sees - via a temporary \code{unlockBinding()}). Without the second
#' step, a bare \code{names(GROWTH_MODELS)} right after
#' \code{\link{add_growth_model}} would confusingly still show the old
#' list, even though \code{fit_growth()} would already be using the new one.
#'
#' @param object_name Character; name of the object to update (e.g.
#'   \code{"GROWTH_MODELS"}).
#' @param value New value.
#' @param pkg Package name.
#' @return \code{value}, invisibly.
#' @keywords internal
.register_package_object <- function(object_name, value, pkg = "MMgrowthR") {
  ns <- tryCatch(asNamespace(pkg), error = function(e) NULL)
  if (!is.null(ns) && exists(object_name, envir = ns, inherits = FALSE)) {
    utils::assignInNamespace(object_name, value, ns = ns)
  }
  search_name <- paste0("package:", pkg)
  if (search_name %in% search()) {
    env <- as.environment(search_name)
    if (exists(object_name, envir = env, inherits = FALSE)) {
      locked <- bindingIsLocked(object_name, env)
      if (locked) unlockBinding(object_name, env)
      assign(object_name, value, envir = env)
      if (locked) lockBinding(object_name, env)
    }
  }
  invisible(value)
}

#' Register a new growth model at runtime
#'
#' Adds a user-defined growth model to \code{\link{GROWTH_MODELS}} (and a
#' colour to \code{\link{MODEL_PALETTE}}), so it becomes usable exactly
#' like any of the package's built-in models - with \code{fit_growth()},
#' \code{fit_multimodel()}, \code{fit_growth_bayes()},
#' \code{fit_multimodel_bayes()}, \code{\link{simulate_growth_data}},
#' \code{plot_growth_fit()}, \code{plot_multimodel()},
#' \code{\link{aic_table}}, and every other function that takes a
#' \code{model}/\code{models} name - for the rest of the current R
#' session. In other words, the catalogue of ~27 built-in curves is a
#' starting point, not a ceiling: if the shape you need is not in it,
#' write it as a plain R function and register it here instead of giving
#' up on the package.
#'
#' @param fn A function of the form \code{function(t, p1, p2, ...) ...}:
#'   \code{t} (age) must be its first argument, and it must return a
#'   numeric vector of predicted sizes the same length as \code{t}
#'   (vectorised - use \code{ifelse()}, never \code{if()}, for any
#'   piecewise/domain logic; see \code{\link{power_growth}}'s source for
#'   the pattern used to return \code{NA} outside a model's valid domain
#'   instead of an error or warning). A short smoke test (calling
#'   \code{fn} once with test values) runs automatically and reports a
#'   clear error immediately if \code{fn} does not follow this contract,
#'   rather than failing obscurely later inside the optimiser.
#' @param name Character; the model's registry key (letters, digits, and
#'   underscores only, starting with a letter) - this is what you will
#'   pass afterwards as \code{model = name} (or inside \code{models =
#'   c(...)}) to \code{fit_growth()} and friends. Cannot be the name of a
#'   built-in model.
#' @param parameters Character vector naming \code{fn}'s parameters (i.e.
#'   its arguments after \code{t}), in the order you want them
#'   reported/fitted. If \code{NULL} (default), inferred automatically as
#'   every argument of \code{fn} after \code{t}.
#' @param label Display label used in plots, tables, and \code{print()}/
#'   \code{summary()} output. Defaults to \code{name}.
#' @param formula Optional character string with the model's formula as
#'   plain text (cosmetic only - stored, never evaluated), shown by
#'   \code{show_formula = TRUE} in \code{\link{plot_growth_fit}} and in
#'   \code{GROWTH_MODELS[[name]]$formula}. Defaults to a generic
#'   placeholder mentioning \code{name}.
#' @param colour Hex colour (or any R colour name) assigned to this model
#'   in \code{\link{MODEL_PALETTE}}, used consistently across every plot
#'   (e.g. \code{\link{plot_multimodel}}).
#' @param start Optional starting-value heuristic, reused automatically
#'   every time this model is fitted without an explicit \code{start}
#'   argument. Either a fixed named list (e.g.
#'   \code{list(a = 1, b = 0.5)}) or a function
#'   \code{function(t, y) list(a = ..., b = ...)} that derives sensible
#'   heuristics from the actual age/size data being fitted (see
#'   \code{\link{initial_values}}'s source for the pattern used for the
#'   built-in models, e.g. a linearising regression). \strong{This is
#'   genuinely optional} - you do not need to already understand your own
#'   model's mathematics to get a working fit. If \code{NULL} (default)
#'   and the caller does not pass \code{start} directly to
#'   \code{fit_growth()} either, an automatic multi-start search runs
#'   instead: many candidate parameter combinations (spanning several
#'   orders of magnitude, scaled to the data's own age/size ranges) are
#'   scored directly against the data, and the best-scoring one is used
#'   to seed the optimiser - see \code{\link{.auto_multistart}} for
#'   details. Registering a \code{start} heuristic here remains useful
#'   for speed, or for an unusually difficult/multi-modal model where the
#'   automatic search is not enough, and always takes priority over it
#'   when supplied.
#' @param overwrite If \code{TRUE}, allows replacing a model previously
#'   added under the same \code{name} in this session. Built-in models
#'   can never be overwritten this way - register under a different name
#'   instead.
#' @return \code{name}, invisibly. Called for its side effect of updating
#'   \code{\link{GROWTH_MODELS}} and \code{\link{MODEL_PALETTE}} for the
#'   rest of the session - this is \strong{not persisted}: re-run
#'   \code{add_growth_model()} again after a fresh \code{library(MMgrowthR)}.
#' @seealso \code{\link{remove_growth_model}} to undo this.
#' @examples
#' # A simple two-parameter curve that is not in the built-in catalogue.
#' # 'start' is left unspecified on purpose: the automatic multi-start
#' # search finds a working starting point on its own, so a user does not
#' # need to derive a heuristic like a = diff(range(y)) / sqrt(max(t)) by hand.
#' add_growth_model(
#'   fn = function(t, a, b) a * sqrt(t) + b,
#'   name = "sqrt_growth",
#'   label = "Square-root growth",
#'   formula = "L(t) = a * sqrt(t) + b",
#'   colour = "#009e73"
#' )
#' names(GROWTH_MODELS)   # "sqrt_growth" is now listed, right alongside the built-ins
#'
#' data(growth_data)
#' fit <- fit_growth(data = growth_data, age = "age", size = "size", model = "sqrt_growth")
#' fit
#'
#' remove_growth_model(name = "sqrt_growth")   # tidy up when done
#'
#' # A 'start' heuristic can still be registered - useful mainly for speed,
#' # or for a model where the automatic search struggles - and always takes
#' # priority over it when supplied:
#' add_growth_model(
#'   fn = function(t, a, b) a * sqrt(t) + b,
#'   name = "sqrt_growth2",
#'   start = function(t, y) list(a = diff(range(y)) / sqrt(max(t)), b = min(y))
#' )
#' remove_growth_model(name = "sqrt_growth2")
#' @export
add_growth_model <- function(fn, name, parameters = NULL, label = NULL,
                              formula = NULL, colour = "gray20", start = NULL,
                              overwrite = FALSE) {
  stopifnot(is.function(fn))
  stopifnot(is.character(name), length(name) == 1, !is.na(name), nzchar(name))
  if (!grepl("^[a-zA-Z][a-zA-Z0-9_]*$", name)) {
    stop("add_growth_model(): 'name' must start with a letter and contain only ",
         "letters, digits, and underscores (got '", name, "').")
  }
  if (name %in% .BUILTIN_GROWTH_MODELS) {
    stop("add_growth_model(): '", name, "' is one of MMgrowthR's built-in models and ",
         "cannot be overwritten - choose a different 'name'.")
  }
  if (name %in% names(GROWTH_MODELS) && !overwrite) {
    stop("add_growth_model(): '", name, "' was already registered with add_growth_model() ",
         "earlier in this session. Pass overwrite = TRUE to replace it, or choose a ",
         "different name.")
  }

  fn_args <- names(formals(fn))
  if (length(fn_args) == 0 || fn_args[1] != "t") {
    stop("add_growth_model(): 'fn's first argument must be named 't' (age), e.g. ",
         "function(t, a, b) a * t + b. Got: function(",
         paste(fn_args, collapse = ", "), if (length(fn_args) > 0) ", " else "", "...).")
  }
  if (is.null(parameters)) {
    parameters <- setdiff(fn_args, "t")
    if (length(parameters) == 0) {
      stop("add_growth_model(): could not infer 'parameters' from 'fn' (it takes only 't') ",
           "- pass 'parameters' explicitly.")
    }
  } else {
    stopifnot(is.character(parameters), length(parameters) > 0)
    unknown_params <- setdiff(parameters, fn_args)
    if (length(unknown_params) > 0) {
      stop("add_growth_model(): 'parameters' name(s) not found among 'fn's arguments: ",
           paste(unknown_params, collapse = ", "), ".")
    }
  }
  # Any OTHER required (no-default) argument of fn, besides t and the
  # declared parameters, would make every fit fail with a missing-argument
  # error deep inside the optimiser - catch it here instead, with a clear
  # message pointing at the fix.
  has_no_default <- vapply(fn_args, function(a) identical(formals(fn)[[a]], quote(expr = )), logical(1))
  missing_required <- setdiff(fn_args[has_no_default], c("t", parameters))
  if (length(missing_required) > 0) {
    stop("add_growth_model(): 'fn' has required argument(s) not covered by 't' or ",
         "'parameters': ", paste(missing_required, collapse = ", "), ". Add them to ",
         "'parameters', or give them a default value in 'fn'.")
  }

  if (!is.null(start)) {
    if (is.function(start)) {
      start_args <- names(formals(start))
      if (!all(c("t", "y") %in% start_args)) {
        stop("add_growth_model(): 'start', when a function, must accept 't' and 'y' ",
             "arguments, e.g. function(t, y) list(",
             paste(sprintf("%s = 1", parameters), collapse = ", "), ").")
      }
    } else if (is.list(start)) {
      missing_p <- setdiff(parameters, names(start))
      if (length(missing_p) > 0) {
        stop("add_growth_model(): 'start' is missing a value for parameter(s): ",
             paste(missing_p, collapse = ", "), ".")
      }
    } else {
      stop("add_growth_model(): 'start' must be NULL, a named list, or a function(t, y).")
    }
  }

  # Smoke test: catch an obviously broken fn (wrong argument names,
  # non-vectorised use of if() instead of ifelse(), etc.) right now, with
  # a clear message, instead of a confusing failure deep inside
  # minpack.lm::nlsLM/bbmle::mle2 much later.
  test_t <- c(1, 2, 5)
  test_start <- if (is.list(start)) {
    start[parameters]
  } else {
    stats::setNames(as.list(rep(1, length(parameters))), parameters)
  }
  test_out <- tryCatch(do.call(fn, c(list(t = test_t), test_start)), error = function(e) e)
  if (inherits(test_out, "error")) {
    stop("add_growth_model(): calling 'fn' with test values (t = c(1, 2, 5), ",
         paste(sprintf("%s = %s", names(test_start), unlist(test_start)), collapse = ", "),
         ") failed: ", conditionMessage(test_out))
  }
  if (!is.numeric(test_out) || length(test_out) != length(test_t)) {
    stop("add_growth_model(): 'fn' must return a numeric vector the same length as 't' ",
         "(is it vectorised over t? use ifelse(), not if(), for any piecewise logic) - got ",
         "an object of length ", length(test_out), " for length(t) = ", length(test_t), ".")
  }

  label <- label %||% name
  formula <- formula %||% sprintf("L(t) = <user-defined model '%s'>", name)

  spec <- list(fn = fn, parameters = parameters, label = label, formula = formula,
               start = start, custom = TRUE)

  gm <- GROWTH_MODELS
  gm[[name]] <- spec
  pal <- MODEL_PALETTE
  pal[name] <- colour

  .register_package_object("GROWTH_MODELS", gm)
  .register_package_object("MODEL_PALETTE", pal)

  message("Model '", name, "' added",
          if (overwrite && name %in% names(GROWTH_MODELS)) " (replacing the previous version)" else "",
          ". Use model = \"", name, "\" with fit_growth(), fit_multimodel(), ",
          "simulate_growth_data(), and the rest of the package.")
  invisible(name)
}

#' Remove a model previously added with add_growth_model()
#'
#' Undoes \code{\link{add_growth_model}}: removes \code{name} from
#' \code{\link{GROWTH_MODELS}} and \code{\link{MODEL_PALETTE}} for the
#' rest of the session. Built-in models cannot be removed this way.
#'
#' @param name Character; the model name to remove.
#' @return \code{TRUE}, invisibly, if a model was removed; \code{FALSE}
#'   (with a warning) if \code{name} was not a currently registered
#'   custom model.
#' @examples
#' add_growth_model(fn = function(t, a) a * t, name = "trivial_linear")
#' "trivial_linear" %in% names(GROWTH_MODELS)
#' remove_growth_model(name = "trivial_linear")
#' "trivial_linear" %in% names(GROWTH_MODELS)
#' @export
remove_growth_model <- function(name) {
  stopifnot(is.character(name), length(name) == 1, !is.na(name))
  if (name %in% .BUILTIN_GROWTH_MODELS) {
    stop("remove_growth_model(): '", name, "' is one of MMgrowthR's built-in models and ",
         "cannot be removed.")
  }
  if (!(name %in% names(GROWTH_MODELS))) {
    warning("remove_growth_model(): '", name, "' is not a currently registered model - ",
            "nothing to remove.")
    return(invisible(FALSE))
  }

  gm <- GROWTH_MODELS
  gm[[name]] <- NULL
  pal <- MODEL_PALETTE
  pal <- pal[names(pal) != name]

  .register_package_object("GROWTH_MODELS", gm)
  .register_package_object("MODEL_PALETTE", pal)

  message("Model '", name, "' removed.")
  invisible(TRUE)
}
