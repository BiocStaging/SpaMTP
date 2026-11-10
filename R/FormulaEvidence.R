# Development assistance: see inst/COPYRIGHTS.

# Cache complete element counts within one query. MetaboCoreUtils can return
# partial counts after dropping unsupported isotope names; these are unavailable
# evidence, just like a formula containing an unspecified R group.
.formula_counts <- function(formula, cache = NULL) {
    key <- paste0("formula:", formula)
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
        return(get(key, envir = cache, inherits = FALSE))
    }
    parser_state <- new.env(parent = emptyenv())
    parser_state$incomplete <- FALSE
    counts <- if (is.na(formula) || !nzchar(formula)) {
        NULL
    } else {
        withCallingHandlers(
            MetaboCoreUtils::countElements(formula)[[1L]],
            warning = function(condition) {
                parser_state$incomplete <- TRUE
                text <- conditionMessage(condition)
                invalid_formula <- startsWith(text, "The given formula '") &&
                    endsWith(text, "' contains invalid symbols.")
                dropped_element <- startsWith(text,
                    "The following names are not valid and are dropped: "
                )
                if (invalid_formula || dropped_element) {
                    invokeRestart("muffleWarning")
                }
            }
        )
    }
    if (parser_state$incomplete || !length(counts) || is.null(names(counts)) ||
        anyNA(counts) || any(!is.finite(counts) | counts < 0)) {
        counts <- NULL
    }
    if (!is.null(cache)) assign(key, counts, envir = cache)
    counts
}

.formula_elements <- function(formula, cache = NULL) {
    unique(names(.formula_counts(formula, cache)))
}

.element_count <- function(formula, element, cache = NULL) {
    counts <- .formula_counts(formula, cache)
    if (is.null(counts)) return(NA_integer_)
    as.integer(sum(counts[names(counts) == element]))
}

.warn_formula_evidence <- function(formulas) {
    formulas <- unique(formulas[!is.na(formulas) & nzchar(formulas)])
    if (length(formulas)) {
        message <- paste0(
            "Formula-based evidence is unavailable for ", length(formulas),
            " distinct matched reference formula(s). Inspect formula_status ",
            "in query results or FormulaStatus in annotationCandidates(); ",
            "exact-mass matching remains available."
        )
        condition <- structure(
            list(message = message, call = NULL, formulas = formulas),
            class = c("spamtp_formula_warning", "warning", "condition")
        )
        warning(condition)
    }
    invisible(NULL)
}
