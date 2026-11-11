#' SpaMTP argument names and migration
#'
#' Public function arguments use lower camel case, such as `assayName`,
#' `minScore` and `databaseLocalDir`. `BPPARAM` follows the BiocParallel
#' convention. This convention applies to ordinary functions and S4 methods.
#'
#' @section Updating named calls:
#' The pre-release cleanup completes the argument rename. Use the current
#' spelling in named calls and in argument lists forwarded to these functions.
#' The complete function-by-function mapping is shipped as
#' `extdata/argument-renames.csv`.
#'
#' | Previous argument | Current argument |
#' |:---|:---|
#' | `SM.assay`, `SM_assay` | `smAssay` |
#' | `ST.assay`, `ST_assay` | `stAssay` |
#' | `min_score` | `minScore` |
#' | `raw.mz.column` | `rawMzColumn` |
#' | `database_local_dir` | `databaseLocalDir` |
#' | `normalisation.type` | `normalizationType` |
#' | `DE.list` | `deList` |
#' | `logFC_threshold` | `foldChangeThreshold` |
#'
#' The migration changes function arguments. Resource columns, saved result
#' fields and documented workflow configuration fields keep their schemas.
#' For example, annotation output still includes `observed_mz`, and the
#' workflow's `structure` list still uses `stability_blocks`.
#'
#' @section Fold-change units:
#' `findAllDEMs()` and `runDE()` accept `foldChangeThreshold` as a ratio of
#' at least 1: 1.2 represents a 20 percent change and 2 represents a twofold
#' change. `demsHeatmap()` accepts `logFCThreshold` on the scale of its `logFC`
#' results. For log2 results, a threshold of 1 selects twofold changes.
#'
#' @return This documentation topic describes the argument convention. The
#'   example reads a migration data frame with `function`, `previous` and
#'   `current` columns.
#' @seealso [experimentAccess], [findAllDEMs], [demsHeatmap]
#' @examples
#' migration <- read.csv(system.file("extdata", "argument-renames.csv",
#'     package = "SpaMTP"), check.names = FALSE)
#' subset(migration, `function` == "findAllDEMs")
#' @name argumentNaming
NULL

.argumentNamingCache <- new.env(parent = emptyenv())

.argumentRenameTable <- function() {
    if (is.null(.argumentNamingCache$renames)) {
        path <- system.file(
            "extdata", "argument-renames.csv",
            package = "SpaMTP"
        )
        .argumentNamingCache$renames <- if (nzchar(path) && file.exists(path)) {
            utils::read.csv(
                path,
                stringsAsFactors = FALSE, check.names = FALSE
            )
        } else {
            data.frame(
                `function` = character(),
                previous = character(),
                current = character(),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        }
    }
    .argumentNamingCache$renames
}

.suggestArgumentName <- function(arg, fn = NULL, candidates = character()) {
    renames <- .argumentRenameTable()
    if (!is.null(fn) && nrow(renames)) {
        rows <- renames[
            renames[["function"]] == fn & renames$previous == arg, ,
            drop = FALSE
        ]
        if (nrow(rows)) {
            return(unique(rows$current)[[1L]])
        }
    }
    canonical <- function(x) tolower(gsub("[._]", "", x))
    candidates <- unique(setdiff(candidates, "..."))
    if (length(candidates)) {
        matched <- candidates[canonical(candidates) == canonical(arg)]
        if (length(matched) == 1L) {
            return(matched)
        }
    }
    if (nrow(renames)) {
        rows <- renames[renames$previous == arg, , drop = FALSE]
        targets <- unique(rows$current)
        if (length(targets) == 1L &&
            (!length(candidates) || targets %in% candidates)) {
            return(targets)
        }
    }
    NULL
}

.checkDots <- function(
    ...,
    allowed = character(),
    fn = NULL,
    formalNames = character()
) {
    if (!...length()) {
        return(invisible(NULL))
    }
    label <- if (is.null(fn)) "`...`" else paste0(fn, "()")
    argNames <- ...names()
    if (is.null(argNames) || any(!nzchar(argNames))) {
        if (!length(allowed)) {
            stop("Unused unnamed argument(s) in ", label, ".", call. = FALSE)
        }
        stop(
            "Arguments in `...` passed to ", label, " must be named.",
            call. = FALSE
        )
    }
    duplicates <- unique(argNames[duplicated(argNames)])
    if (length(duplicates)) {
        stop(
            "Duplicate argument name(s) in ", label, ": ",
            paste(duplicates, collapse = ", "), ".",
            call. = FALSE
        )
    }
    invalid <- setdiff(argNames, allowed)
    if (!length(invalid)) {
        return(invisible(NULL))
    }
    candidates <- c(formalNames, allowed)
    hints <- vapply(invalid, function(arg) {
        replacement <- .suggestArgumentName(
            arg,
            fn = fn, candidates = candidates
        )
        if (is.null(replacement) || identical(replacement, arg)) {
            ""
        } else {
            paste0("Use `", replacement, "` instead of `", arg, "`.")
        }
    }, character(1))
    hints <- hints[nzchar(hints)]
    suffix <- if (length(hints)) {
        paste0(" ", paste(hints, collapse = " "))
    } else {
        ""
    }
    stop(
        "Unused argument(s) in ", label, ": ",
        paste(invalid, collapse = ", "), ".", suffix,
        call. = FALSE
    )
}
