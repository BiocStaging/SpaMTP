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
