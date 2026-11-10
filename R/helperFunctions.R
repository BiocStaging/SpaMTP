#' Helper function for suppressing function progress messages
#'
#' @param message_text Character string containing the message being shown
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the messsage will be suppressed (default =
#'   TRUE).
#' @noRd
#' @noRd
#' @noRd
verbose_message <- function(message_text, verbose) {
    if (verbose) {
        message(message_text)
    }
}

# Text progress bar shown only when verbose is TRUE.
.progressBar <- function(n, verbose) {
    if (!isTRUE(verbose) || n < 1L) {
        return(NULL)
    }
    utils::txtProgressBar(min = 0, max = n, style = 3)
}

.progressTick <- function(bar, value) {
    if (!is.null(bar)) utils::setTxtProgressBar(bar, value)
}

.progressDone <- function(bar) {
    if (!is.null(bar)) close(bar)
}

#' Locale-independent order of character labels
#'
#' sort(), order(), table() and factor() follow the collation locale, and
#' R CMD build and check run with LC_COLLATE=C, so the same labels can sort
#' differently on different machines. Case-insensitive radix order is the same
#' everywhere and matches English collation for plain ASCII labels.
#'
#' @param x Character vector of labels.
#' @return .labelOrder returns an ordering permutation of x; .sortLabels
#'   returns the unique labels in that order.
#' @noRd
.labelOrder <- function(x) order(tolower(x), x, method = "radix")

.sortLabels <- function(x) {
    x <- unique(x)
    x[.labelOrder(x)]
}


#' @importFrom magrittr %>% %<>%
NULL

#' Subset features and pixels through Bioconductor containers
#'
#' Uses the container's standard bracket method so alternative experiments,
#' reduced dimensions and spatial edges remain synchronized with pixels.
#'
#' Function params/args:
#' @param object A SummarizedExperiment, including SingleCellExperiment or
#'   SpatialExperiment.
#' @param subset Logical expression indicating features/variables to keep
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param idents A vector of identity classes to keep
#' @param features A vector of feature names or indices to keep
#' @param ... Reserved; additional arguments are rejected.
#'
#' @return A subset of the same Bioconductor container class.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' subsetSPM(spe, cells = colnames(spe)[1:6], features = c("ion1", "ion3"))
subsetSPM <- function(
    object = NULL,
    subset = NULL,
    cells = NULL,
    idents = NULL,
    features = NULL,
    ...
) {
    .requireExperiment(object, "SummarizedExperiment")
    if (length(list(...))) stop("Unused arguments in ...", call. = FALSE)
    keep <- rep(TRUE, ncol(object))
    if (!missing(subset)) {
        selection <- eval(
            substitute(subset),
            .cellMetadata(object),
            parent.frame()
        )
        if (!is.null(selection)) {
            if (!is.logical(selection) || length(selection) != ncol(object)) {
                stop(
                    "subset must evaluate to one logical value per pixel.",
                    call. = FALSE
                )
            }
            keep <- keep & !is.na(selection) & selection
        }
    }
    if (!is.null(idents)) {
        labels <- if (methods::is(object, "SingleCellExperiment")) {
            SingleCellExperiment::colLabels(object)
        } else {
            NULL
        }
        if (is.null(labels)) {
            stop(
                "Set colLabels(object), or use subset with a colData column.",
                call. = FALSE
            )
        }
        keep <- keep & labels %in% idents
    }
    indices <- if (is.null(cells)) {
        seq_len(ncol(object))
    } else {
        if (is.character(cells)) match(cells, colnames(object)) else cells
    }
    if (anyNA(indices) || any(indices < 1 | indices > ncol(object)) ||
        any(indices != floor(indices))) {
        stop("cells must identify existing pixels.", call. = FALSE)
    }
    rows <- if (is.null(features)) seq_len(nrow(object)) else features
    return(object[rows, indices[keep[indices]], drop = FALSE])
}


#' Check Cardinal Version
#'
#' Checks if the Cardinal Package installed is >= version 3.6.0
#'
#' @return Boolean indicating if the version of Cardinal is >= 3.6.0
#' @noRd
check_cardinal_version <- function() {
    if (numeric_version(utils::packageVersion("Cardinal")) >= numeric_version("3.6")) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


########################################################################
.requireSuggested <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("This operation requires the suggested package '", package, "'.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}
