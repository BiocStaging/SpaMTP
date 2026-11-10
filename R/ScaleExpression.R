#' Centre and scale features within a Bioconductor experiment
#'
#' Each feature (row) is centred and/or scaled to unit variance across pixels.
#' Zero-variance features become zero. The input assay is preserved and the
#' result is stored in `scaled` by default.
#'
#' Centring makes sparse data dense, so a fully materialised scaled matrix
#' needs 8 bytes per feature-pixel entry (for example, 2,000 features by
#' 100,000 pixels need about 1.6 GB). By default the result is therefore a
#' `DelayedMatrix`: row means and standard deviations are computed by block
#' processing and the arithmetic is applied lazily when values are read.
#' This avoids storing the full scaled matrix during this step; memory also
#' includes the input, block-processing buffers and row statistics.
#'
#' Downstream analyses can still allocate dense matrices.
#' [runSpatialGraphPCA()] materialises the selected expression assay and its
#' dense solve result. [createMergedModalityAssay()] also materialises each
#' modality that still needs scaling. Subset features or pixels before those
#' calls when needed, and allow for temporary copies in addition to the
#' 8 bytes per entry. Setting `delayed = FALSE` materialises scaling
#' immediately.
#'
#' @param data A SingleCellExperiment, including SpatialExperiment.
#' @param assay Primary (`"main"`) or alternative experiment name.
#' @param assayName Input expression assay, usually `"logcounts"`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param center,scale Logical; subtract row means and/or divide by row
#'   standard deviations, as in base `scale()` applied to features.
#' @param outputAssay Name of a new assay for scaled values. Must differ from
#'   the input assay and must not already exist.
#' @param delayed Logical; store a lazily evaluated `DelayedMatrix` (default)
#'   or an in-memory `matrix`.
#'
#' @return The input container with a scaled assay in the selected experiment.
#' @export
#' @examples
#' counts <- matrix(seq_len(20),
#'     nrow = 4,
#'     dimnames = list(paste0("mz-", 101:104), paste0("p", 1:5))
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce <- normalizeSMData(sce, "LogNormalize", verbose = FALSE)
#' sce <- scaleSMData(sce)
#' SummarizedExperiment::assay(sce, "scaled")
scaleSMData <- function(
    data, assay = "main", assayName = "logcounts",
    center = TRUE, scale = TRUE, outputAssay = "scaled",
    delayed = TRUE
) {
    .requireExperiment(data, "SingleCellExperiment")
    flags <- list(center, scale, delayed)
    if (!all(vapply(flags, function(flag) {
        is.logical(flag) && length(flag) == 1L &&
            !is.na(flag)
    }, logical(1)))) {
        stop(
            "center, scale and delayed must each be TRUE or FALSE.",
            call. = FALSE
        )
    }
    experiment <- .experimentForAssay(data, assay)
    if (length(outputAssay) != 1L || is.na(outputAssay) || !nzchar(outputAssay) ||
        outputAssay %in% SummarizedExperiment::assayNames(experiment)) {
        stop("outputAssay must name a new assay; existing values are not overwritten.",
            call. = FALSE
        )
    }
    values <- DelayedArray::DelayedArray(.assayData(data, assay, assayName))
    if (!nrow(values) || ncol(values) < 2L || !.allFiniteAssay(values)) {
        stop(
            "Scaling requires finite values and at least two pixels.",
            call. = FALSE
        )
    }
    centers <- if (center) DelayedArray::rowMeans(values) else rep(0, nrow(values))
    scales <- if (!scale) {
        rep(1, nrow(values))
    } else if (center) {
        sqrt(MatrixGenerics::rowVars(values))
    } else {
        # base::scale() without centring divides by the root mean square.
        sqrt(DelayedArray::rowSums(values^2) / (ncol(values) - 1L))
    }
    names(centers) <- names(scales) <- rownames(values)
    divisor <- ifelse(scales > 0, scales, 1)
    scaled <- (values - centers) / divisor
    if (!delayed) {
        scaled <- as.matrix(scaled)
    }
    SummarizedExperiment::assay(experiment, outputAssay) <- scaled
    S4Vectors::metadata(experiment)$spamtp_scaling <- list(
        input_assay = assayName, output_assay = outputAssay,
        center = if (center) centers, scale = if (scale) scales, delayed = delayed
    )
    .replaceExperiment(data, experiment, assay)
}
