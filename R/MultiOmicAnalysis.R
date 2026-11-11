#' Integrate paired modalities using native PCA embeddings
#'
#' SpatialExperiment input uses paired alternative experiments and combines
#' equal-weight PCA embeddings generated with [scrapper::runPca()]. This is not
#' Seurat WNN;
#' the historical WNN workflow remains in the published-workflow branch.
#'
#' @param multiomicData A SpatialExperiment with a primary modality and paired
#'   alternative experiments in `SingleCellExperiment::altExps()`.
#' @param reductionList Reduction names for the primary MSI and alternative
#'   modalities (default = list("spm.pca", "spt.pca")).
#' @param dimsList List containing the numeric range of principle component
#'   dimension to include for each modality (default = list(1:30,1:30)).
#' @param returnIntermediate Retain per-modality PCA results in reducedDims().
#' @param verbose Show progress messages (default = FALSE).
#' @param modalities Primary experiment (`"main"`) followed by alternative
#'   experiments to integrate. Derived pathway and merged experiments are not
#'   included automatically.
#' @param ... Additional arguments passed to [scrapper::runPca()], such as
#'   `scale`.
#'
#' @return The input container with an integrated representation. For
#'   SpatialExperiment this is `reducedDim(x, "integrated")`.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' integrated <- multiOmicIntegration(spe, dimsList = list(1:2, 1:2))
#' head(SingleCellExperiment::reducedDim(integrated, "integrated"))
methods::setGeneric(
    "multiOmicIntegration",
    function(
        multiomicData,
                     reductionList = list("spm.pca", "spt.pca"),
                     dimsList = list(1:30, 1:30),
                     returnIntermediate = FALSE,
                     verbose = FALSE,
                     modalities = c("main", "transcriptome"),
                     ...
    ) {
        methods::standardGeneric("multiOmicIntegration")
    }
)

.integrationAssay <- function(experiment) {
    available <- SummarizedExperiment::assayNames(experiment)
    selected <- intersect(c("logcounts", "normcounts", "counts"), available)
    if (!length(selected)) {
        stop(
            "Each modality needs one of these assays: logcounts, normcounts, counts.",
            call. = FALSE
        )
    }
    selected[[1L]]
}

.modalityPca <- function(experiment, reductionName, dimensions, ...) {
    if (!length(dimensions) || anyNA(dimensions) || any(dimensions < 1L) ||
        anyDuplicated(dimensions)) {
        stop(
            "Each dimsList entry must contain unique positive dimensions.",
            call. = FALSE
        )
    }
    if (!methods::is(experiment, "SingleCellExperiment")) {
        experiment <- methods::as(experiment, "SingleCellExperiment")
    }
    available <- SingleCellExperiment::reducedDimNames(experiment)
    if (!reductionName %in% available) {
        maximum <- min(nrow(experiment), ncol(experiment)) - 1L
        if (maximum < 1L) {
            stop(
                "Each modality needs at least two features and pixels.",
                call. = FALSE
            )
        }
        ncomponents <- min(max(as.integer(dimensions)), maximum)
        values <- SummarizedExperiment::assay(
            experiment,
            .integrationAssay(experiment)
        )
        SingleCellExperiment::reducedDim(experiment, reductionName) <-
            .scrapperPca(values, ncomponents, ntop = 500, ...)
    }
    embedding <- SingleCellExperiment::reducedDim(experiment, reductionName)
    retained <- dimensions[dimensions <= ncol(embedding)]
    if (!length(retained)) {
        stop(
            "Requested dimensions are absent from reduction `",
            reductionName,
            "`."
        )
    }
    list(
        experiment = experiment,
        embedding = embedding[, retained, drop = FALSE]
    )
}

#' @rdname multiOmicIntegration
#' @export
methods::setMethod(
    "multiOmicIntegration",
    "SpatialExperiment",
    function(
        multiomicData,
                     reductionList = list("spm.pca", "spt.pca"),
                     dimsList = list(1:30, 1:30),
                     returnIntermediate = FALSE,
                     verbose = FALSE,
                     modalities = c("main", "transcriptome"),
                     ...
    ) {
        if (length(modalities) < 2L || !identical(modalities[[1L]], "main") ||
            anyNA(modalities) || anyDuplicated(modalities)) {
            stop("modalities must start with main and contain distinct altExp names.",
                call. = FALSE
            )
        }
        alternatives <- modalities[-1L]
        if (!all(alternatives %in% SingleCellExperiment::altExpNames(multiomicData))) {
            stop(
                "Requested alternative modality was not found. Use addTranscriptome() ",
                "or supply modalities matching altExpNames(x).",
                call. = FALSE
            )
        }
        modalityNames <- c("main", alternatives)
        experiments <- c(
            list(main = multiomicData),
            stats::setNames(
                lapply(alternatives, function(name) {
                    SingleCellExperiment::altExp(multiomicData, name)
                }),
                alternatives
            )
        )
        if (length(reductionList) != length(experiments)) {
            stop(
                "Provide one reductionList entry per modality.",
                call. = FALSE
            )
        }
        if (length(dimsList) != length(experiments)) {
            stop("Provide one dimsList entry per modality.", call. = FALSE)
        }

        results <- lapply(seq_along(experiments), function(index) {
            .modalityPca(
                experiments[[index]],
                reductionName = as.character(reductionList[[index]])[[1L]],
                dimensions = as.integer(dimsList[[index]]),
                ...
            )
        })
        embeddings <- lapply(seq_along(results), function(index) {
            embedding <- scale(results[[index]]$embedding)
            embedding[!is.finite(embedding)] <- 0
            embedding <- embedding / sqrt(ncol(embedding))
            colnames(embedding) <- paste(
                modalityNames[[index]],
                colnames(embedding) %||% seq_len(ncol(embedding)),
                sep = "_"
            )
            embedding
        })
        integrated <- do.call(cbind, embeddings)
        SingleCellExperiment::reducedDim(
            multiomicData,
            "integrated"
        ) <- integrated

        if (isTRUE(returnIntermediate)) {
            mainReductions <- SingleCellExperiment::reducedDims(results[[1L]]$experiment)
            SingleCellExperiment::reducedDims(multiomicData) <- mainReductions
            SingleCellExperiment::reducedDim(
                multiomicData,
                "integrated"
            ) <- integrated
            for (index in seq_along(alternatives)) {
                SingleCellExperiment::altExp(
                    multiomicData,
                    alternatives[[index]]
                ) <-
                    results[[index + 1L]]$experiment
            }
        }
        S4Vectors::metadata(multiomicData)$spamtp_integration <- list(
            method = "equal-weight concatenated PCA",
            modalities = modalityNames,
            reductions = unlist(reductionList),
            dimensions = dimsList
        )
        verbose_message(
            paste0(
                "Stored a ", ncol(integrated),
                "-dimensional joint embedding in reducedDim(x, `integrated`)."
            ),
            verbose = verbose
        )
        multiomicData
    }
)

#' @rdname multiOmicIntegration
#' @export
methods::setMethod(
    "multiOmicIntegration",
    "ANY",
    function(
        multiomicData,
                     reductionList = list("spm.pca", "spt.pca"),
                     dimsList = list(1:30, 1:30),
                     returnIntermediate = FALSE,
                     verbose = FALSE,
                     modalities = c("main", "transcriptome"),
                     ...
    ) {
        .requireExperiment(multiomicData, "SpatialExperiment")
    }
)


#' Combine paired modalities in a scaled alternative experiment
#'
#' Combines selected modalities as a scaled alternative experiment.
#' Accepts paired transcriptomic, proteomic or metabolomic measurements and
#' reuses existing scaled assays where available.
#'
#' @param object A SpatialExperiment containing at least two paired modalities.
#' @param assaysToMerge At least two distinct modality names: `"main"` and/or
#'   names in `SingleCellExperiment::altExpNames(SpaMTP)`. Multiple aliases for
#' the same primary experiment are not distinct modalities. See
#'   [experimentAccess].
#' @param outputExperiment Name of the output alternative experiment, not an expression
#'   matrix (default = `"merged"`).
#' @param returnOriginal TRUE adds an altExp; FALSE returns only the merged
#'   feature space.
#'
#' @return If `returnOriginal = TRUE`, the input SpatialExperiment with
#'   `altExp(SpaMTP, outputExperiment)` containing a `scaled` assay. Otherwise a new
#'   SpatialExperiment with merged rows, the same pixels and a `scaled` assay.
#' @export
#'
#' @details
#' Prefers the native `scaled` assay, then the legacy `scale.data` assay;
#' otherwise centres and scales `logcounts`, `normcounts` or `counts` across
#' pixels. Constant features become zero.
#' The merged values are stored as scaled, never relabelled as counts.
#'
#' A modality without a scaled assay is converted to a dense matrix before
#' centring and scaling, even when its input is sparse or delayed. Each dense
#' matrix needs 8 bytes per feature-pixel entry; the combined matrix and
#' temporary copies require additional memory. For large inputs, prepare
#' delayed scaled assays with [scaleSMData()] for each modality before merging,
#' or subset the features and pixels first. Downstream consumers may still
#' materialise the merged values; see [runSpatialGraphPCA()].
#'
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = rbind(a = 1:4, b = 4:1)),
#'     spatialCoords = cbind(x = 1:4, y = 0),
#'     colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4))
#' )
#' rna <- matrix(c(1, 3, 2, 4),
#'     nrow = 1,
#'     dimnames = list("gene1", colnames(x))
#' )
#' x <- addTranscriptome(x, rna)
#' x <- createMergedModalityAssay(x, c("main", "transcriptome"))
#' merged <- SingleCellExperiment::altExp(x, "merged")
#' SummarizedExperiment::assay(merged, "scaled")
createMergedModalityAssay <- function(
    object,
    assaysToMerge,
    outputExperiment = "merged",
    returnOriginal = TRUE) {
    .requireExperiment(object, "SpatialExperiment")
    modalityNames <- SingleCellExperiment::altExpNames(object)
    if (!is.character(assaysToMerge) || length(assaysToMerge) < 2L ||
        anyNA(assaysToMerge) || any(!nzchar(assaysToMerge)) ||
        anyDuplicated(assaysToMerge)) {
        stop(
            "At least two distinct modalities must be supplied.",
            call. = FALSE
        )
    }
    modalities <- lapply(assaysToMerge, function(name) {
        .experimentForAssay(object, name)
    })
    if (anyDuplicated(match(assaysToMerge, modalityNames, nomatch = 0L))) {
        stop("Supply distinct modalities, not aliases of the same experiment.",
            call. = FALSE
        )
    }
    matrices <- Map(
        function(modality, name) {
            available <- SummarizedExperiment::assayNames(modality)
            preferred <- c(
                "scaled",
                "scale.data",
                "logcounts",
                "normcounts",
                "counts"
            )
            selected <- intersect(preferred, available)
            if (!length(selected)) stop(
                "No supported expression assay in modality ",
                name,
                call. = FALSE
            )
            selected <- selected[[1L]]
            matrix <- SummarizedExperiment::assay(modality, selected)
            if (!selected %in% c("scale.data", "scaled")) {
                matrix <- t(scale(t(as.matrix(matrix))))
                matrix[!is.finite(matrix)] <- 0
            }
            rownames(matrix) <- make.unique(paste(name, rownames(matrix), sep = "::"))
            matrix
        },
        modalities,
        assaysToMerge
    )
    merged <- do.call(rbind, matrices)
    featureData <- S4Vectors::DataFrame(
        modality = rep(assaysToMerge, vapply(matrices, nrow, integer(1))),
        original_feature = unlist(lapply(
            matrices,
            function(x) sub("^[^:]+::", "", rownames(x))
        ))
    )
    rownames(featureData) <- rownames(merged)
    mergedExperiment <- SingleCellExperiment::SingleCellExperiment(
        assays = list(scaled = merged),
        rowData = featureData
    )
    if (isTRUE(returnOriginal)) {
        SingleCellExperiment::altExp(object, outputExperiment) <- mergedExperiment
        return(object)
    }
    return(SpatialExperiment::SpatialExperiment(
        assays = list(scaled = merged),
        rowData = featureData,
        colData = SummarizedExperiment::colData(object),
        spatialCoords = SpatialExperiment::spatialCoords(object),
        imgData = SpatialExperiment::imgData(object),
        metadata = S4Vectors::metadata(object)
    ))
}
