#' Extract regional markers from a SpaMTP workflow
#'
#' @param fit Result of [runSpaMTPWorkflow()].
#' @param modality Modality name supplied to the workflow, default `"main"`.
#' @return The descriptive marker result, usable by [markerRanks()] and
#'   [demsHeatmap()]. An error identifies missing or uncomputed results.
#' @export
#' @examples
#' data(spamtp_example)
#' fit <- runSpaMTPWorkflow(spamtp_example$experiment,
#'     modalities = list(main = list(type = "metabolomics", layer = "counts")),
#'     group = "region", npcs = 2, clusters = 2, verbose = FALSE)
#' head(workflowMarkers(fit)$DEMs)
workflowMarkers <- function(fit, modality = "main") {
    if (length(modality) != 1L || is.na(modality) || !is.list(fit)) {
        stop("Supply a workflow result and one modality name.", call. = FALSE)
    }
    result <- fit$region_analysis$modalities[[modality]]$markers
    if (is.null(result)) stop("No regional markers were computed for modality `", modality, "`.", call. = FALSE)
    result
}

#' Inspect pathway identity mapping and coverage
#'
#' @param object A container produced by [createPathwayAssay()] or
#'   [createPathwayObject()], or the returned alternative experiment.
#' @param assay Experiment to inspect: `"main"` or an `altExp()` name.
#' @return The mapping audit, including coverage where pathway scores were
#'   computed. An error is raised if no pathway mapping has been stored.
#' @export
#' @examples
#' data(spamtp_example)
#' x <- createPathwayAssay(spamtp_example$experiment, analyte_type = "genes",
#'     assay = "transcriptome", database = spamtp_example$database,
#'     gene_mapping = "ramp", organism = "custom", verbose = FALSE)
#' names(pathwayMapping(x, "pathway"))
pathwayMapping <- function(object, assay = "main") {
    result <- S4Vectors::metadata(.experimentForAssay(object, assay))$pathway_mapping
    if (is.null(result)) stop("No pathway mapping is stored for experiment `", assay, "`.", call. = FALSE)
    result
}

#' Summarize feature intensities by group
#'
#' Reports descriptive summaries on the supplied expression scale. These are
#' pixel summaries, not estimates across independent biological replicates.
#' A named list of experiments compares the same features across processed
#' datasets, such as alternative normalization panels.
#'
#' @param data A Bioconductor experiment or named list of experiments.
#' @param features Feature IDs to summarize; NULL selects all rows.
#' @param groupBy A colData column defining groups. Missing labels are excluded.
#' @param assay Main or alternative experiment name.
#' @param assayName Expression assay name. See [experimentAccess] for the
#'   distinction between an experiment and its expression matrices.
#' @return A data frame with dataset, feature, group, n, mean, median and detection
#'   (the fraction of pixels with nonzero intensity). No P values are computed.
#' @export
#' @examples
#' data(spamtp_example)
#' summarizeFeatures(spamtp_example$experiment, groupBy = "region")
summarizeFeatures <- function(data, features = NULL, groupBy,
                              assay = "main", assayName = "counts") {
    datasets <- if (.isSummarizedExperiment(data)) list(data = data) else data
    if (!is.list(datasets) || !length(datasets) || is.null(names(datasets)) ||
        anyNA(names(datasets)) || any(!nzchar(names(datasets))) || anyDuplicated(names(datasets))) {
        stop("Supply an experiment or a uniquely named list of experiments.", call. = FALSE)
    }
    out <- lapply(names(datasets), function(name) {
        object <- datasets[[name]]
        metadata <- .cellMetadata(object)
        if (length(groupBy) != 1L || is.na(groupBy) || !groupBy %in% names(metadata)) {
            stop("groupBy must name a colData column.", call. = FALSE)
        }
        values <- .assayData(object, assay, assayName)
        selected <- features %||% rownames(values)
        if (!length(selected) || anyNA(selected) || anyDuplicated(selected) ||
            any(!selected %in% rownames(values))) {
            stop("features must name existing, unique rows in every dataset.", call. = FALSE)
        }
        groups <- split(seq_len(ncol(object)), as.character(metadata[[groupBy]]))
        if (!length(groups)) stop("No non-missing groups to summarize.", call. = FALSE)
        do.call(rbind, lapply(selected, function(feature) {
            v <- as.numeric(values[feature, ])
            if (!.allFiniteAssay(v)) stop("Selected intensities must be finite.", call. = FALSE)
            do.call(rbind, lapply(names(groups), function(group) {
                x <- v[groups[[group]]]
                data.frame(dataset = name, feature = feature, group = group,
                    n = length(x), mean = mean(x), median = stats::median(x),
                    detection = mean(x != 0), stringsAsFactors = FALSE)
            }))
        }))
    })
    result <- do.call(rbind, out)
    rownames(result) <- NULL
    result
}
