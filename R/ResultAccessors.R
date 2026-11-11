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
#'     group = "region", npcs = 2, clusters = 2, associationFeatures = 0,
#'     structure = list(spatial = FALSE, umap = FALSE, stability_blocks = 0))
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
#' x <- createPathwayAssay(spamtp_example$experiment, analyteType = "genes",
#'     assay = "transcriptome", database = spamtp_example$database,
#'     geneMapping = "ramp", organism = "custom", verbose = FALSE)
#' names(pathwayMapping(x, "pathway"))
pathwayMapping <- function(object, assay = "main") {
    result <- S4Vectors::metadata(.experimentForAssay(object, assay))$pathway_mapping
    if (is.null(result)) stop("No pathway mapping is stored for experiment `", assay, "`.", call. = FALSE)
    result
}

#' Extract representation metrics from a SpaMTP workflow
#'
#' @param fit Result of [runSpaMTPWorkflow()] with structure analysis enabled.
#' @return A data frame with representation, cluster count, silhouette,
#'   spatial-neighbour agreement and reference agreement metrics.
#' @export
#' @examples
#' data(spamtp_example)
#' fit <- runSpaMTPWorkflow(spamtp_example$experiment,
#'     modalities = list(main = list(type = "metabolomics", layer = "counts")),
#'     group = "region", npcs = 2, clusters = 2, associationFeatures = 0,
#'     structure = list(spatial = FALSE, umap = FALSE, stability_blocks = 0))
#' workflowMetrics(fit)
workflowMetrics <- function(fit) {
    if (!is.list(fit) || !is.list(fit$structure) ||
        !is.data.frame(fit$structure$metrics)) {
        stop("Supply a workflow result with computed structure metrics.",
            call. = FALSE)
    }
    fit$structure$metrics
}

#' Inspect the expected ions in an annotation index
#'
#' @param index An index returned by [buildMZAnnotationIndex()].
#' @return A data frame in increasing expected m/z order, with one row per
#'   indexed compound/adduct combination. Columns are `adduct`, `expected_mz`,
#'   `formula`, `metabolite_names`, `ramp_ids` and `charge`.
#' @export
#' @examples
#' data(spamtp_example)
#' index <- buildMZAnnotationIndex(spamtp_example$metabolites, adducts = "M+H")
#' annotationIndexIons(index)
annotationIndexIons <- function(index) {
    if (!inherits(index, "spamtp_mz_index")) {
        stop("index must be created by buildMZAnnotationIndex().", call. = FALSE)
    }
    compounds <- index$compounds[index$compound_idx, , drop = FALSE]
    rules <- index$rules[index$rule_idx, , drop = FALSE]
    data.frame(adduct = rules$name, expected_mz = index$expected_mz,
        formula = compounds$formula, metabolite_names = compounds$isomers_names,
        ramp_ids = compounds$ramp_ids, charge = rules$charge,
        stringsAsFactors = FALSE)
}

#' Summarize feature intensities by group
#'
#' Reports descriptive summaries on the supplied expression scale. These are
#' pixel summaries, not estimates across independent biological replicates.
#' A named list of experiments compares the same features across processed
#' datasets, such as alternative normalization panels.
#'
#' @param data A Bioconductor experiment or named list of experiments.
#' @param features Character feature IDs to summarize; NULL selects all rows.
#' @param groupBy A colData column defining groups. Missing labels are excluded.
#' @param assay Main or alternative experiment name. With `NULL` (default),
#'   explicit feature IDs are found in the main experiment first and then in
#'   paired alternative experiments, as in [plotSpatialFeature()].
#' @param assayName Expression assay name. See [experimentAccess] for the
#'   distinction between an experiment and its expression matrices.
#' @return A data frame with dataset, feature, group, n, mean, median and detection
#'   (the fraction of pixels with nonzero intensity). No P values are computed.
#' @export
#' @examples
#' data(spamtp_example)
#' summarizeFeatures(spamtp_example$experiment, groupBy = "region")
summarizeFeatures <- function(data, features = NULL, groupBy,
                              assay = NULL, assayName = "counts") {
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
        selected <- features %||% rownames(.experimentForAssay(object, assay))
        if (!is.character(selected) || !length(selected) || anyNA(selected) || anyDuplicated(selected) ||
            any(!nzchar(selected))) {
            stop("features must contain unique, non-empty feature IDs.", call. = FALSE)
        }
        if (is.null(assay)) {
            values <- .spatialFeatureValues(object, selected, assayName)
        } else {
            expression <- .assayData(object, assay, assayName)
            if (any(!selected %in% rownames(expression))) {
                stop("features must name existing rows in the selected experiment.",
                    call. = FALSE)
            }
            values <- stats::setNames(lapply(selected, function(feature) {
                as.numeric(expression[feature, ])
            }), selected)
        }
        groups <- split(seq_len(ncol(object)), as.character(metadata[[groupBy]]))
        if (!length(groups)) stop("No non-missing groups to summarize.", call. = FALSE)
        do.call(rbind, lapply(selected, function(feature) {
            v <- values[[feature]]
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
