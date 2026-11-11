.annotationStatisticsInput <- function(
    data, mzAssay, mzAssayName, pathwayAssay,
    pathwayAssayName, database, threshold,
    corrWeight, nWeight
) {
    if (!methods::is(data, "SingleCellExperiment")) {
        stop("Use SingleCellExperiment/SpatialExperiment; convert Seurat input explicitly.",
            call. = FALSE
        )
    }
    if (length(threshold) != 1L || !is.finite(threshold) || abs(threshold) > 1) {
        stop("corrThreshold must be between -1 and 1.", call. = FALSE)
    }
    if (length(corrWeight) != 1L || length(nWeight) != 1L ||
        any(!is.finite(c(corrWeight, nWeight))) || min(corrWeight, nWeight) < 0 ||
        corrWeight + nWeight == 0) {
        stop(
            "Weights must be finite, non-negative and not both zero.",
            call. = FALSE
        )
    }
    expression <- .nativeExpression(data, mzAssay, mzAssayName)
    pathways <- .nativeExpression(data, pathwayAssay, pathwayAssayName)
    if (ncol(expression) < 3L) {
        stop("At least three paired pixels are required.", call. = FALSE)
    }
    if (!all(c("rampId", "pathwayRampId") %in% names(database$analytehaspathway)) ||
        !all(c("rampId", "commonName") %in% names(database$source_df))) {
        stop(
            "Database requires pathway membership and RaMP/common-name columns.",
            call. = FALSE
        )
    }
    members <- unique(as.data.frame(database$analytehaspathway)[
        ,
        c("rampId", "pathwayRampId")
    ])
    members <- members[!is.na(members$rampId) & !is.na(members$pathwayRampId), ]
    sets <- lapply(
        split(as.character(members$pathwayRampId), members$rampId),
        unique
    )
    candidates <- .annotationStatisticsCandidates(
        data,
        mzAssay,
        database$source_df
    )
    list(
        expression = expression, pathways = pathways, sets = sets,
        candidates = candidates, source = database$source_df
    )
}

.annotationStatisticsCandidates <- function(object, assay, source) {
    object <- .pathwayAnnotationObject(object, assay)
    features <- .featureMetadata(object)
    current <- .storedData(object, "mz_annotation")
    compatible <- .storedData(object, "db_3")
    if (.annotation_has_current_schema(current$results) ||
        .annotation_has_current_schema(compatible)) {
        value <- .resolve_pathway_metabolite_annotations(object,
            annotationSource = "current", scoreThreshold = 0
        )
        value <- unique(value[, c("mz_name", "ramp_id")])
        return(split(value$ramp_id, value$mz_name))
    }
    column <- intersect(
        c("all_Ramp_IDs", "Ramp_IDs", "all_Isomers_IDs"),
        names(features)
    )
    if (!length(column)) {
        stop(
            "No RaMP candidates in current annotations or rowData.",
            call. = FALSE
        )
    }
    column <- column[1L]
    values <- lapply(as.character(features[[column]]), function(value) {
        ids <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
        ids <- ids[!is.na(ids) & nzchar(ids)]
        if (column == "all_Isomers_IDs") {
            if (!"sourceId" %in% names(source)) {
                stop(
                    "Legacy annotation mapping requires sourceDf$sourceId.",
                    call. = FALSE
                )
            }
            ids <- source$rampId[source$sourceId %in% ids]
        }
        sort(unique(as.character(ids[!is.na(ids)])))
    })
    stats::setNames(values, rownames(features))
}

.annotationStatisticsName <- function(id, source) {
    matched <- as.character(source$commonName[source$rampId %in% id])
    matched <- matched[!is.na(matched) & nzchar(matched)]
    if (!length(matched)) {
        return(NA_character_)
    }
    names(sort(table(matched), decreasing = TRUE))[1L]
}

.annotationZ <- function(values) {
    deviation <- stats::sd(values)
    if (!is.finite(deviation) || deviation == 0) {
        return(rep(0, length(values)))
    }
    (values - mean(values)) / deviation
}

.scoreAnnotationFeature <- function(
    feature,
    input,
    threshold,
    corrWeight,
    nWeight
) {
    ids <- sort(unique(input$candidates[[feature]]))
    if (length(ids) < 2L) {
        return(NULL)
    }
    selected <- sort(intersect(
        unique(unlist(input$sets[ids])),
        rownames(input$pathways)
    ))
    if (!length(selected)) {
        return(NULL)
    }
    target <- as.numeric(input$expression[feature, ])
    if (any(!is.finite(target))) stop("m/z expression must be finite.", call. = FALSE)
    # One target and only its candidate pathways are materialized, not a
    # feature-by-pathway all-pairs matrix or a synthetic MSI object.
    correlations <- vapply(selected, function(pathway) {
        values <- as.numeric(input$pathways[pathway, ])
        if (!.allFiniteAssay(values)) {
            stop("Pathway expression must be finite.", call. = FALSE)
        }
        if (stats::sd(target) == 0 || stats::sd(values) == 0) {
            return(NA_real_)
        }
        stats::cor(target, values, method = "pearson")
    }, numeric(1))
    maximum <- vapply(ids, function(id) {
        values <- correlations[intersect(input$sets[[id]], selected)]
        values <- values[is.finite(values)]
        if (!length(values)) {
            return(NA_real_)
        }
        values[which.max(abs(values))]
    }, numeric(1))
    number <- vapply(ids, function(id) {
        values <- correlations[intersect(input$sets[[id]], selected)]
        as.integer(sum(values > threshold, na.rm = TRUE))
    }, integer(1))
    valid <- is.finite(maximum)
    scores <- pvalue <- adjusted <- rep(NA_real_, length(ids))
    if (sum(valid) > 1L) {
        scores[valid] <- corrWeight * .annotationZ(abs(maximum[valid])) +
            nWeight * .annotationZ(number[valid])
        pvalue[valid] <- stats::pnorm(scores[valid], lower.tail = FALSE)
        adjusted[valid] <- stats::p.adjust(pvalue[valid], method = "BH")
    }
    result <- data.frame(
        metabolite = unname(vapply(ids, .annotationStatisticsName, character(1),
            source = input$source
        )),
        ramp_id = ids, n_sig_path = unname(number), max_cor = unname(maximum),
        z_score = scores, pval = pvalue, pval_adj = adjusted
    )
    result[order(-result$z_score, result$ramp_id, na.last = TRUE), ]
}

#' Rank candidate metabolite annotations using paired pathway expression
#'
#' Pearson correlations are calculated directly between a selected m/z assay
#' and pathway expression in altExp, aligned by pixel IDs. This replaces the
#' former temporary Seurat assay and synthetic-mass Cardinal conversion.
#'
#' @details For each candidate, max_cor is the signed correlation of the
#' pathway with largest absolute correlation; ties use pathway ID order.
#' n_sig_path counts finite correlations strictly greater than `corrThreshold`
#' (zero therefore counts positive correlations, not every pathway).
#' The score is corrWeight * scale(abs(max_cor)) +
#' nWeight * scale(n_sig_path), across candidates for one feature.
#' A constant component contributes zero. Constant expression and unmeasured
#' pathways have undefined correlations. Scores require two evaluable
#' candidates; otherwise scores and probabilities are NA.
#'
#' pval is the standard-normal upper tail of this weighted score; pval_adj
#' applies BH within each feature. These are legacy heuristic ranking measures,
#' not calibrated identification p-values: the weighted components are not
#' independent, the combined score need not be standard normal, and pixels are
#' not biological replicates. Correlation is across all supplied pixels;
#' subset samples/conditions beforehand when appropriate.
#'
#' Current scored RaMP annotations take precedence over rowData RaMP IDs.
#' Older all_Isomers_IDs columns can be mapped via source_df$sourceId.
#' Duplicate candidate IDs are counted once.
#'
#' @param mz One numeric m/z query or exact character feature ID.
#' @param mzPpm Maximum ppm distance for a numeric m/z query (default 5).
#'   Equally close features require an exact character feature ID.
#' @param data A SingleCellExperiment or SpatialExperiment with paired altExp.
#' @param mzAssay Primary experiment ("main") or alternative MSI experiment.
#' @param pathwayAssay Alternative experiment containing pathway-level scores.
#' @param mzAssayName Expression assay in the MSI experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param pathwayAssayName Expression assay in the pathway experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param corrThreshold Minimum signed pathway correlation.
#' @param corrWeight,nWeight Non-negative weights, not both zero.
#' @param database Optional named list from loadSpaMTPDatabase.
#' @param databaseVersion,databaseSource,databaseLocalDir Database
#'   selection;
#'   see loadSpaMTPDatabase.
#' @return A data frame ordered by decreasing score, containing metabolite,
#'   ramp_id, n_sig_path, max_cor, z_score, pval and pval_adj.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' spe <- createPathwayObject(spe, assay = "transcriptome", pathwayIndex = index)
#' calculateSingleAnnotationStatistics(149.04445, spe,
#'     database = spamtp_example$database
#' )
calculateSingleAnnotationStatistics <- function(
    mz, data, mzAssay = "main", pathwayAssay = "pathway",
    mzAssayName = "counts", pathwayAssayName = "pathwayScores", corrThreshold = 0,
    corrWeight = 1, nWeight = 1, database = NULL,
    databaseVersion = "latest", databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL, mzPpm = 5
) {
    resources <- .spamtp_db_bundle(c("source_df", "analytehaspathway"),
        database = database, version = databaseVersion,
        source = match.arg(databaseSource), localDir = databaseLocalDir
    )
    input <- .annotationStatisticsInput(
        data, mzAssay, mzAssayName, pathwayAssay,
        pathwayAssayName, resources, corrThreshold, corrWeight, nWeight
    )
    if (length(mz) != 1L || is.na(mz)) stop("Supply one m/z query.", call. = FALSE)
    if (is.numeric(mz)) {
        if (!is.finite(mz) || mz <= 0) stop("m/z must be positive and finite.", call. = FALSE)
        if (!is.numeric(mzPpm) || length(mzPpm) != 1L || !is.finite(mzPpm) || mzPpm < 0) {
            stop("mzPpm must be one finite non-negative number.", call. = FALSE)
        }
        masses <- .massValues(
            .featureMetadata(data, mzAssay),
            rownames(input$expression)
        )
        distance <- abs(masses - mz) / mz * 1e6
        nearest <- which(is.finite(distance) & distance <= mzPpm)
        if (!length(nearest)) {
            stop("No measured feature matches m/z ", mz, " within ", mzPpm, " ppm.", call. = FALSE)
        }
        nearest <- nearest[distance[nearest] == min(distance[nearest])]
        if (length(nearest) != 1L) {
            stop("Several features are equally close to m/z ", mz,
                "; supply an exact character feature ID.", call. = FALSE)
        }
        mz <- rownames(input$expression)[nearest]
    }
    if (!mz %in% rownames(input$expression)) stop("Unknown m/z feature.", call. = FALSE)
    result <- .scoreAnnotationFeature(
        mz,
        input,
        corrThreshold,
        corrWeight,
        nWeight
    )
    if (is.null(result)) {
        stop(
            "The feature needs at least two candidate IDs and measured pathways.",
            call. = FALSE
        )
    }
    result
}

#' Rank annotations for all MSI features
#'
#' When `pathwayAssay` contains RaMP analytes rather than pathway scores,
#' createPathwayObject is reused to calculate scores in a temporary altExp.
#' Set pathwayScores=TRUE to use an existing pathway-level assay directly.
#' Input containers are never modified.
#' @details Uses the score definition and missing-data rules documented in
#' [calculateSingleAnnotationStatistics()]. The reported pval and pval_adj
#' columns are heuristic ranking measures, not calibrated identification
#' p-values or biological-replicate inference.
#' @inheritParams calculateSingleAnnotationStatistics
#' @param pathwayAssay Alternative experiment of RaMP analytes, or pathway
#'   scores when pathwayScores=TRUE.
#' @param pathwayAssayName Input expression assay in pathway.assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param returnTop Return one best evaluable candidate per feature; FALSE
#'   returns a named list of all candidate tables. Features without a score
#' remain as NA rows in top results and NULL entries when no table is possible.
#' @param pathwayScores Whether `pathwayAssay` already contains pathway scores.
#' @return A data.frame preserving feature order and rowData, with ranked
#'   annotation columns; or a named list when returnTop=FALSE.
#' @inheritParams createPathwayAssay
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' spe <- createPathwayObject(spe, assay = "transcriptome", pathwayIndex = index)
#' result <- calculateAnnotationStatistics(spe,
#'     pathwayAssayName = "pathwayScores",
#'     pathwayScores = TRUE, pathwayIndex = index, database = spamtp_example$database
#' )
#' head(result)
calculateAnnotationStatistics <- function(
    data, mzAssay = "main", pathwayAssay = "pathway",
    mzAssayName = "counts", pathwayAssayName = "counts", returnTop = TRUE,
    corrThreshold = 0, corrWeight = 1, nWeight = 1,
    database = NULL, databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"), databaseLocalDir = NULL,
    pathwayScores = FALSE, pathwayIndex = NULL,
    geneMapping = c("auto", "hgnc", "ramp"), geneReference = NULL, geneIndex = NULL,
    geneReferenceVersion = "latest", geneReferenceLocalDir = NULL,
    organism = "Homo sapiens", duplicateGenes = c("error", "mean", "sum")
) {
    source <- match.arg(databaseSource)
    database <- .pathway_database(database, pathwayIndex)
    required <- c(
        "source_df",
        "analytehaspathway",
        if (!pathwayScores) "pathway"
    )
    resources <- .spamtp_db_bundle(required,
        database = database,
        version = databaseVersion, source = source, localDir = databaseLocalDir
    )
    if (!methods::is(data, "SingleCellExperiment")) {
        stop("Use SingleCellExperiment/SpatialExperiment; convert Seurat input explicitly.",
            call. = FALSE
        )
    }
    if (!pathwayScores) {
        .nativeExpression(data, pathwayAssay, pathwayAssayName)
        newName <- tail(make.unique(c(
            SingleCellExperiment::altExpNames(data),
            ".annotationPathways"
        )), 1L)
        data <- createPathwayObject(data,
            assay = pathwayAssay, assayName = pathwayAssayName,
            outputExperiment = newName, database = resources, databaseVersion = databaseVersion,
            databaseSource = source, databaseLocalDir = databaseLocalDir,
            pathwayIndex = pathwayIndex, geneMapping = match.arg(geneMapping),
            geneReference = geneReference %||% database$gene_reference, geneIndex = geneIndex,
            geneReferenceVersion = geneReferenceVersion, geneReferenceLocalDir = geneReferenceLocalDir,
            organism = organism, duplicateGenes = match.arg(duplicateGenes)
        )
        pathwayAssay <- newName
        pathwayAssayName <- "pathwayScores"
    }
    input <- .annotationStatisticsInput(
        data, mzAssay, mzAssayName, pathwayAssay,
        pathwayAssayName, resources, corrThreshold, corrWeight, nWeight
    )
    ids <- rownames(input$expression)
    result <- stats::setNames(lapply(ids, .scoreAnnotationFeature,
        input = input,
        threshold = corrThreshold, corrWeight = corrWeight, nWeight = nWeight
    ), ids)
    audit <- S4Vectors::metadata(.experimentForAssay(data, pathwayAssay))$pathway_mapping
    if (!is.null(pathwayIndex) && !is.null(audit) && !identical(audit$provenance, pathwayIndex$provenance)) {
        stop(
            "Stored pathway scores use a different pathway_index.",
            call. = FALSE
        )
    }
    if (!returnTop) {
        attr(result, "pathway_coverage") <- audit$coverage
        return(result)
    }
    empty <- data.frame(
        metabolite = NA_character_, ramp_id = NA_character_,
        n_sig_path = NA_integer_, max_cor = NA_real_, z_score = NA_real_,
        pval = NA_real_, pval_adj = NA_real_
    )
    top <- do.call(rbind, lapply(result, function(table) {
        if (is.null(table) || !any(is.finite(table$z_score))) {
            return(empty)
        }
        as.data.frame(table[which(is.finite(table$z_score))[1L], ])
    }))
    metadata <- .featureMetadata(data, mzAssay)
    metadata$mz_names <- ids
    metadata$mz <- .massValues(metadata, ids)
    # A repeated call must not create duplicate output column names.
    for (column in names(top)) metadata[[column]] <- top[[column]]
    rownames(metadata) <- ids
    attr(metadata, "pathway_coverage") <- audit$coverage
    metadata
}
