# Development assistance: see inst/COPYRIGHTS.

#### Inspect stored annotation candidates ######################################

# The row key survives ordinary SummarizedExperiment subsetting and renaming.
.annotationFeatureIds <- function(metadata, initialize = FALSE) {
    ids <- metadata[[".spamtp_feature_id"]]
    if (is.null(ids) && initialize) {
        ids <- as.character(seq_len(nrow(metadata)))
    }
    if (!is.character(ids) || length(ids) != nrow(metadata) ||
        anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop("Annotation feature keys must be present and unique. ",
            "Restore rowData(x)$.spamtp_feature_id, or remove it and re-run ",
            "annotateSM() or addFMP10Annotations().",
            call. = FALSE
        )
    }
    ids
}

# Current display names for stored candidates; NA marks removed features.
.candidateFeatures <- function(x, candidates) {
    if (!nrow(x) || !nrow(candidates)) {
        return(rep(NA_character_, nrow(candidates)))
    }
    labels <- rownames(x)
    if (is.null(labels) || anyNA(labels) || any(!nzchar(labels)) ||
        anyDuplicated(labels)) {
        stop("Use unique, non-empty feature names to inspect annotations.",
            call. = FALSE
        )
    }
    metadata <- as.data.frame(SummarizedExperiment::rowData(x))
    if ("feature_id" %in% names(candidates)) {
        ids <- .annotationFeatureIds(metadata)
        if (anyNA(candidates$feature_id) ||
            any(!nzchar(candidates$feature_id))) {
            stop("Stored annotation feature keys are missing. ",
                "Re-run annotation.",
                call. = FALSE
            )
        }
        return(rownames(x)[match(candidates$feature_id, ids)])
    }
    # Older stores have no preserved row key. Recover only unambiguous masses;
    # equal-mass rows can exchange places without changing their mass vector.
    declared <- attr(candidates, "annotation_metadata")$raw_mz_column
    masses <- lapply(
        intersect(c(declared, "raw_mz", "mz"), names(metadata)),
        function(column) as.numeric(metadata[[column]])
    )
    if (!length(masses)) {
        masses <- list(.massValues(metadata, rownames(x)))
    }
    positions <- lapply(masses, function(m) match(candidates$observed_mz, m))
    matched <- vapply(positions, function(p) sum(!is.na(p)), integer(1))
    best <- which.max(matched)
    position <- positions[[best]]
    conflicting <- vapply(positions, function(other) {
        any(!is.na(position) & !is.na(other) & position != other)
    }, logical(1))
    mass <- masses[[best]]
    duplicated_mass <- mass[
        duplicated(mass) | duplicated(mass, fromLast = TRUE)
    ]
    identity <- intersect(c("feature", "ID"), names(candidates))
    if (length(identity)) {
        pairs <- unique(candidates[, c(identity[1], "observed_mz")])
        duplicated_mass <- c(
            duplicated_mass,
            pairs$observed_mz[duplicated(pairs$observed_mz)]
        )
    }
    if (any(conflicting) ||
        any(!is.na(position) & candidates$observed_mz %in% duplicated_mass)) {
        stop("Older annotation records have ambiguous feature masses. ",
            "Re-run annotation to preserve feature identities.",
            call. = FALSE
        )
    }
    rownames(x)[position]
}

# All candidate-based summaries use the same retained-feature view. The full
# annotation store remains available as provenance in metadata(x).
.currentAnnotationCandidates <- function(x) {
    declared <- .storedData(x, "mz_annotation")$metadata$assay
    if (inherits(x, "SingleCellExperiment") &&
        length(declared) == 1L &&
        declared %in% SingleCellExperiment::altExpNames(x)) {
        x <- .pathwayAnnotationObject(x, "main")
    }
    candidates <- .get_mz_annotations(x)
    candidates$feature <- .candidateFeatures(x, candidates)
    candidates[!is.na(candidates$feature), , drop = FALSE]
}

# Distinct RaMP compound identities per feature among the supplied candidates.
.candidateIdentityCounts <- function(feature, rampIds) {
    ids <- strsplit(rampIds, ";[[:space:]]*")
    pairs <- unique(data.frame(
        feature = rep(feature, lengths(ids)),
        compound = unlist(ids, use.names = FALSE)
    ))
    pairs <- pairs[grepl("^RAMP_C_", pairs$compound), , drop = FALSE]
    table(factor(pairs$feature, levels = unique(feature)))
}

#' Ranked annotation candidates for measured features
#'
#' Returns the candidate identities stored by [annotateSM()] (with
#' `save.intermediate = TRUE`), ranked within each feature by decreasing score
#' and then increasing absolute ppm error. Candidates remain hypotheses: the
#' `identities` column counts the distinct RaMP compound identities that remain
#' for each feature at `minScore`, so values above one mark ambiguous features.
#' After row subsetting, only candidates for retained features are returned.
#' New annotations keep a unique `.spamtp_feature_id` column in `rowData(x)`;
#' preserve this column when editing feature metadata. Candidates follow their
#' data rows when rows are reordered or renamed. Older stores without this key
#' are matched by observed mass and require reannotation if those masses
#' are ambiguous.
#' Display names must remain unique and non-empty.
#' For a container with alternative experiments, these helpers inspect the
#' primary experiment's annotation store.
#'
#' @param x A SummarizedExperiment-derived object annotated with [annotateSM()].
#' @param features Optional feature names, or numeric m/z values matched to
#'   the nearest measured feature with [findNearestMZ()]. NULL returns all
#'   annotated features.
#' @param minScore Minimum candidate score retained.
#' @param n Maximum number of candidates returned per feature.
#' @param distinct Keep only the best-ranked row for each combination of RaMP
#'   identities and adduct within a feature.
#'
#' @return A data frame with `feature`, `rank` and `identities` followed by
#'   the stored candidate columns (for example `observed_mz`, `IsomerNames`,
#'   `Adduct`, `Error` in ppm, `Score` and `Ramp_IDs`), ordered by feature and
#'   rank.
#' @seealso [summarizeAnnotationScores()], [plotAnnotationScores()]
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' annotationCandidates(spe, features = "ion1", n = 3)
annotationCandidates <- function(
    x, features = NULL, minScore = 0, n = Inf,
    distinct = TRUE
) {
    .requireExperiment(x, "SummarizedExperiment")
    if (length(minScore) != 1L || !is.finite(minScore)) {
        stop("minScore must be one finite number.", call. = FALSE)
    }
    if (length(n) != 1L || is.na(n) || n < 1) {
        stop("n must be a positive number.", call. = FALSE)
    }
    candidates <- .currentAnnotationCandidates(x)
    if (!is.null(features)) {
        if (is.numeric(features)) {
            features <- findNearestMZ(x, features)
        }
        unknown <- setdiff(features, rownames(x))
        if (length(unknown)) {
            stop(
                "Unknown features: ",
                paste(unknown, collapse = ", "),
                call. = FALSE
            )
        }
        candidates <- candidates[candidates$feature %in% features, , drop = FALSE]
    }
    featureLevels <- if (is.null(features)) rownames(x) else unique(features)
    candidates <- candidates[candidates$Score >= minScore, , drop = FALSE]
    candidates <- candidates[order(
        match(candidates$feature, featureLevels),
        -candidates$Score, abs(candidates$Error)
    ), , drop = FALSE]
    if (distinct) {
        candidates <- candidates[!duplicated(candidates[, c("feature", "Ramp_IDs", "Adduct")]), ,
            drop = FALSE
        ]
    }
    candidates$rank <- stats::ave(seq_len(nrow(candidates)), candidates$feature,
        FUN = seq_along
    )
    identities <- .candidateIdentityCounts(
        candidates$feature,
        candidates$Ramp_IDs
    )
    candidates$identities <- as.integer(identities[candidates$feature])
    candidates <- candidates[candidates$rank <= n, , drop = FALSE]
    first <- c("feature", "rank", "identities")
    candidates <- candidates[, c(first, setdiff(names(candidates), first)), drop = FALSE]
    rownames(candidates) <- NULL
    candidates
}

#' Summarize annotation candidates across score thresholds
#'
#' For each score threshold, counts the retained candidates, the features with
#' at least one candidate, and how many of those features have exactly one or
#' several distinct RaMP compound identities.
#' Counts refer only to features currently present in `x`.
#'
#' @param x A SummarizedExperiment-derived object annotated with [annotateSM()].
#' @param thresholds Numeric score thresholds.
#'
#' @return A data frame with one row per threshold and columns
#'   `score_threshold`, `candidates`, `matched_features`,
#'   `single_identity_features` and `ambiguous_features`.
#' @seealso [annotationCandidates()]
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' summarizeAnnotationScores(spe, thresholds = c(0, 0.2, 0.5))
summarizeAnnotationScores <- function(x, thresholds = c(0, 0.05, 0.2, 0.5)) {
    .requireExperiment(x, "SummarizedExperiment")
    if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds)) {
        stop("thresholds must be numeric.", call. = FALSE)
    }
    candidates <- .currentAnnotationCandidates(x)
    rows <- lapply(thresholds, function(threshold) {
        retained <- candidates[candidates$Score >= threshold, , drop = FALSE]
        identities <- .candidateIdentityCounts(
            retained$feature,
            retained$Ramp_IDs
        )
        data.frame(
            score_threshold = threshold,
            candidates = nrow(retained),
            matched_features = length(unique(retained$feature)),
            single_identity_features = sum(identities == 1L),
            ambiguous_features = sum(identities > 1L)
        )
    })
    do.call(rbind, rows)
}

#' Plot the best annotation score of each measured feature
#'
#' Draws a histogram of the highest candidate score per annotated feature,
#' with a dashed line at a chosen score threshold. Features without a
#' ppm-valid candidate are not shown.
#' Candidates for features removed by row subsetting are excluded.
#'
#' @param x A SummarizedExperiment-derived object annotated with [annotateSM()].
#' @param threshold Score threshold to mark, or NULL for none.
#' @param binwidth Histogram bin width on the score scale.
#'
#' @return A `ggplot`.
#' @seealso [summarizeAnnotationScores()]
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' plotAnnotationScores(spe, threshold = 0.2)
plotAnnotationScores <- function(x, threshold = 0.05, binwidth = 0.025) {
    .requireExperiment(x, "SummarizedExperiment")
    candidates <- .currentAnnotationCandidates(x)
    best <- data.frame(score = as.numeric(tapply(
        candidates$Score, candidates$feature, max
    )))
    plot <- ggplot2::ggplot(best, ggplot2::aes(x = .data$score)) +
        ggplot2::geom_histogram(
            binwidth = binwidth, boundary = 0,
            fill = "#567D9A", colour = "white"
        ) +
        ggplot2::labs(
            x = "Best candidate score per feature",
            y = "Mass features"
        ) +
        ggplot2::theme_classic(base_size = 11)
    if (!is.null(threshold)) {
        plot <- plot + ggplot2::geom_vline(
            xintercept = threshold, linetype = 2,
            colour = "#8A4D34"
        )
    }
    plot
}
