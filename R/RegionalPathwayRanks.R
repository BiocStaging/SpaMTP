.regionalPathwaysFromRanks <- function(
    ranks, index, minSize, maxSize, rankScale,
    nPermSimple, regions = NULL,
    BPPARAM = BiocParallel::SerialParam()
) {
    .pathway_validate_index(index)
    .validateFeatureCount(nPermSimple, "nPermSimple")
    if (!is.list(ranks) || !length(ranks) || is.null(names(ranks)) ||
        anyDuplicated(names(ranks)) || any(!names(ranks) %in% c("genes", "metabolites"))) {
        stop(
            "ranks must be a named list of genes/metabolites feature-by-region matrices.",
            call. = FALSE
        )
    }
    selected_regions <- regions
    mapped <- list()
    mapping <- list()
    regions <- NULL
    for (type in names(ranks)) {
        E <- ranks[[type]]
        if (length(dim(E)) != 2L || !is.numeric(E) || is.null(rownames(E)) ||
            is.null(colnames(E)) || anyNA(rownames(E)) || anyNA(colnames(E)) ||
            anyDuplicated(rownames(E)) || anyDuplicated(colnames(E))) {
            stop(
                "Rank matrices require distinct feature and region names.",
                call. = FALSE
            )
        }
        if (is.null(regions)) regions <- colnames(E)
        if (!setequal(regions, colnames(E))) stop(
            "Rank matrices must describe the same regions.",
            call. = FALSE
        )
        if (type == "metabolites" && !all(grepl("^RAMP_C_", rownames(E)))) {
            stop(
                "Metabolite ranks require explicit RaMP compound identities; use createPathwayAssay first.",
                call. = FALSE
            )
        }
        m <- .pathway_expression(E[, regions, drop = FALSE], index, "error")
        ids <- rownames(m$expression)
        prefix <- if (type == "genes") "^RAMP_G_" else "^RAMP_C_"
        if (any(!grepl(prefix, ids))) stop(
            "Rank identities do not match analyte type.",
            call. = FALSE
        )
        mapped[[type]] <- as.matrix(m$expression)
        mapping[[type]] <- m$mapping
    }
    if (!is.null(selected_regions)) {
        if (!length(selected_regions) || !all(selected_regions %in% regions)) {
            stop(
                "regions must name columns of the rank matrices.",
                call. = FALSE
            )
        }
        regions <- unique(selected_regions)
    }
    measured <- unique(unlist(lapply(mapped, rownames)))
    coverage <- used_ranks <- tables <- list()
    for (region in regions) {
        blocks <- lapply(mapped, function(expression) {
            v <- stats::setNames(expression[, region], rownames(expression))
            v <- v[is.finite(v)]
            if (rankScale == "rms" && length(v)) {
                divisor <- sqrt(mean(v^2))
                if (divisor > 0) v <- v / divisor
            }
            v
        })
        stat <- unlist(blocks, use.names = FALSE)
        names(stat) <- unlist(lapply(blocks, names), use.names = FALSE)
        if (anyDuplicated(names(stat))) stop(
            "An identity was repeated across rank blocks.",
            call. = FALSE
        )
        used_ranks[[region]] <- sort(stat, decreasing = TRUE)
        c <- .pathway_coverage(index, measured,
            usedIds = names(stat),
            types = ifelse(names(ranks) == "genes", "G", "C"), minSize = minSize, maxSize = maxSize
        )
        coverage[[region]] <- c
        sets <- .pathway_sets_from_coverage(c)
        if (!length(sets) || length(stat) < 2L || all(stat == 0)) next
        t <- as.data.frame(fgsea::fgseaMultilevel(
            pathways = sets, stats = used_ranks[[region]],
            minSize = minSize, maxSize = maxSize, nPermSimple = nPermSimple,
            eps = 1e-10, BPPARAM = BPPARAM
        ))
        if (!nrow(t)) next
        # Include all eligible pathways in the declared adjustment family, including
        # pathways for which the Monte Carlo algorithm could not estimate a P value.
        t$padj <- stats::p.adjust(t$pval, method = "BH", n = length(sets))
        names(t)[names(t) == "pathway"] <- "pathwayRampId"
        t$Cluster_id <- region
        t$tested_family <- length(sets)
        members <- c$used_members[match(t$pathwayRampId, c$pathwayRampId)]
        t$measured_genes <- vapply(
            members,
            function(ids) sum(grepl("^RAMP_G_", ids)),
            integer(1)
        )
        t$measured_metabolites <- vapply(
            members,
            function(ids) sum(grepl("^RAMP_C_", ids)),
            integer(1)
        )
        t$leading_genes <- vapply(
            t$leadingEdge,
            function(ids) sum(grepl("^RAMP_G_", ids)),
            integer(1)
        )
        t$leading_metabolites <- vapply(
            t$leadingEdge,
            function(ids) sum(grepl("^RAMP_C_", ids)),
            integer(1)
        )
        t <- merge(
            t,
            index$metadata,
            by = "pathwayRampId",
            all.x = TRUE,
            sort = FALSE
        )
        tables[[region]] <- t
    }
    answer <- if (length(tables)) {
        do.call(rbind, tables)
    } else {
        data.frame(
            pathwayRampId = character(), pathwayName = character(), Cluster_id = character(),
            NES = numeric(), pval = numeric(), padj = numeric(), size = integer()
        )
    }
    rownames(answer) <- NULL
    attr(answer, "pathway_coverage") <- coverage
    attr(answer, "pathway_index") <- index$provenance
    attr(answer, "gene_mapping") <- mapping
    attr(answer, "ranks") <- used_ranks
    attr(answer, "rank_analysis") <- list(
        input = "complete unfiltered feature rankings",
        rank_scale = rankScale, analyte_types = names(ranks), nPermSimple = nPermSimple,
        family = "All measured/used-size-eligible pathways, BH within region and configured modality combination",
        interpretation = "Competitive feature-set enrichment conditional on supplied ranks; P values do not test between-specimen regional reproducibility. Joint ranks are RMS-scaled per modality; pathways with many measured genes can still be gene-dominated."
    )
    answer
}

#' Convert regional marker statistics to a rank matrix
#'
#' Arranges one statistic from a [findAllDEMs()] result as a feature-by-region
#' matrix, the complete-ranks input of [findRegionalPathways()]. Features
#' without a statistic in a region are `NA` and are omitted from that region's
#' ranking.
#'
#' @param markers A [findAllDEMs()] result with `DEMs` (one row per feature
#'   and region, with `gene` and `cluster` columns) and `expression` (feature
#'   by region means, defining the row and column order).
#' @param statistic Column of `markers$DEMs` to arrange, e.g. `"logFC"` or
#'   `"mean_auc"`.
#'
#' @return A numeric matrix with features in rows and regions in columns.
#' @export
#' @examples
#' markers <- list(
#'     expression = matrix(0, 2, 2, dimnames = list(c("f1", "f2"), c("A", "B"))),
#'     DEMs = data.frame(
#'         gene = c("f1", "f2", "f1"), cluster = c("A", "A", "B"),
#'         logFC = c(1.5, -0.5, 0.2)
#'     )
#' )
#' markerRanks(markers)
markerRanks <- function(markers, statistic = "logFC") {
    if (!is.list(markers) || !is.data.frame(markers$DEMs) || is.null(markers$expression)) {
        stop(
            "markers must be a findAllDEMs() result with DEMs and expression.",
            call. = FALSE
        )
    }
    dems <- markers$DEMs
    if (length(statistic) != 1L || !statistic %in% names(dems) || !is.numeric(dems[[statistic]])) {
        stop(
            "statistic must name one numeric column of markers$DEMs.",
            call. = FALSE
        )
    }
    features <- rownames(markers$expression)
    regions <- colnames(markers$expression)
    ranks <- matrix(NA_real_, length(features), length(regions),
        dimnames = list(features, regions)
    )
    row <- match(dems$gene, features)
    column <- match(as.character(dems$cluster), regions)
    keep <- !is.na(row) & !is.na(column)
    ranks[cbind(row[keep], column[keep])] <- dems[[statistic]][keep]
    ranks
}
