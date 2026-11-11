.pathway_links <- function(table) {
    required <- c("rampId", "pathwayRampId")
    if (is.null(table)) {
        return(data.frame(rampId = character(), pathwayRampId = character()))
    }
    if (!all(required %in% names(table))) stop(
        "Pathway links require rampId and pathwayRampId.",
        call. = FALSE
    )
    table <- as.data.frame(table)[, required, drop = FALSE]
    table[] <- lapply(table, as.character)
    table <- unique(table[!is.na(table$rampId) & nzchar(table$rampId) &
        !is.na(table$pathwayRampId) & nzchar(table$pathwayRampId), ])
    table <- table[order(table$pathwayRampId, table$rampId, method = "radix"), ]
    rownames(table) <- NULL
    table
}

.pathway_resource_fingerprint <- function(resources) {
    meta <- resources$pathway
    fields <- intersect(
        c(
            "pathwayRampId",
            "pathwayName",
            "sourceId",
            "type",
            "pathwayCategory"
        ),
        names(meta)
    )
    meta <- as.data.frame(meta)[, fields, drop = FALSE]
    meta[] <- lapply(meta, as.character)
    if (nrow(meta)) meta <- unique(meta[do.call(
        order,
        c(unname(meta), list(na.last = TRUE, method = "radix"))
    ), , drop = FALSE])
    rownames(meta) <- NULL
    rlang::hash(list(
        links = .pathway_links(resources$analytehaspathway),
        metadata = meta
    ))
}

.pathway_raw_fingerprint <- function(resources) {
    rlang::hash(resources[c("analytehaspathway", "pathway")])
}

.pathway_validate_index <- function(
    index,
    resources = NULL,
    geneIndex = NULL
) {
    if (!inherits(index, "spamtp_pathway_index") || !identical(index$provenance$schema_version, 1L)) {
        stop(
            "pathwayIndex must be built by buildPathwayIndex().",
            call. = FALSE
        )
    }
    if (!is.null(resources)) {
        # Identical inputs need no sorting or duplicate removal. Changed inputs
        # still use the order-independent identity check, including old indexes
        # without a recorded raw fingerprint.
        unchanged <- identical(
            .pathway_raw_fingerprint(resources),
            index$provenance$raw_membership_fingerprint
        )
        if (!unchanged && !identical(
            .pathway_resource_fingerprint(resources),
            index$provenance$membership_fingerprint
        )) {
            stop(
                "pathwayIndex does not match pathway memberships or metadata; rebuild it.",
                call. = FALSE
            )
        }
        if (!is.null(index$gene_index) && !is.null(resources$source_df)) {
            .gene_check_source_index(resources$source_df, index$gene_index)
        }
    }
    if (!is.null(geneIndex) && !identical(geneIndex, index$gene_index)) {
        stop(
            "geneIndex and pathwayIndex use different gene identities.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.pathway_database <- function(database, pathwayIndex) {
    if (is.null(pathwayIndex)) {
        return(database)
    }
    .pathway_validate_index(pathwayIndex)
    if (!is.null(database)) {
        return(database)
    }
    database <- pathwayIndex$raw_resources
    if (is.null(database$source_df)) {
        database$source_df <- data.frame(
            rampId = character(), sourceId = character()
        )
    }
    database
}

.pathway_context <- function(
    resources, database = resources,
    geneMapping = "auto", geneReference = NULL, geneIndex = NULL,
    geneReferenceVersion = "latest", geneReferenceLocalDir = NULL,
    organism = "Homo sapiens", pathwayIndex = NULL
) {
    geneMapping <- match.arg(geneMapping, c("auto", "hgnc", "ramp"))
    if (!is.null(pathwayIndex)) {
        .pathway_validate_index(pathwayIndex, resources, geneIndex)
        if (!is.null(pathwayIndex$gene_index) &&
            !identical(
                organism,
                pathwayIndex$gene_index$provenance$organism
            )) {
            stop(
                "pathwayIndex organism does not match organism.",
                call. = FALSE
            )
        }
        if (identical(geneMapping, "ramp") && !is.null(pathwayIndex$gene_index)) {
            stop(
                "geneMapping = 'ramp' conflicts with the HGNC pathway_index.",
                call. = FALSE
            )
        }
        if (identical(geneMapping, "hgnc") && is.null(pathwayIndex$gene_index)) {
            stop(
                "geneMapping = 'hgnc' requires an HGNC pathway_index.",
                call. = FALSE
            )
        }
        resources[names(pathwayIndex$resources)] <- pathwayIndex$resources
        return(list(
            resources = resources,
            index = pathwayIndex$gene_index, pathway_index = pathwayIndex
        ))
    }
    raw <- resources
    if (is.null(resources$source_df)) {
        if (!is.null(geneIndex) || !is.null(geneReference) || identical(geneMapping, "hgnc")) {
            stop(
                "Gene harmonization requires sourceDf or a prebuilt pathway_index.",
                call. = FALSE
            )
        }
        view <- list(resources = resources, index = NULL)
    } else {
        view <- .gene_pathway_view(
            resources, database, geneMapping, geneReference,
            geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism
        )
    }
    links <- .pathway_links(view$resources$analytehaspathway)
    raw_links <- .pathway_links(raw$analytehaspathway)
    ids <- sort(unique(c(
        as.character(raw$pathway$pathwayRampId),
        raw_links$pathwayRampId
    )))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    sets <- function(table) {
        answer <- split(table$rampId, table$pathwayRampId)
        answer <- answer[match(ids, names(answer))]
        answer[lengths(answer) == 0L] <- list(character())
        stats::setNames(answer, ids)
    }
    meta <- as.data.frame(raw$pathway)
    if (!"pathwayRampId" %in% names(meta)) meta <- data.frame(pathwayRampId = character())
    meta <- meta[match(ids, meta$pathwayRampId), , drop = FALSE]
    meta$pathwayRampId <- ids
    if (!"pathwayName" %in% names(meta)) meta$pathwayName <- ids
    absent <- is.na(meta$pathwayName) | !nzchar(meta$pathwayName)
    meta$pathwayName[absent] <- ids[absent]
    rownames(meta) <- NULL
    conflicts <- if (is.null(view$index)) {
        character()
    } else {
        view$index$nodes$rampId[view$index$nodes$status == "conflicting_record"]
    }
    index <- structure(list(
        members = sets(links), raw_members = sets(raw_links), metadata = meta,
        conflicting_members = sets(raw_links[raw_links$rampId %in% conflicts, ]),
        gene_index = view$index,
        resources = view$resources[intersect(
            c("source_df", "analyte", "analytehaspathway", "pathway"),
            names(view$resources)
        )],
        raw_resources = raw[intersect(
            c("source_df", "analyte", "analytehaspathway", "pathway"),
            names(raw)
        )],
        provenance = list(
            schema_version = 1L,
            membership_fingerprint = .pathway_resource_fingerprint(raw),
            raw_membership_fingerprint = .pathway_raw_fingerprint(raw),
            gene_mapping = if (is.null(view$index)) "ramp" else "hgnc",
            gene = view$index$provenance,
            resources = lapply(
                raw[intersect(
                    c("source_df", "analytehaspathway", "pathway"),
                    names(raw)
                )],
                attr,
                which = "spamtp_database"
            )
        )
    ), class = "spamtp_pathway_index")
    list(resources = view$resources, index = view$index, pathway_index = index)
}

#' Build a shared, auditable pathway membership index
#'
#' Unites pathway memberships of non-conflicting RaMP records belonging to one
#' HGNC gene. Retains raw memberships and excluded conflicting records. Reuse
#' this index in enrichment, expression aggregation, scoring and named pathway
#' plots to obtain the same identities and membership sets. Pathways are keyed
#' by pathwayRampId; identical display names never merge distinct pathways.
#'
#' @param database A resource bundle with analytehaspathway and pathway, and
#'   source_df for gene harmonization. NULL retrieves these from SpaMTPdb.
#' @inheritParams fishersPathwayAnalysis
#' @inheritParams buildGeneMappingIndex
#' @return A spamtp_pathway_index containing canonical and raw member sets,
#'   conflicting members, metadata, gene_index and resource fingerprints.
#'   Database size means unique analyte identities after reconciliation and
#'   conflict exclusion, before intersecting with measured features.
#' @export
#' @examples
#' db <- list(analytehaspathway = data.frame(
#'     rampId = "RAMP_C_1",
#'     pathwayRampId = "P1"
#' ), pathway = data.frame(pathwayRampId = "P1"))
#' buildPathwayIndex(db)$members
buildPathwayIndex <- function(
    database = NULL, geneMapping = c("auto", "hgnc", "ramp"),
    geneReference = NULL, geneIndex = NULL, geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL, organism = "Homo sapiens",
    databaseVersion = "latest", databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL
) {
    needed <- c("analytehaspathway", "pathway")
    if (is.null(database) || "source_df" %in% names(database)) needed <- c(needed, "source_df")
    resources <- .spamtp_db_bundle(needed, database,
        version = databaseVersion,
        source = match.arg(databaseSource), localDir = databaseLocalDir
    )
    if (!is.null(database$analyte)) resources$analyte <- database$analyte
    if (is.null(resources$analyte)) {
        resources$analyte <- data.frame(
            rampId = unique(c(
                resources$source_df$rampId,
                resources$analytehaspathway$rampId
            ))
        )
    }
    .pathway_context(
        resources, database, match.arg(geneMapping), geneReference,
        geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism
    )$pathway_index
}

.pathway_coverage <- function(
    index, measuredIds, usedIds = measuredIds,
    foregroundIds = character(), types = c("G", "C"), minSize = 1, maxSize = Inf
) {
    members <- .filterPathwayMembers(index$members, types = types)
    raw <- .filterPathwayMembers(index$raw_members, types = types)
    measured <- .filterPathwayMembers(members, measuredIds)
    used <- .filterPathwayMembers(measured, usedIds)
    foreground <- .filterPathwayMembers(measured, foregroundIds)
    coverage <- index$metadata
    coverage$database_size <- lengths(members)
    coverage$database_raw_size <- lengths(raw)
    coverage$measured_size <- lengths(measured)
    coverage$coverage_fraction <- ifelse(
        lengths(members) > 0,
        lengths(measured) / lengths(members),
        NA_real_
    )
    coverage$used_size <- lengths(used)
    coverage$eligible <- lengths(used) >= minSize & lengths(used) <= maxSize & lengths(used) > 0
    coverage$database_members <- I(unname(members))
    coverage$measured_members <- I(unname(measured))
    coverage$used_members <- I(unname(used))
    coverage$foreground_members <- I(unname(foreground))
    coverage$excluded_conflict_records <- I(unname(
        .filterPathwayMembers(index$conflicting_members, types = types)
    ))
    coverage$excluded_conflict_count <- lengths(coverage$excluded_conflict_records)
    attr(coverage, "pathway_index") <- index$provenance
    coverage
}

# Match the complete membership vector once instead of rebuilding a hash of
# measured identifiers for every pathway. Empty pathways retain their positions.
.filterPathwayMembers <- function(sets, allowed = NULL, types = NULL) {
    if (!length(sets)) {
        return(sets)
    }
    values <- unlist(sets, use.names = FALSE)
    keep <- rep(TRUE, length(values))
    if (!is.null(types)) {
        keep <- grepl(
            paste0("^RAMP_(", paste(types, collapse = "|"), ")_"), values
        )
    }
    if (!is.null(allowed)) keep <- keep & values %in% allowed
    group <- factor(rep(seq_along(sets), lengths(sets))[keep],
        levels = seq_along(sets)
    )
    stats::setNames(unname(split(values[keep], group)), names(sets))
}

.pathway_sets_from_coverage <- function(coverage, eligibleOnly = TRUE) {
    rows <- if (eligibleOnly) which(coverage$eligible) else seq_len(nrow(coverage))
    stats::setNames(
        unclass(coverage$used_members[rows]),
        coverage$pathwayRampId[rows]
    )
}

.pathway_select <- function(index, pathways) {
    if (!is.character(pathways) || !length(pathways) || anyNA(pathways)) {
        stop(
            "pathways must contain pathway IDs or unambiguous names.",
            call. = FALSE
        )
    }
    vapply(pathways, function(query) {
        if (query %in% index$metadata$pathwayRampId) {
            return(query)
        }
        hit <- index$metadata$pathwayRampId[tolower(index$metadata$pathwayName) == tolower(query)]
        if (length(hit) != 1L) {
            stop("Unknown or ambiguous pathway name: ", query,
                ". Supply an exact pathwayRampId.",
                call. = FALSE
            )
        }
        hit
    }, character(1), USE.NAMES = FALSE)
}

.pathway_expression <- function(
    expression, index, duplicateGenes = c("error", "mean", "sum"),
    featureIds = rownames(expression), verbose = TRUE
) {
    duplicateGenes <- match.arg(duplicateGenes)
    if (is.null(featureIds) || length(featureIds) != nrow(expression) || anyNA(featureIds) ||
        any(!nzchar(featureIds))) {
        stop("Expression features need non-missing identifiers.", call. = FALSE)
    }
    ids <- as.character(featureIds)
    canonical <- ids
    genes <- !grepl("^RAMP_C_", ids, ignore.case = TRUE)
    audit <- NULL
    if (!is.null(index$gene_index) && any(genes)) {
        audit <- .map_gene_index(
            ids[genes],
            index$gene_index,
            verbose = verbose
        )
        canonical[genes] <- audit$canonical_ramp_id
    }
    if (is.null(index$gene_index) && any(genes) && !is.null(index$raw_resources$source_df)) {
        source <- index$raw_resources$source_df
        if (!"sourceId" %in% names(source)) source$sourceId <- NA_character_
        mapped <- .fisher_map_ids(
            ids,
            "G",
            source,
            unique(c(source$rampId, unlist(index$members)))
        )
        by_input <- split(mapped$rampId, mapped$input)
        canonical[genes] <- NA_character_
        for (i in which(genes)) {
            match_ids <- unique(by_input[[tolower(ids[i])]])
            if (length(match_ids) > 1L) stop(
                "Ambiguous RaMP-only expression mapping: ",
                ids[i],
                call. = FALSE
            )
            if (length(match_ids)) canonical[i] <- match_ids
        }
        audit <- data.frame(
            input = ids[genes], canonical_ramp_id = canonical[genes],
            status = ifelse(
                is.na(canonical[genes]),
                "unmapped",
                "mapped_ramp_only"
            ),
            method = "ramp_source"
        )
        if (verbose && anyNA(canonical[genes])) {
            warning(
                sum(is.na(canonical[genes])), " gene input(s) unresolved in RaMP-only mapping; inspect the mapping audit.",
                call. = FALSE
            )
        }
    }
    keep <- !is.na(canonical)
    values <- expression[keep, , drop = FALSE]
    canonical <- canonical[keep]
    if (!length(canonical)) stop(
        "No expression features have usable identities.",
        call. = FALSE
    )
    if (!.allFiniteAssay(values)) stop(
        "Pathway expression requires finite values.",
        call. = FALSE
    )
    groups <- split(
        seq_along(canonical),
        factor(canonical, levels = unique(canonical))
    )
    choose <- lapply(groups, function(rows) {
        if (length(rows) > 1L && duplicateGenes == "error") {
            if (!all(vapply(rows[-1L], function(i) {
                isTRUE(all.equal(
                    as.numeric(values[i, ]), as.numeric(values[rows[1L], ]),
                    tolerance = 0
                ))
            }, logical(1)))) {
                stop("Multiple expression features map to ", canonical[rows[1L]],
                    "; choose duplicateGenes = 'mean' or 'sum' after reviewing the assay scale.",
                    call. = FALSE
                )
            }
            rows <- rows[1L]
        }
        rows
    })
    j <- unlist(choose, use.names = FALSE)
    i <- rep(seq_along(choose), lengths(choose))
    weight <- if (duplicateGenes == "mean") rep(1 / lengths(choose), lengths(choose)) else rep(1, length(j))
    W <- Matrix::sparseMatrix(
        i = i,
        j = j,
        x = weight,
        dims = c(length(groups), nrow(values))
    )
    result <- W %*% values
    dimnames(result) <- list(names(groups), colnames(expression))
    list(
        expression = result, mapping = audit,
        original_features = stats::setNames(
            lapply(groups, function(rows) ids[keep][rows]),
            names(groups)
        ),
        duplicate_genes = duplicateGenes
    )
}

.pathway_gene_ids <- function(ids, index) {
    if (!is.null(index$gene_index)) {
        mapped <- .map_gene_index(ids, index$gene_index, verbose = FALSE)$canonical_ramp_id
        return(unique(mapped[!is.na(mapped)]))
    }
    source <- index$raw_resources$source_df
    if (is.null(source)) {
        return(unique(ids))
    }
    if (!"sourceId" %in% names(source)) source$sourceId <- NA_character_
    unique(.fisher_map_ids(
        ids,
        "G",
        source,
        unique(c(source$rampId, unlist(index$members)))
    )$rampId)
}

.pathway_score_matrix <- function(expression, sets, standardize = TRUE) {
    scores <- lapply(sets, function(ids) {
        ids <- intersect(ids, rownames(expression))
        if (!length(ids)) {
            return(rep(NA_real_, ncol(expression)))
        }
        value <- as.numeric(Matrix::colSums(expression[ids, , drop = FALSE])) / sqrt(length(ids))
        value <- value - mean(value)
        if (standardize) {
            deviation <- stats::sd(value)
            value <- if (is.finite(deviation) && deviation > 0) value / deviation else rep(0, length(value))
        }
        value
    })
    if (!length(scores)) {
        return(matrix(
            numeric(),
            0L,
            ncol(expression),
            dimnames = list(character(), colnames(expression))
        ))
    }
    answer <- do.call(rbind, scores)
    dimnames(answer) <- list(names(sets), colnames(expression))
    answer
}

.pathway_score_context <- function(
    object, assay, assayName, database, pathwayIndex,
    geneMapping, geneReference, geneIndex, geneReferenceVersion,
    geneReferenceLocalDir, organism, databaseVersion, databaseSource,
    databaseLocalDir, duplicateGenes, minPathSize = 1, maxPathSize = Inf
) {
    # Validate the assay before touching external resources.
    E <- .assayData(object, assay, assayName)
    database <- .pathway_database(database, pathwayIndex)
    if (is.null(pathwayIndex)) {
        pathwayIndex <- buildPathwayIndex(
            database, geneMapping, geneReference,
            geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism,
            databaseVersion, databaseSource, databaseLocalDir
        )
    } else {
        .pathway_context(
            database, database, geneMapping, geneReference, geneIndex,
            geneReferenceVersion, geneReferenceLocalDir, organism, pathwayIndex
        )
    }
    if (!is.null(pathwayIndex$gene_index)) .gene_check_experiment_species(object, assay, organism)
    prior <- S4Vectors::metadata(.experimentForAssay(object, assay))$pathway_mapping$provenance
    if (!is.null(prior) && identical(prior$gene_mapping, "hgnc") &&
        !identical(prior, pathwayIndex$provenance)) {
        stop(
            "The input assay was mapped with a different pathway index. Rebuild it from original features.",
            call. = FALSE
        )
    }
    mapped <- .pathway_expression(E, pathwayIndex, duplicateGenes)
    types <- if (all(grepl("^RAMP_[GC]_", rownames(mapped$expression)))) {
        unique(sub("^RAMP_([GC])_.*", "\\1", rownames(mapped$expression)))
    } else {
        NULL
    }
    if (length(minPathSize) != 1L || !is.finite(minPathSize) || minPathSize < 1 ||
        minPathSize != floor(minPathSize) || length(maxPathSize) != 1L ||
        is.na(maxPathSize) || maxPathSize < minPathSize) {
        stop("Invalid pathway size limits.", call. = FALSE)
    }
    coverage <- .pathway_coverage(pathwayIndex, rownames(mapped$expression),
        types = types,
        minSize = minPathSize, maxSize = maxPathSize
    )
    list(
        index = pathwayIndex, expression = mapped$expression, coverage = coverage,
        audit = list(
            provenance = pathwayIndex$provenance, inputs = mapped$mapping,
            original_features = mapped$original_features, duplicate_genes = mapped$duplicate_genes,
            coverage = coverage
        )
    )
}
