#' Fisher Exact Tests for Pathway Enrichment
#'
#' @param analytes Foreground analytes as a named list with elements `genes`,
#'   `metabolites` and/or `mzs`. Identifiers are matched case-insensitively to
#'   source IDs (e.g. `gene_symbol:TP53`, `hmdb:HMDB0000122`), known RaMP IDs,
#'   then common names for otherwise unmatched queries. Duplicate mappings to
#'   the same RaMP ID count once. m/z values may be numeric or numeric strings,
#'   optionally prefixed with `mz-` or `mz_`.
#' @param maxPathSize Maximum pathway size within the analysis universe
#'   (inclusive; default 500).
#' @param minPathSize Minimum pathway size within the analysis universe
#'   (inclusive; positive integer, default 5).
#' @param alternative Fisher test alternative: `"greater"`, `"less"` or
#'   `"two.sided"`.
#' @param pathwayAllInfo Include names, input IDs and m/z adduct information
#'   for foreground members of each pathway. Does not change the tests.
#' @param pvalCutoff Optional BH-adjusted p-value (FDR) cutoff in `[0, 1]`.
#'   Applied only after correcting over all eligible pathways.
#' @param verbose Whether to print progress messages. Mapping warnings are
#'   always emitted, including when `verbose = FALSE`.
#' @param database Optional named resource list from [loadSpaMTPDatabase()].
#'   Identifier analysis requires `source_df`, `analyte`, `analytehaspathway`
#'   and `pathway`; m/z analysis additionally uses `chem_props` unless `db`
#'   or `index` is supplied through `...`.
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param universe Optional measured background, in the same named-list format
#'   as `analytes`. Supply all analytes eligible for foreground selection, not
#'   just significant hits. Foreground and universe must specify the same
#'   biological modalities (genes and/or compounds; `mzs` and `metabolites`
#'   both denote compounds). All mapped foreground IDs must belong to the
#'   mapped universe, otherwise an error is raised. Mapped background IDs
#'   without any pathway membership still contribute to the universe size.
#'   When `NULL`, use all pathway-linked IDs of the specified modalities in
#'   the database, before pathway-size filtering; mapped foreground IDs without
#'   pathway membership are then excluded with a warning.
#' @param geneMapping `"auto"` uses HGNC identity harmonization for official
#'   SpaMTPdb resources and for custom bundles with `gene_reference` or
#'   `geneIndex`. `"hgnc"` requires harmonization; `"ramp"` explicitly retains
#'   RaMP-only mapping for historical reproduction or a curated non-human
#'   database. Custom bundles without a reference use RaMP-only mapping in auto
#'   mode and never trigger a reference download.
#' @param geneIndex Optional index built from the same source table by
#'   [buildGeneMappingIndex()]. Multiple RaMP records of a resolved HGNC gene
#'   have their pathway memberships united and count once. Conflicting RaMP
#'   records and ambiguous symbols are excluded and recorded in the audit.
#' @inheritParams buildGeneMappingIndex
#' @param pathwayIndex A reusable index from buildPathwayIndex(). When
#'   supplied,
#'   its gene identities and memberships must match any explicit database or
#'   gene_index. Fingerprints prevent incompatible reuse.
#' @param ... Arguments passed to the indexed `annotateTable()` pipeline for
#'   m/z inputs, such as `ppm`, `adducts`, `db` or `index`.
#'
#' @details
#' With universe size U, foreground size K, pathway size M and overlap a, the
#' test table is `matrix(c(a, M - a, K - a, U - K - M + a), nrow = 2)`.
#' Memberships and identifiers are deduplicated by RaMP ID. Pathway sizes are
#' computed after intersecting with the universe; size filtering never changes
#' that universe. All pathways passing the size limits are tested, including
#' zero-overlap pathways. BH correction uses this complete family, before any
#' output cutoff. Missing display metadata does not remove a test.
#'
#' Unmapped input values are excluded with warnings and recorded in the
#' `enrichment` attribute. An empty mapped universe is an error. An empty
#' foreground gives p-values of 1 for all eligible pathways. If no pathways
#' pass the size limits, a zero-row data frame is returned.
#'
#' For a targeted panel, provide the measured panel as `universe`; the default
#' database background is not a substitute for the detection background.
#' m/z foreground and background values are annotated together with identical
#' settings. All mapped candidate RaMP IDs are retained; ambiguous mass matches
#' are not confirmed compound identifications. In a combined analysis, genes
#' and compounds are counted as individual RaMP IDs in a single pooled test.
#'
#' @return A data frame sorted by `p_val`, retaining the columns `pathway_name`,
#'   `pathway_id`, `type`, `pathwayCategory`, `p_val`, `fdr`, `ratio`,
#'   `analytes_in_pathways` (overlap) and `total_in_pathways` (size within the
#'   universe). Also includes `pathwayRampId`, `foreground_analytes_number`
#'   and `background_analytes_number`. `pathwayAllInfo = TRUE` adds member
#'   information. The `enrichment` attribute records mapped foreground and
#'   universe IDs, unmapped inputs, excluded foreground IDs, universe source,
#'   alternative, size limits and the number of tests before output filtering.
#'   Its `geneMapping` entry contains per-input mapping status, gene groups,
#'   excluded conflicting RaMP records and versioned reference provenance.
#' @export
#' @examples
#' data(spamtp_example)
#' db <- spamtp_example$database
#' fishersPathwayAnalysis(list(genes = c("GENE1", "GENE2")),
#'     universe = list(genes = paste0("GENE", 1:8)), database = db,
#'     minPathSize = 2, verbose = FALSE
#' )
fishersPathwayAnalysis <- function(
    analytes,
    maxPathSize = 500,
    minPathSize = 5,
    alternative = "greater",
    pathwayAllInfo = FALSE,
    pvalCutoff = NULL,
    verbose = TRUE,
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    universe = NULL,
    geneMapping = c("auto", "hgnc", "ramp"),
    geneReference = NULL,
    geneIndex = NULL,
    geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL,
    organism = "Homo sapiens",
    pathwayIndex = NULL,
    ...
) {
    analytes <- .fisher_validate_input(analytes, "analytes")
    if (!is.null(universe)) {
        universe <- .fisher_validate_input(universe, "universe")
        if (!setequal(
            .fisher_modalities(analytes),
            .fisher_modalities(universe)
        )) {
            stop("analytes and universe must specify the same biological modalities (genes and/or compounds).",
                call. = FALSE
            )
        }
    }
    alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
    valid_size <- function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
    }
    if (!valid_size(minPathSize) || !valid_size(maxPathSize) ||
        minPathSize > maxPathSize) {
        stop("Pathway size limits must be positive finite integers with minPathSize <= max_path_size.",
            call. = FALSE
        )
    }
    if (!is.null(pvalCutoff) && (!is.numeric(pvalCutoff) ||
        length(pvalCutoff) != 1L || !is.finite(pvalCutoff) ||
        pvalCutoff < 0 || pvalCutoff > 1)) {
        stop(
            "pvalCutoff must be NULL or one finite number in [0, 1].",
            call. = FALSE
        )
    }
    for (flag in c("verbose", "pathwayAllInfo")) {
        value <- get(flag)
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop(flag, " must be TRUE or FALSE.", call. = FALSE)
        }
    }
    args <- list(...)
    if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))) ||
        anyDuplicated(names(args)))) {
        stop(
            "Annotation arguments in ... must have unique names.",
            call. = FALSE
        )
    }
    needed <- c("source_df", "analyte", "analytehaspathway", "pathway")
    if (length(c(analytes$mzs, universe$mzs)) && is.null(args$db) && is.null(args$index)) {
        needed <- c(needed, "chem_props")
    }
    database <- .pathway_database(database, pathwayIndex)
    resources <- .spamtp_db_bundle(
        needed,
        database = database, version = databaseVersion,
        source = match.arg(databaseSource), localDir = databaseLocalDir
    )
    required <- list(
        source_df = c("sourceId", "rampId"), analyte = "rampId",
        analytehaspathway = c("rampId", "pathwayRampId"),
        pathway = "pathwayRampId"
    )
    for (resource in names(required)) {
        if (!all(required[[resource]] %in% names(resources[[resource]]))) {
            stop(resource, " is missing required column(s): ",
                paste(
                    setdiff(required[[resource]], names(resources[[resource]])),
                    collapse = ", "
                ),
                call. = FALSE
            )
        }
    }
    geneMapping <- match.arg(geneMapping)
    gene_view <- .pathway_context(
        resources, database,
        if ("genes" %in% names(analytes) || !is.null(pathwayIndex)) geneMapping else "ramp",
        geneReference, geneIndex, geneReferenceVersion, geneReferenceLocalDir,
        organism, pathwayIndex
    )
    resources <- gene_view$resources
    mapping <- .fisher_map_inputs(
        analytes,
        universe,
        resources,
        args,
        verbose,
        gene_view$index
    )
    foreground_ids <- unique(mapping$mapped$Analyte$rampId)
    links <- as.data.frame(resources$analytehaspathway)
    links <- unique(links[, c("rampId", "pathwayRampId"), drop = FALSE])
    pattern <- paste0(
        "^RAMP_(",
        paste(.fisher_modalities(analytes), collapse = "|"),
        ")_"
    )
    links <- links[!is.na(links$rampId) & grepl(pattern, links$rampId) &
        !is.na(links$pathwayRampId) & nzchar(links$pathwayRampId), , drop = FALSE]
    universe_ids <- if (is.null(universe)) {
        unique(links$rampId)
    } else {
        unique(mapping$mapped$universe$rampId)
    }
    if (!length(universe_ids)) stop("The mapped universe is empty.", call. = FALSE)
    outside <- setdiff(foreground_ids, universe_ids)
    if (length(outside)) {
        if (!is.null(universe)) {
            stop("Mapped foreground IDs are outside universe: ", paste(outside, collapse = ", "),
                ". Supply the complete measured background.",
                call. = FALSE
            )
        }
        warning(length(outside), " mapped foreground ID(s) without pathway membership excluded from the default universe.",
            call. = FALSE
        )
        foreground_ids <- intersect(foreground_ids, universe_ids)
    }
    coverage <- .pathway_coverage(gene_view$pathway_index, universe_ids,
        foregroundIds = foreground_ids, types = .fisher_modalities(analytes),
        minSize = minPathSize, maxSize = maxPathSize
    )
    sets <- .pathway_sets_from_coverage(coverage)
    M <- lengths(sets)
    K <- length(foreground_ids)
    U <- length(universe_ids)
    a <- vapply(sets, function(ids) sum(ids %in% foreground_ids), integer(1))
    verbose_message(paste0(
        "Testing ", length(sets), " pathways; foreground = ", K,
        ", universe = ", U, "."
    ), verbose)
    p <- vapply(seq_along(sets), function(i) {
        stats::fisher.test(
            matrix(c(
                a[i], M[i] - a[i], K - a[i],
                U - K - M[i] + a[i]
            ), nrow = 2),
            alternative = alternative
        )$p.value
    }, numeric(1))
    metadata <- as.data.frame(resources$pathway)
    matched <- match(names(sets), metadata$pathwayRampId)
    field <- function(column, fallback = NA_character_) {
        value <- if (column %in% names(metadata)) {
            as.character(metadata[[column]][matched])
        } else {
            rep(NA_character_, length(sets))
        }
        missing <- is.na(value) | !nzchar(value)
        value[missing] <- rep_len(fallback, length(sets))[missing]
        value
    }
    result <- data.frame(
        pathway_name = field("pathwayName", names(sets)),
        pathway_id = field("sourceId", names(sets)),
        type = field("type"), pathwayCategory = field("pathwayCategory"),
        p_val = p, fdr = stats::p.adjust(p, method = "BH"), ratio = unname(a / M),
        analytes_in_pathways = unname(a), total_in_pathways = unname(M),
        pathwayRampId = names(sets), foreground_analytes_number = rep(K, length(sets)),
        background_analytes_number = rep(U, length(sets)), stringsAsFactors = FALSE
    )
    coverage_fields <- c(
        "database_size", "database_raw_size", "measured_size",
        "coverage_fraction", "used_size", "excluded_conflict_count"
    )
    result <- cbind(
        result,
        coverage[match(result$pathwayRampId, coverage$pathwayRampId), coverage_fields, drop = FALSE]
    )
    if (pathwayAllInfo) {
        result <- cbind(
            result,
            .fisher_pathway_details(
                sets,
                foreground_ids,
                mapping$mapped$Analyte
            )
        )
    }
    if (!is.null(pvalCutoff)) result <- result[result$fdr <= pvalCutoff, , drop = FALSE]
    result <- result[order(result$p_val, result$pathwayRampId), , drop = FALSE]
    rownames(result) <- NULL
    attr(result, "pathway_coverage") <- coverage
    attr(result, "pathway_index") <- gene_view$pathway_index$provenance
    attr(result, "enrichment") <- list(
        universe_source = if (is.null(universe)) "database_pathway_members" else "measured",
        universe_ids = universe_ids, foreground_ids = foreground_ids,
        unmapped = mapping$unmapped, excluded_foreground_ids = outside,
        n_tested = length(sets), alternative = alternative,
        min_path_size = minPathSize, max_path_size = maxPathSize,
        gene_mapping = if (is.null(gene_view$index)) {
            list(mode = "ramp")
        } else {
            list(
                mode = "hgnc", provenance = gene_view$index$provenance, inputs = mapping$gene_audit,
                groups = gene_view$index$nodes[gene_view$index$nodes$canonical_ramp_id %in% universe_ids, ],
                conflicting_records = gene_view$index$nodes[gene_view$index$nodes$status == "conflicting_record", ]
            )
        }
    )
    result
}

.expand_pathway_annotation_ids <- function(db3) {
    if (!"Isomers_IDs" %in% names(db3)) {
        stop("The stored legacy db3 result must contain an Isomers_IDs column.")
    }
    db3 <- tidyr::separate_rows(
        as.data.frame(db3, stringsAsFactors = FALSE),
        dplyr::all_of("Isomers_IDs"),
        sep = "\\s*;\\s*"
    )
    db3$Isomers_IDs <- trimws(as.character(db3$Isomers_IDs))
    db3[
        !is.na(db3$Isomers_IDs) & nzchar(db3$Isomers_IDs), ,
        drop = FALSE
    ]
}


#' Regional Pathway Enrichment
#'
#' This the function used to compute the gene/metabolites set enrichment for
#' multi-omics spatial data
#'
#' @param object A Bioconductor experiment contains spatial
#'   metabolomics(SM)/transcriptomics(ST) data or both, if contains SM data, it
#'   should be annotated via SpaMTP::annotateSM function.
#' @param ident A character scalar specifying the region column in cell
#'   metadata.
#' @param deList A list consisting of differential expression data.frames for
#'   each input modality. Within each data.frame column names MUST include
#'   'cluster', 'gene', ('avg_log2FC' or 'logFC') and ('p_val_adj' or 'FDR').
#' @param analyteTypes Vector of character strings defining which analyte types
#'   to use. Options can be c("genes"), c("metabolites") or both (default =
#'   c("genes", "metabolites")).
#' @param smAssay Primary MSI experiment name (default = "main").
#' @param stAssay Paired transcriptome altExp name (default = "transcriptome").
#' @param smAssayName Expression assay name within the metabolomics experiment
#'   selected by `smAssay` (default = `"counts"`).
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param stAssayName Expression assay name within the transcriptomics
#'   experiment
#'   selected by `stAssay` (default = `"counts"`).
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param minPathSize The min number of metabolites in a specific pathway
#'   (default = 5).
#' @param maxPathSize The max number of metabolites in a specific pathway
#'   (default = 500).
#' @param pvalCutoffMets Adjusted p-value cutoff used when constructing
#'   metabolite ranks. Set to `1` to include annotated metabolites regardless
#'   of DE significance (default = 0.05).
#' @param pvalCutoffGenes A numerical value defining the adjusted p value
#'   cutoff for significant differentially expressed genes. If `NULL` cutoff =
#'   `0.05` (default = 0.05).
#' @param annotationScoreThreshold Minimum indexed annotation score used for
#'   pathway mapping. It can be changed without re-running `annotateSM()` when
#'   the object was annotated with `minScore = 0` (default = 0.05).
#' @param annotationSource Metabolite annotation provenance. The default,
#'   `"current"`, requires scored RaMP IDs from the indexed annotation
#'   pipeline. `"auto"` permits a warned fallback to a legacy stored `db_3`,
#'   while `"legacy"` explicitly requests that compatibility path.
#' @param verbose Boolean indicating whether to show informative messages. If
#'   FALSE these messages will be suppressed (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#'
#' @return A data frame of regional pathway enrichment results, not a modified
#'   experiment. The complete-ranks interface adds the number of measured and
#'   leading-edge genes and metabolites per pathway (`measured_genes`,
#'   `measured_metabolites`, `leading_genes`, `leading_metabolites`) and
#'   attaches the ranks, pathway coverage, mapping and multiple-testing
#'   provenance as attributes.
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' ranks <- matrix(c(3, 2, 1, -1, -2, -3, 0.5, -0.5),
#'     dimnames = list(paste0("GENE", 1:8), "left")
#' )
#' enrichment <- findRegionalPathways(spe,
#'     ident = "region", ranks = list(genes = ranks),
#'     pathwayIndex = index, minPathSize = 2, maxPathSize = 6, verbose = FALSE
#' )
#' enrichment[, c("pathwayName", "NES", "padj")]
#' @inheritParams fishersPathwayAnalysis
#' @details Gene mapping uses the same identity index as
#'   [fishersPathwayAnalysis()]. Multiple RaMP records of a gene have their
#'   memberships united. If multiple input features map to that gene with
#'   different differential statistics in one cluster, resolve the duplicate
#'   features before differential analysis; the function does not select an
#'   arbitrary statistic.
#' @param ranks Alternative to DE.list: a list named `genes` and/or
#'   `metabolites`, each either a feature-by-region matrix of complete,
#'   unfiltered ranks or a [findAllDEMs()] result (converted with
#'   [markerRanks()]). Metabolite rows must be explicit RaMP compound IDs. Map
#' expression with createPathwayAssay before calculating identity-level effects.
#' @param regions Optional subset of regions (rank columns) to test with the
#'   ranks interface; NULL tests all regions.
#' @param seed Random seed for fgsea's permutations. It is applied with
#'   [withr::local_seed()], so the caller's random number stream is unchanged.
#' @param rankScale RMS-scale each modality's rank vector before concatenating,
#' or none. Pathways can still be dominated by the more densely measured
#'   modality.
#' @param nPermSimple Preliminary permutations for the ranks interface.
#' @param BPPARAM BiocParallel execution parameters. The default uses one
#'   worker and shows progress only when `verbose = TRUE`. An explicit
#'   parameter object controls its own workers and progress display.
findRegionalPathways <- function(
    object,
    ident,
    deList = NULL,
    analyteTypes = c("genes", "metabolites"),
    smAssay = "main",
    stAssay = "transcriptome",
    smAssayName = "counts",
    stAssayName = "counts",
    minPathSize = 5,
    maxPathSize = 500,
    pvalCutoffMets = 0.05,
    pvalCutoffGenes = 0.05,
    annotationScoreThreshold = 0.05,
    annotationSource = c("current", "auto", "legacy"),
    verbose = TRUE,
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    geneMapping = c("auto", "hgnc", "ramp"),
    geneReference = NULL, geneIndex = NULL,
    geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL,
    organism = "Homo sapiens", pathwayIndex = NULL,
    ranks = NULL, rankScale = c("rms", "none"), nPermSimple = 1000,
    regions = NULL, seed = 1234,
    BPPARAM = BiocParallel::SerialParam(progressbar = verbose)
) {
    withr::local_seed(seed)
    if (!is.null(ranks)) {
        if (!is.null(deList)) stop("Supply ranks or deList, not both.", call. = FALSE)
        if (is.null(pathwayIndex)) {
            pathwayIndex <- buildPathwayIndex(
                database,
                match.arg(geneMapping), geneReference, geneIndex, geneReferenceVersion,
                geneReferenceLocalDir, organism, databaseVersion, match.arg(databaseSource), databaseLocalDir
            )
        }
        ranks <- lapply(
            ranks,
            function(r) if (is.list(r) && !is.null(r$DEMs)) markerRanks(r) else r
        )
        return(.regionalPathwaysFromRanks(
            ranks, pathwayIndex, minPathSize,
            maxPathSize, match.arg(rankScale), nPermSimple, regions, BPPARAM
        ))
    }
    database <- .pathway_database(database, pathwayIndex)
    annotationSource <- match.arg(annotationSource)
    database_resources <- .spamtp_db_bundle(
        c("chem_props", "source_df", "analytehaspathway", "pathway"),
        database = database,
        version = databaseVersion,
        source = match.arg(databaseSource),
        localDir = databaseLocalDir
    )
    gene_view <- .pathway_context(
        database_resources, database,
        if ("genes" %in% analyteTypes || !is.null(pathwayIndex)) geneMapping else "ramp",
        geneReference, geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism, pathwayIndex
    )
    database_resources <- gene_view$resources
    if (!is.null(gene_view$index)) .gene_check_experiment_species(object, stAssay, organism)
    chem_props <- database_resources$chem_props
    source_df <- database_resources$source_df
    analytehaspathway <- database_resources$analytehaspathway
    pathway <- database_resources$pathway

    ## Checks for ident in SpaMTP Object
    cellMetadata <- .cellMetadata(object)
    if (!(ident %in% colnames(cellMetadata))) {
        stop(
            "Ident: ",
            ident,
            " not found in the SpaMTP object's cell metadata. Make sure the column exists and is a factor."
        )
    }
    cluster_vector <- as.factor(cellMetadata[[ident]])
    cluster <- levels(cluster_vector)
    # Match each selected assay to the observation labels.
    if ("genes" %in% analyteTypes) {
        st_obj <- .assayData(object, stAssay, stAssayName)
        if (length(cluster_vector) != ncol(st_obj)) {
            stop(
                "ident must have one label per observation in the gene assay.",
                call. = FALSE
            )
        }
    }
    if ("metabolites" %in% analyteTypes) {
        sm_obj <- .assayData(object, smAssay, smAssayName)
        if (length(cluster_vector) != ncol(sm_obj)) {
            stop(
                "ident must have one label per observation in the metabolite assay.",
                call. = FALSE
            )
        }
    }
    annotation_metadata <- NULL
    db_3 <- NULL
    if ("metabolites" %in% analyteTypes) {
        # (2) Resolve annotations. Scored Ramp_IDs from the indexed pipeline are
        # the primary key; source-ID joins are only available in legacy mode.
        verbose_message(
            messageText = "Resolving current RaMP metabolite annotations ... ",
            verbose = verbose
        )
        db_3 <- .resolve_pathway_metabolite_annotations(
            .pathwayAnnotationObject(object, smAssay),
            annotationSource = annotationSource,
            scoreThreshold = annotationScoreThreshold,
            chemicalProperties = chem_props
        )
        annotation_metadata <- attr(db_3, "annotation_metadata")
        annotation_metadata$pathway_pval_cutoff_mets <- pvalCutoffMets
        if (isTRUE(verbose) && length(annotation_metadata$engine)) {
            ramp_label <- annotation_metadata$ramp_version
            if (is.null(ramp_label) || is.na(ramp_label)) ramp_label <- "unknown"
            message(
                "Annotation engine: ", annotation_metadata$engine,
                "; RaMP: ", ramp_label
            )
        }
    }

    ### Adding DE Results
    if (length(deList) != length(analyteTypes)) {
        stop(
            "Number of DE data.frames provided does not match the number of analyte types specified. Please make sure a DE dataframe is provided for each analyte type"
        )
    }
    verbose_message(
        messageText = "Constructing DE dataframes.... ",
        verbose = verbose
    )
    if (!is.null(names(deList))) {
        if (anyDuplicated(names(deList)) || !setequal(names(deList), analyteTypes)) {
            stop("Named deList entries must match analyteTypes exactly.", call. = FALSE)
        }
        deList <- deList[analyteTypes]
    }
    source_de <- deList
    deList <- stats::setNames(vector("list", length(analyteTypes)), analyteTypes)
    for (i in seq_along(analyteTypes)) {
        verbose_message(
            messageText = paste0(
                "Assuming DE.list[",
                i,
                "] contains ",
                analyteTypes[i],
                " results .... "
            ),
            verbose = verbose
        )
        DE <- .standardizeRegionalDE(source_de[[i]], i)
        if (analyteTypes[i] == "metabolites") {
            DE <- DE %>% dplyr::rename(mz_name = "gene")
            db_3 <- merge(db_3, DE, by = "mz_name")
            deList[[analyteTypes[i]]] <- db_3
        } else {
            if (!is.null(gene_view$index)) {
                source_gene <- .gene_prepare_de(DE, gene_view$index)
            } else {
                DE <- DE %>% dplyr::mutate(commonName = toupper(.data$gene))
                source_gene <- merge(
                    DE,
                    source_df[which(grepl(source_df$rampId, pattern = "RAMP_G")), ],
                    by = "commonName"
                )
            }
            deList[[analyteTypes[i]]] <- source_gene
        }
    }
    # Get pathway db
    verbose_message(
        messageText = "Constructing pathway database ...",
        verbose = verbose
    )
    measured_ids <- character()
    if ("genes" %in% analyteTypes) {
        measured_ids <- .pathway_gene_ids(
            rownames(.assayData(object, stAssay, stAssayName)),
            gene_view$pathway_index
        )
    }
    if ("metabolites" %in% analyteTypes) measured_ids <- unique(c(measured_ids, deList[["metabolites"]]$ramp_id))
    coverage_by_cluster <- list()

    gc()
    pb3 <- .progressBar(length(cluster), verbose)
    per_cluster <- lapply(cluster, function(i) {
        i <- as.character(i)
        ranks <- c()
        if ("metabolites" %in% analyteTypes) {
            ## metabolites
            DE_met <- deList[["metabolites"]]
            sub_db3 <- DE_met[which(as.character(DE_met$cluster) == i), ] %>%
                dplyr::filter(.data$p_val_adj <= pvalCutoffMets %||% 0.05) %>%
                dplyr::filter(!duplicated(.data$ramp_id))
            met_ranks <- scale(sub_db3$avg_log2FC, center = 0)
            names(met_ranks) <- sub_db3$ramp_id
            ranks <- c(ranks, met_ranks)
        }
        if ("genes" %in% names(deList)) {
            ## genes
            DE_rna <- deList[["genes"]]
            sub_de_gene <- DE_rna[which(as.character(DE_rna$cluster) == i), ] %>%
                dplyr::filter(.data$p_val_adj <= pvalCutoffGenes %||% 0.05) %>%
                dplyr::filter(!duplicated(.data$rampId))
            ranks_gene_vector <- scale(sub_de_gene$avg_log2FC, center = 0)
            names(ranks_gene_vector) <- sub_de_gene$rampId
            # Genes and metabolites
            ranks <- c(ranks, ranks_gene_vector)
        }

        ranks <- ranks[which(!duplicated(names(ranks)))]
        cluster_ranks <- ranks[is.finite(ranks)]

        coverage <- .pathway_coverage(gene_view$pathway_index,
            measured_ids,
            usedIds = names(cluster_ranks),
            types = ifelse(analyteTypes == "genes", "G", "C"),
            minSize = minPathSize, maxSize = maxPathSize
        )
        pathway_db <- .pathway_sets_from_coverage(coverage)
        gsea_result <- c()
        if (length(cluster_ranks) > 0) {
            suppressWarnings({
                gsea_result <- fgsea::fgsea(
                    pathways = pathway_db,
                    stats = cluster_ranks,
                    minSize = minPathSize,
                    maxSize = maxPathSize,
                    BPPARAM = BPPARAM
                ) %>% dplyr::mutate(Cluster_id = i)
            })
        } else {
            gsea_result <- data.table::data.table(
                pathway = character(0),
                pval = numeric(0),
                padj = numeric(0),
                log2err = numeric(0),
                ES = numeric(0),
                NES = numeric(0),
                size = integer(0),
                leadingEdge = list(),
                Cluster_id = i
            )
        }
        gsea_result <- na.omit(gsea_result) %>% dplyr::filter(!duplicated(.data$pathway))
        short_source <- source_df[which((source_df$rampId %in% names(cluster_ranks)) &
            !duplicated(source_df$rampId)), ]

        if (!nrow(gsea_result)) {
            return(list(coverage = coverage, result = NULL))
        }
        additional_entry <- do.call(rbind, lapply(
            gsea_result$leadingEdge, .regionalLeadingEdge,
            analyteTypes = analyteTypes,
            metaboliteDe = if ("metabolites" %in% analyteTypes) sub_db3 else NULL,
            geneSource = short_source, ranks = ranks,
            geneResults = "genes" %in% names(deList)
        ))
        gsea_result <- cbind(gsea_result, additional_entry)
        .progressTick(pb3, as.numeric(which(cluster == i)))
        list(coverage = coverage, result = gsea_result)
    })
    .progressDone(pb3)
    coverage_by_cluster[cluster] <- lapply(per_cluster, `[[`, "coverage")
    gsea_all_cluster <- do.call(
        rbind,
        c(list(data.frame()), lapply(per_cluster, `[[`, "result"))
    )

    if (!nrow(gsea_all_cluster)) {
        result <- data.frame(
            pathwayName = character(), pval = numeric(), padj = numeric(),
            NES = numeric(), Cluster_id = character()
        )
        attr(result, "pathway_coverage") <- coverage_by_cluster
        attr(result, "pathway_index") <- gene_view$pathway_index$provenance
        attr(result, "annotation_metadata") <- annotation_metadata
        attr(result, "gene_mapping") <- if (!is.null(gene_view$index)) attr(deList[["genes"]], "gene_mapping") else list(mode = "ramp")
        return(result)
    }
    gsea_all_cluster <- na.omit(gsea_all_cluster) %>%
        dplyr::mutate(group_importance = sum(abs(.data$NES)))
    colnames(gsea_all_cluster)[1] <- "pathwayRampId"
    gsea_all_cluster <- merge(
        gsea_all_cluster,
        gene_view$pathway_index$metadata,
        by = "pathwayRampId"
    )
    attr(gsea_all_cluster, "pathway_coverage") <- coverage_by_cluster
    attr(gsea_all_cluster, "pathway_index") <- gene_view$pathway_index$provenance
    attr(gsea_all_cluster, "annotation_metadata") <- annotation_metadata
    attr(gsea_all_cluster, "gene_mapping") <- if (!is.null(gene_view$index)) {
        attr(deList[["genes"]], "gene_mapping")
    } else {
        list(mode = "ramp")
    }
    return(gsea_all_cluster)
}


#' Runs multilevel Monte-Carlo variant for performing gene sets co-regulation
#' analysis using the RAMP_DB metabolite/gene database.
#'
#' This function is adapted from the [fgsea::geseca](https://github.com/alserglab/fgsea/blob/master/R/geseca-multilevel.R) package to identify significantly expressed RAMP_DB pathways based on an expression/feature embedding matrix.
#'
#' @param expression Expression matrix with gene identifiers or RaMP IDs as rows and
#'   observations as columns.
#' @param minSize Minimal size of a gene set to test. All pathways below the
#'   threshold are excluded (default = 1).
#' @param maxSize Maximal size of a gene set to test. All pathways above the
#'   threshold are excluded (default = `nrow(expression) - 1`).
#' @param center a logical value indicating whether the gene expression should
#'   be centered to have zero mean before the analysis takes place (default =
#'   TRUE).
#' @param scale a logical value indicating whether the gene expression should be
#'   scaled to have unit variance before the analysis takes place (default =
#'   FALSE).
#' @param sampleSize sample size for conditional sampling (default = 101).
#' @param eps This parameter sets the boundary for calculating P-values (default
#'   = 1e-50).
#' @param verbose Show progress with the default BPPARAM (default FALSE).
#' @param BPPARAM BiocParallel execution parameters. The default runs serially
#'   with progress controlled by `verbose`; explicit parameters control their
#'   own progress display.
#' @param nPermSimple Number of permutations in the simple geseca implementation
#'   for preliminary estimation of P-values (default = 1000).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#'
#' @inheritParams createPathwayAssay
#' @param organism Species of the input. HGNC mapping requires `"Homo sapiens"`.
#'   RaMP-only mapping (`geneMapping = "ramp"`, or a pathway index without an
#'   HGNC index) does not perform orthology conversion. Non-human gene analysis
#'   requires a separately curated species-specific database.
#' @return A table with GESECA results. Each row corresponds to a tested RAMP_DB
#'   pathway.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' expression <- SummarizedExperiment::assay(
#'     SingleCellExperiment::altExp(spe, "transcriptome"), "logcounts"
#' )
#' result <- withr::with_seed(1, runRAMPGeseca(expression,
#'     pathwayIndex = index,
#'     minSize = 2, maxSize = 6, nPermSimple = 100, eps = 0.01
#' ))
#' head(result)
runRAMPGeseca <- function(
    expression,
    minSize = 1,
    maxSize = nrow(expression) - 1,
    center = TRUE,
    scale = FALSE,
    sampleSize = 101,
    eps = 1e-50,
    BPPARAM = BiocParallel::SerialParam(progressbar = verbose),
    nPermSimple = 1000,
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL, pathwayIndex = NULL,
    geneMapping = c("auto", "hgnc", "ramp"), geneReference = NULL,
    geneIndex = NULL, geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL, organism = "Homo sapiens",
    duplicateGenes = c("error", "mean", "sum"),
    verbose = FALSE
) {
    database <- .pathway_database(database, pathwayIndex)
    if (is.null(pathwayIndex)) {
        pathwayIndex <- buildPathwayIndex(
            database,
            match.arg(geneMapping), geneReference, geneIndex, geneReferenceVersion,
            geneReferenceLocalDir, organism, databaseVersion, match.arg(databaseSource), databaseLocalDir
        )
    } else {
        .pathway_context(
            database, database, match.arg(geneMapping), geneReference,
            geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism, pathwayIndex
        )
    }
    mapped <- .pathway_expression(expression, pathwayIndex, match.arg(duplicateGenes))
    expression <- as.matrix(mapped$expression)
    if (missing(maxSize)) maxSize <- nrow(expression) - 1L
    types <- if (all(grepl("^RAMP_[GC]_", rownames(expression)))) {
        unique(sub("^RAMP_([GC])_.*", "\\1", rownames(expression)))
    } else {
        NULL
    }
    coverage <- .pathway_coverage(
        pathwayIndex,
        rownames(expression),
        types = types,
        minSize = minSize,
        maxSize = maxSize
    )
    pathway_db <- .pathway_sets_from_coverage(coverage)
    gesecaRes <- fgsea::geseca(pathway_db, expression,
        minSize = minSize, maxSize = maxSize,
        center = center, scale = scale, sampleSize = sampleSize, eps = eps,
        BPPARAM = BPPARAM, nPermSimple = nPermSimple
    )

    attr(gesecaRes, "pathway_coverage") <- coverage
    attr(gesecaRes, "pathway_index") <- pathwayIndex$provenance
    attr(gesecaRes, "gene_mapping") <- mapped$mapping
    return(gesecaRes)
}


#' Create a RaMP identity alternative experiment
#'
#' Maps genes or mass features to RaMP identities in a new `altExp()`.
#' Rows are gene/compound identities, not pathway scores. Use this identity
#' expression for enrichment or co-regulation analysis, or call
#' [createPathwayObject()] to construct pathway-level scores.
#'
#' @param object A SingleCellExperiment, including SpatialExperiment, containing
#'   metabolomic or transcriptomic expression.
#' @param analyteType Character string specifying the type of analytes to
#'   process.Must be either "genes" or "metabolites" (default = "metabolites").
#' @param assay Primary (`"main"`) or alternative experiment name. The default
#'   is `"main"` for the primary experiment.
#' @param assayName Expression assay name within the selected experiment
#'   (default = `"counts"`), read with `SummarizedExperiment::assay()`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param outputExperiment Name for the identity-level `altExp()` (default =
#'   `"pathway"`); this is not an expression assay name.
#' @param annotationScoreThreshold Minimum indexed annotation score used to
#'   map m/z features to RaMP compounds (default = 0.05).
#' @param annotationSource Metabolite annotation provenance. `"current"`
#'   requires the indexed, scored RaMP output; `"auto"` and `"legacy"` enable
#'   compatibility with older serialized SpaMTP objects.
#' @param verbose Boolean logical value indicating whether to print verbose
#'   messages during execution. (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#'
#' @inheritParams fishersPathwayAnalysis
#' @inheritParams buildGeneMappingIndex
#' @param duplicateGenes How to aggregate expression rows of one gene.
#'   The default errors for differing rows and deduplicates identical rows.
#'   Explicit mean or sum requires choosing an appropriate expression scale;
#'   sum is generally inappropriate for log-transformed expression.
#' @inheritParams annotateGeneIdentifiers
#' @details Gene expression uses the shared HGNC identity index. Non-conflicting
#'   RaMP records of one gene count once. Different expression rows for a gene
#'   require an explicit aggregation choice. In custom RaMP-only compatibility
#'   mode the historical default is mean. The selected expression assay name
#'   is preserved for genes, including logcounts. Inputs, member coverage,
#'   duplicate handling and index provenance are stored in
#'   metadata(altExp(x, outputExperiment))$pathway_mapping.
#' @return The input container with a new
#'   `SingleCellExperiment::altExp(SpaMTP, outputExperiment)`. Its expression assay
#'   retains the name supplied in `assayName`, and its `rowData()` describes the
#'   RaMP identities and original measured features.
#' @export
#'
#' @importFrom dplyr mutate group_by summarise ungroup
#' @importFrom tidyr separate_rows
#' @importFrom data.table as.data.table
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' mapped <- createPathwayAssay(spe,
#'     analyteType = "genes", assay = "transcriptome",
#'     outputExperiment = "gene_ids", pathwayIndex = index, verbose = FALSE
#' )
#' rownames(SingleCellExperiment::altExp(mapped, "gene_ids"))
#'
#' ## Create a pathway assay from gene data with verbose output
#' # spamtp_obj <- createPathwayAssay(spamtp_obj, analyte_type = "genes", assay = "SPT", outputExperiment = "gene_pathway", verbose = TRUE)
#' @param metaboliteAmbiguity Expand a mass feature into candidate compounds
#' (historical behavior), or exclude features mapping to more than one compound.
#'   The workflow explicitly defaults to exclude. Candidate pairs, excluded
#'   features and exact aggregation weights are retained in metadata.
createPathwayAssay <- function(
    object, analyteType = "metabolites", assay = "main", assayName = "counts", outputExperiment = "pathway", annotationScoreThreshold = 0.05, annotationSource = c("current", "auto", "legacy"), verbose = TRUE, database = NULL, databaseVersion = "latest", databaseSource = c("auto", "spamtpdb"), databaseLocalDir = NULL,
    geneMapping = c("auto", "hgnc", "ramp"), geneReference = NULL, geneIndex = NULL,
    geneReferenceVersion = "latest", geneReferenceLocalDir = NULL,
    organism = "Homo sapiens", pathwayIndex = NULL,
    duplicateGenes = c("error", "mean", "sum"), idColumn = NULL,
    metaboliteAmbiguity = c("expand", "exclude")
) {
    .requireExperiment(object, "SingleCellExperiment")
    geneMapping <- match.arg(geneMapping)
    database <- .pathway_database(database, pathwayIndex)
    if (identical(analyteType, "genes")) {
        needed <- "source_df"
        if (is.null(database) || "analytehaspathway" %in% names(database)) needed <- c(needed, "analytehaspathway", "pathway")
        resources <- .spamtp_db_bundle(needed, database,
            version = databaseVersion,
            source = match.arg(databaseSource), localDir = databaseLocalDir
        )
        view <- .pathway_context(
            resources, database, geneMapping, geneReference, geneIndex,
            geneReferenceVersion, geneReferenceLocalDir, organism, pathwayIndex
        )
        if (!is.null(view$index)) .gene_check_experiment_species(object, assay, organism)
        E <- .assayData(object, assay, assayName)
        ids <- rownames(E)
        if (!is.null(idColumn)) {
            rd <- .featureMetadata(object, assay)
            if (length(idColumn) != 1L || !idColumn %in% names(rd)) stop("idColumn was not found.", call. = FALSE)
            ids <- as.character(rd[[idColumn]])
        }
        policy <- if (is.null(view$index) && missing(duplicateGenes)) "mean" else match.arg(duplicateGenes)
        mapped <- .pathway_expression(
            E,
            view$pathway_index,
            policy,
            ids,
            verbose
        )
        coverage <- .pathway_coverage(
            view$pathway_index,
            rownames(mapped$expression),
            types = "G"
        )
        ids <- rownames(mapped$expression)
        node <- if (is.null(view$index)) NULL else view$index$nodes[match(ids, view$index$nodes$rampId), ]
        rd <- S4Vectors::DataFrame(rampId = ids, row.names = ids)
        if (!is.null(node)) {
            rd$gene_id <- node$gene_id
            rd$commonName <- node$symbol
        }
        rd$original_features <- unname(mapped$original_features)
        target <- SingleCellExperiment::SingleCellExperiment(
            assays = stats::setNames(list(mapped$expression), assayName), rowData = rd
        )
        S4Vectors::metadata(target)$pathway_mapping <- list(
            provenance = view$pathway_index$provenance, inputs = mapped$mapping,
            original_features = mapped$original_features, duplicate_genes = policy, coverage = coverage
        )
        SingleCellExperiment::altExp(object, outputExperiment) <- target
        return(object)
    }


    annotationSource <- match.arg(annotationSource)
    database_resources <- .spamtp_db_bundle(
        c("chem_props", "source_df"),
        database = database,
        version = databaseVersion,
        source = match.arg(databaseSource),
        localDir = databaseLocalDir
    )
    chem_props <- database_resources$chem_props
    source_df <- database_resources$source_df

    if (!analyteType %in% c("genes", "metabolites")) {
        stop("Incorrect `analyteType` provided! must be either 'genes' or 'metabolites'. Please provided the correct analyte matching the selected assay data.")
    }

    assayMatrix <- .assayData(object, assay, assayName)
    annotations <- .resolve_pathway_metabolite_annotations(
        .pathwayAnnotationObject(object, assay),
        annotationSource = annotationSource,
        scoreThreshold = annotationScoreThreshold, chemicalProperties = chem_props
    )
    annotation_pairs <- unique(annotations[c("mz_name", "ramp_id")])
    annotation_pairs <- annotation_pairs[annotation_pairs$mz_name %in% rownames(assayMatrix) &
        !is.na(annotation_pairs$ramp_id) & grepl("^RAMP_C_", annotation_pairs$ramp_id), , drop = FALSE]
    counts <- table(annotation_pairs$mz_name)
    conflicts <- names(counts)[counts > 1L]
    ambiguity <- match.arg(metaboliteAmbiguity)
    pairs <- annotation_pairs
    if (ambiguity == "exclude") pairs <- pairs[!pairs$mz_name %in% conflicts, , drop = FALSE]
    identifiers <- unique(pairs$ramp_id)
    if (!length(identifiers)) stop(
        "No features map to RaMP IDs after the annotation/ambiguity policy.",
        call. = FALSE
    )
    group <- match(pairs$ramp_id, identifiers)
    weights <- Matrix::sparseMatrix(
        i = group, j = match(pairs$mz_name, rownames(assayMatrix)),
        x = 1 / tabulate(group, length(identifiers))[group],
        dims = c(length(identifiers), nrow(assayMatrix)),
        dimnames = list(identifiers, rownames(assayMatrix))
    )
    pathwayMatrix <- weights %*% assayMatrix
    dimnames(pathwayMatrix) <- list(identifiers, colnames(assayMatrix))
    original_features <- stats::setNames(
        lapply(identifiers, function(id) pairs$mz_name[pairs$ramp_id == id]),
        identifiers
    )
    pathwayMetadata <- S4Vectors::DataFrame(
        rampId = identifiers,
        mz_name = vapply(
            original_features,
            paste,
            collapse = "; ",
            character(1)
        ),
        row.names = identifiers
    )
    pathwayExperiment <- SingleCellExperiment::SingleCellExperiment(
        assays = stats::setNames(list(pathwayMatrix), assayName), rowData = pathwayMetadata
    )
    S4Vectors::metadata(pathwayExperiment)$pathway_mapping <- list(
        inputs = annotation_pairs, ambiguity = ambiguity,
        excluded_conflicts = if (ambiguity == "exclude") conflicts else character(),
        candidate_conflicts = conflicts, aggregation = "Mean of measured features per explicit RaMP compound",
        original_features = original_features, weights = weights,
        annotation = attr(annotations, "annotation_metadata"),
        provenance = if (is.null(pathwayIndex)) list(database_version = databaseVersion) else pathwayIndex$provenance
    )
    SingleCellExperiment::altExp(object, outputExperiment) <- pathwayExperiment
    object
}


############################# Annotation Estimation using Pathway Results ########################################


#' Store pathway scores in an alternative experiment
#'
#' This function computes pathway-level scores from analyte-level expression
#' data
#' and stores them in altExp for SingleCellExperiment-derived objects.
#' Each score is the sum of matched expression divided by the square root
#' of the matched feature count, then centred and scaled across pixels.
#' Constant scores are zero; unmatched pathways are removed or stored as NA.
#'
#' @param object A SingleCellExperiment (including SpatialExperiment) with
#'   expression
#'   features indexed by gene identifiers or RaMP IDs.
#' @param assay Primary or alternative experiment to score.
#' @param assayName Expression assay. If omitted for a Bioconductor container,
#'   prefers logcounts, normcounts, then counts in the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param outputExperiment Name for the pathway-level `altExp()` (default =
#'   `"pathway"`).
#'   Scores are stored in that experiment's `pathwayScores` expression assay.
#' @param removeNans Remove pathways without matched analytes (default TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#'
#' @inheritParams createPathwayAssay
#' @inheritParams fishersPathwayAnalysis
#' @param minPathSize Minimum measured unique members to score (default 1).
#' @param maxPathSize Maximum measured unique members to score (default Inf).
#' @details Accepts gene identifiers directly or a harmonized RaMP assay.
#'   min_path_size and max_path_size are applied to measured unique members.
#'   rowData contains database_size (after identity reconciliation),
#'   database_raw_size, measured_size, coverage_fraction, used_size,
#' database_members, measured_members, used_members and
#'   excluded_conflict_records.
#'   The complete coverage table, including excluded and unmeasured pathways,
#'   is retained in metadata(altExp(x, outputExperiment))$pathway_mapping$coverage.
#'   Scores summarize observed expression; they are not pathway activity tests.
#' @return The input with pathwayScores in a new alternative experiment.
#'   Pathway identifiers are preserved.
#'
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildPathwayIndex(spamtp_example$database)
#' spe <- createPathwayObject(spe, assay = "transcriptome", pathwayIndex = index)
#' SummarizedExperiment::assay(SingleCellExperiment::altExp(spe, "pathway"), "pathwayScores")[, 1:4]
createPathwayObject <- function(
    object, assay = NULL, assayName = "logcounts",
    outputExperiment = "pathway", removeNans = TRUE, database = NULL,
    databaseVersion = "latest", databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL, geneMapping = c("auto", "hgnc", "ramp"),
    geneReference = NULL, geneIndex = NULL, geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL, organism = "Homo sapiens", pathwayIndex = NULL,
    duplicateGenes = c("error", "mean", "sum"), minPathSize = 1, maxPathSize = Inf
) {
    .requireExperiment(object, "SingleCellExperiment")
    if (missing(assayName)) assayName <- .integrationAssay(.experimentForAssay(object, assay))
    context <- .pathway_score_context(
        object, assay, assayName, database, pathwayIndex,
        match.arg(geneMapping), geneReference, geneIndex, geneReferenceVersion,
        geneReferenceLocalDir, organism, databaseVersion, match.arg(databaseSource),
        databaseLocalDir, match.arg(duplicateGenes), minPathSize, maxPathSize
    )
    coverage <- context$coverage
    keep <- coverage$eligible | (!removeNans & coverage$used_size == 0)
    selected <- coverage[keep, , drop = FALSE]
    if (!nrow(selected)) stop(
        "No pathways contain features from the selected assay within the size limits.",
        call. = FALSE
    )
    sets <- .pathway_sets_from_coverage(selected, eligibleOnly = FALSE)
    scores <- .pathway_score_matrix(context$expression, sets)
    target <- SingleCellExperiment::SingleCellExperiment(
        assays = list(pathwayScores = scores),
        rowData = S4Vectors::DataFrame(
            selected,
            row.names = selected$pathwayRampId
        )
    )
    S4Vectors::metadata(target)$pathway_mapping <- context$audit
    SingleCellExperiment::altExp(object, outputExperiment) <- target
    object
}


############################# PATHWAY HELPER FUNCTIONS ########################################

#' @noRd
get_analytes_db <- function(inputId, analytehaspathway, chemProps, pathway) {
    rampid <- unique(chemProps$ramp_id[which(chemProps$chem_source_id %in% unique(inputId))])
    #
    pathway_ids <- unique(analytehaspathway$pathwayRampId[which(analytehaspathway$rampId %in% rampid)])

    analytes_db <- lapply(pathway_ids, function(x) {
        content <- analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == x)]
        content <- content[which(grepl(content, pattern = "RAMP_C"))]
        return(content)
    })
    analytes_db_name <- unlist(lapply(pathway_ids, function(x) {
        name <- pathway$pathwayName[which(pathway$pathwayRampId == x)]
        return(name)
    }))
    names(analytes_db) <- analytes_db_name
    return(analytes_db)
}
