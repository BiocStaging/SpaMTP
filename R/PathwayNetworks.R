#' Construct an interactive pathway network visualization
#'
#' Builds a self-contained data payload and an interactive D3 visualization for
#' selected pathways across spatial clusters. Topology lookup, edge preparation,
#' differential-expression matching, and matrix extraction are indexed so they
#' are not repeated for every cluster-pathway pair.
#'
#' @param object A SpatialExperiment with paired transcriptomics in altExp.
#'   Subset to one sample before exporting. Metabolomics must be annotated by
#'   [annotateSM()], or supplied as explicit compound IDs from
#'   [createPathwayAssay()].
#' @param ident Metadata column used to identify spatial clusters or regions.
#' @param regpathway Data frame returned by [findRegionalPathways()].
#' @param deList One differential-expression data frame per requested analyte
#'   type. Data frames must contain `cluster`, `gene`, `avg_log2FC` (or
#'   `logFC`), and `p_val_adj` (or `FDR`). Descriptive native marker tables
#'   containing `mean_auc` and `mean_cohen` are also supported; their P values
#'   remain unavailable. A named list is recommended.
#' @param selectedPathways Optional pathway names or source IDs. Matching is
#'   case-insensitive. When `NULL`, the most important pathways are selected by
#'   summed absolute NES.
#' @param path Output directory for the generated HTML file. Defaults to the
#'   session's temporary directory; supply a persistent directory to retain it.
#' @param smAssay Primary (`"main"`) or alternative experiment containing MSI.
#' @param stAssay Paired transcriptome alternative experiment name.
#' @param smAssayName Expression assay name within the `smAssay` experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param stAssayName Expression assay name within the `stAssay` experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param colourPalette Colours used for the spatial abundance raster.
#' @param analyteTypes One or both of `"genes"` and `"metabolites"`.
#' @param annotationSource Metabolite annotation provenance. The default,
#'   `"current"`, requires the indexed, scored RaMP annotation output.
#'   `"auto"` permits a warned legacy fallback and `"legacy"` requests it
#'   explicitly.
#' @param annotationScoreThreshold Minimum indexed annotation score used for
#'   metabolite matching. When `NULL`, the threshold recorded by
#'   `findRegionalPathways()` is reused (default = `NULL`).
#'   An explicit compound-identity assay is already filtered; this argument
#'   does not filter it again and annotation-score controls are hidden.
#' @param annotationScoreFloor Lowest annotation score embedded in the HTML
#'   for interactive filtering. It is automatically lowered when the initial
#'   `annotationScoreThreshold` is smaller (default = `0.01`).
#'   Only applies to original mass-feature input with scored candidates.
#' @param metaboliteDetection Which metabolite nodes receive detected styling.
#'   `"leading_edge"` uses only GSEA leading-edge metabolites;
#'   `"annotated"` also includes retained compound identities, or score-filtered
#'   annotations for original mass-feature input, belonging to the selected
#'   pathway. These compounds need not be in the GSEA leading edge.
#' @param image Optional image ID validated against imgData. Coordinates always
#'   come from spatialCoords; the optical raster is not embedded in this viewer.
#' @param verbose Display progress messages.
#' @param topNPathways Number of pathways selected when `selectedPathways` is
#'   `NULL`.
#' @param maxNodes Target maximum number of reference nodes per
#'   cluster-pathway view. Leading-edge analytes and annotated compounds are
#'   retained, followed by their neighbours and high-degree nodes. Required
#'   members and isolated measured nodes can exceed this limit; other measured
#'   genes can be pruned. Use `Inf` to disable pruning.
#' @param labelMode Initial node-label mode: detected nodes only, all nodes, or
#'   no labels. It can also be changed interactively.
#' @param layoutMode Initial network layout. `"repulsion"` maximises spacing,
#'   `"force"` uses a balanced force-directed layout, `"radial"` separates
#'   analyte classes into rings, and `"bipartite"` separates genes and
#'   metabolites vertically. It can also be changed interactively.
#' @param maxSpatialPoints Maximum number of spatial points embedded in the
#'   HTML. Larger datasets are deterministically downsampled for responsive
#'   browser rendering. Use `Inf` to retain all points.
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @inheritParams fishersPathwayAnalysis
#' @param organism Species of the input. HGNC mapping requires `"Homo sapiens"`.
#'   RaMP-only mapping (`geneMapping = "ramp"`, or a pathway index without an
#'   HGNC index) does not perform orthology conversion. Non-human gene analysis
#'   requires a separately curated species-specific database; using human
#' topology for compound-only context does not establish conserved interactions.
#' @details Gene differential values use the shared HGNC mapping index. For
#'   display, the same resolved gene value is attached to its equivalent raw
#'   RaMP topology nodes; this does not introduce additional enrichment tests.
#'
#' @return Invisibly returns the generated HTML file path. Its
#'   `network_summary` attribute is a data frame with one row per non-empty
#'   region and pathway network (`region`, `pathway`, `nodes`, `edges`,
#'   `measured_compounds`); `pathway_coverage` records measured pathway
#'   coverage.
#' @seealso [embedPathwayNetwork()] to include the offline viewer in an HTML
#'   document with a standalone download.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' db <- spamtp_example$database
#' db$ramp_db_metadata <- list(ramp_version = "example")
#' db$pathway$sourceId <- paste0("WP", 1:3)
#' db$pathway$type <- "WikiPathways"
#' db$ramp_wikipathway <- list(list(
#'     id = "WP1", title = "Example pathway 1",
#'     protEdges = data.frame(
#'         src = "RAMP_G_1", dest = "RAMP_G_2",
#'         directed = 1L, reaction_type = 1L
#'     )
#' ))
#' regional <- data.frame(
#'     pathwayRampId = "P1", pathwayName = "Example pathway 1",
#'     sourceId = "WP1", type = "WikiPathways", Cluster_id = "left",
#'     NES = 1.5, padj = 0.03, leadingEdge = I(list(c("RAMP_G_1", "RAMP_G_2")))
#' )
#' de <- data.frame(
#'     gene = c("GENE1", "GENE2"), cluster = "left",
#'     avg_log2FC = c(1, 0.5), p_val_adj = c(0.01, 0.02)
#' )
#' folder <- tempfile("network-")
#' dir.create(folder)
#' output <- pathwayNetworkPlots(spe, "region", regional, list(genes = de),
#'     analyteTypes = "genes", database = db, path = folder, verbose = FALSE
#' )
#' basename(output)
#' unlink(folder, recursive = TRUE)
pathwayNetworkPlots <- function(
    object,
    ident,
    regpathway,
    deList,
    selectedPathways = NULL,
    path = tempdir(),
    smAssayName = "counts",
    stAssayName = "counts",
    colourPalette = NULL,
    smAssay = "main",
    stAssay = "transcriptome",
    analyteTypes = c("genes", "metabolites"),
    annotationSource = c("current", "auto", "legacy"),
    annotationScoreThreshold = NULL,
    annotationScoreFloor = 0.01,
    metaboliteDetection = c("leading_edge", "annotated"),
    image = NULL,
    verbose = TRUE,
    topNPathways = 10L,
    maxNodes = 500L,
    labelMode = c("detected", "all", "none"),
    maxSpatialPoints = 50000L,
    layoutMode = c("repulsion", "force", "radial", "bipartite"),
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    geneMapping = c("auto", "hgnc", "ramp"),
    geneReference = NULL, geneIndex = NULL,
    geneReferenceVersion = "latest",
    geneReferenceLocalDir = NULL,
    organism = "Homo sapiens", pathwayIndex = NULL
) {
    database <- .pathway_database(database, pathwayIndex)
    object <- .nativeSingleSample(object)
    analyteTypes <- match.arg(
        analyteTypes, c("genes", "metabolites"),
        several.ok = TRUE
    )
    analyteTypes <- unique(analyteTypes)
    annotationSource <- match.arg(annotationSource)
    metaboliteDetection <- match.arg(metaboliteDetection)
    labelMode <- match.arg(labelMode)
    layoutMode <- match.arg(layoutMode)
    databaseSource <- match.arg(databaseSource)
    database_resources <- .spamtp_db_bundle(
        c(
            "chem_props", "source_df", "analytehaspathway", "pathway",
            "ramp_db_metadata"
        ),
        database = database,
        version = databaseVersion,
        source = databaseSource,
        localDir = databaseLocalDir
    )
    chem_props <- database_resources$chem_props
    source_df <- database_resources$source_df
    analytehaspathway <- database_resources$analytehaspathway
    pathway <- database_resources$pathway
    if (!is.character(ident) || length(ident) != 1L || !nzchar(ident)) {
        stop("ident must be one metadata column name.")
    }
    if (!ident %in% names(.cellMetadata(object))) {
        stop(
            "Metadata column '",
            ident,
            "' was not found in SpaMTP cell metadata."
        )
    }
    if (!dir.exists(path)) stop("Output directory does not exist: ", path)
    if (file.access(path, mode = 2) != 0) {
        stop("Cannot write to output directory: ", path)
    }
    if (!is.numeric(maxNodes) || length(maxNodes) != 1L ||
        is.na(maxNodes) || maxNodes <= 0) {
        stop("maxNodes must be positive or Inf.")
    }
    if (!is.numeric(maxSpatialPoints) || length(maxSpatialPoints) != 1L ||
        is.na(maxSpatialPoints) || maxSpatialPoints <= 0) {
        stop("maxSpatialPoints must be positive or Inf.")
    }

    if ("genes" %in% analyteTypes) {
        .assayData(object, stAssay, stAssayName)
    }
    if ("metabolites" %in% analyteTypes) {
        .assayData(object, smAssay, smAssayName)
    }
    regpathway_annotation <- attr(regpathway, "annotation_metadata", exact = TRUE)
    if (is.null(annotationScoreThreshold)) {
        annotationScoreThreshold <- regpathway_annotation$analysis_score_threshold %||% 0.05
    }
    if (!is.numeric(annotationScoreThreshold) ||
        length(annotationScoreThreshold) != 1L ||
        !is.finite(annotationScoreThreshold) || annotationScoreThreshold < 0) {
        stop("annotationScoreThreshold must be one non-negative finite number.")
    }
    if (!is.numeric(annotationScoreFloor) ||
        length(annotationScoreFloor) != 1L ||
        !is.finite(annotationScoreFloor) || annotationScoreFloor < 0) {
        stop("annotationScoreFloor must be one non-negative finite number.")
    }
    annotationScoreFloor <- min(
        annotationScoreFloor, annotationScoreThreshold
    )
    compound_ids <- if ("metabolites" %in% analyteTypes) {
        rownames(.assayData(object, smAssay, smAssayName))
    } else {
        character()
    }
    identity_input <- length(compound_ids) > 0L && all(grepl("^RAMP_C_", compound_ids))
    if (identity_input) {
        identity_mapping <- S4Vectors::metadata(.experimentForAssay(object, smAssay))$pathway_mapping
        rank_index <- attr(regpathway, "pathway_index", exact = TRUE)
        if (!is.null(identity_mapping$provenance) && !is.null(rank_index) &&
            !identical(identity_mapping$provenance, rank_index)) {
            stop(
                "Compound identity assay and regional enrichment use different pathway indexes.",
                call. = FALSE
            )
        }
    }
    if ("metabolites" %in% analyteTypes && annotationSource == "current" && !identity_input) {
        if (!.pn_is_current_annotation_metadata(regpathway_annotation)) {
            stop(
                "regpathway was not generated with a current scored annotation ",
                "pipeline. Re-run findRegionalPathways() before plotting, or ",
                "explicitly select annotationSource = 'legacy'."
            )
        }
        current_ramp <- as.character(
            database_resources$ramp_db_metadata$ramp_version %||%
                .annotation_ramp_version()
        )
        if (!identical(regpathway_annotation$ramp_version, current_ramp)) {
            stop(
                "regpathway uses RaMP ", regpathway_annotation$ramp_version %||% "unknown",
                " but this SpaMTP build contains RaMP ", current_ramp,
                ". Re-run findRegionalPathways()."
            )
        }
    }

    verbose_message("Selecting enriched pathways ... ", verbose)
    pathway_rows <- .pn_select_pathway_rows(
        as.data.frame(regpathway), selectedPathways, topNPathways
    )
    requested_sources <- unique(vapply(pathway_rows$type, .pn_type_key, character(1)))
    if (anyNA(requested_sources)) {
        requested_sources <- c("wiki", "reactome", "kegg", "hmdb")
    }
    topology_resource_names <- c(
        wiki = "ramp_wikipathway",
        reactome = "ramp_reactome",
        kegg = "ramp_kegg",
        hmdb = "ramp_hmdb"
    )[requested_sources]
    topology_resources <- .spamtp_db_bundle(
        unname(topology_resource_names),
        database = database,
        version = databaseVersion,
        source = databaseSource,
        localDir = databaseLocalDir
    )
    topology_sources <- list()
    for (source_name in names(topology_resource_names)) {
        topology_sources[[source_name]] <- topology_resources[[topology_resource_names[[source_name]]]]
    }
    catalog <- .pn_topology_catalog(topology_sources)
    resolved <- .pn_resolve_topologies(pathway_rows, catalog, topology_sources)
    if (any(!resolved$found)) {
        missing <- resolved$pathways$pathwayName[!resolved$found]
        verbose_message(
            paste0(
                "No stored topology was found for: ",
                paste(missing, collapse = "; ")
            ),
            verbose
        )
    }
    topologies <- resolved$topologies[resolved$found]
    if (!length(topologies)) {
        stop("None of the selected pathways has a stored topological structure.")
    }
    pathway_names <- names(topologies)
    pathway_rows <- pathway_rows[
        tolower(pathway_rows$pathwayName) %in% tolower(pathway_names), ,
        drop = FALSE
    ]

    verbose_message("Preparing pathway edges once per topology ... ", verbose)
    styles <- new.env(parent = emptyenv())
    utils::data("reaction_type", package = "SpaMTP", envir = styles)
    prepared_edges <- lapply(
        topologies,
        .pn_prepare_edges,
        reactionTypes = styles$reaction_type
    )
    usable <- lengths(lapply(prepared_edges, function(x) unique(c(x$src, x$dest)))) > 0L
    if (any(!usable)) {
        verbose_message(
            paste0(
                "No valid edges were found for: ",
                paste(names(prepared_edges)[!usable], collapse = "; ")
            ),
            verbose
        )
    }
    prepared_edges <- prepared_edges[usable]
    pathway_names <- names(prepared_edges)
    if (!length(pathway_names)) stop("Selected topologies contain no valid network edges.")
    pathway_rows <- pathway_rows[
        tolower(pathway_rows$pathwayName) %in% tolower(pathway_names), ,
        drop = FALSE
    ]

    verbose_message("Indexing differential results ... ", verbose)
    differential <- .pn_normalise_de_list(deList, analyteTypes)
    gene_de <- NULL
    metabolite_de <- NULL
    metabolite_annotations <- NULL
    annotation_metadata <- NULL
    gene_view <- .pathway_context(
        database_resources, database,
        if ("genes" %in% analyteTypes || !is.null(pathwayIndex)) geneMapping else "ramp",
        geneReference, geneIndex, geneReferenceVersion, geneReferenceLocalDir, organism, pathwayIndex
    )
    if ("genes" %in% analyteTypes) {
        if (!is.null(gene_view$index)) .gene_check_experiment_species(object, stAssay, organism)
        gene_de <- .pn_prepare_gene_de(
            differential[["genes"]],
            source_df,
            gene_view$index
        )
    }
    if ("metabolites" %in% analyteTypes) {
        ids <- rownames(.assayData(object, smAssay, smAssayName))
        if (length(ids) && all(grepl("^RAMP_C_", ids))) {
            # A createPathwayAssay result already has explicit identities. Reusing the
            # original mass-candidate expansion here would change the member universe.
            db_3 <- data.frame(
                ramp_id = ids, mz_name = ids, Adduct = "",
                IsomerNames = NA_character_, observed_mz = NA_real_, chem_source_id = ""
            )
            attr(db_3, "annotation_metadata") <- list(
                engine = "explicit-compound-identities",
                pathway_index = identity_mapping$provenance, source_annotation = identity_mapping$annotation,
                ambiguity = identity_mapping$ambiguity, excluded_conflicts = identity_mapping$excluded_conflicts,
                original_features = identity_mapping$original_features
            )
        } else {
            db_3 <- .resolve_pathway_metabolite_annotations(
                .pathwayAnnotationObject(object, smAssay),
                annotationSource = annotationSource,
                scoreThreshold = annotationScoreFloor,
                chemicalProperties = chem_props
            )
        }
        metabolite_annotations <- db_3
        annotation_metadata <- attr(db_3, "annotation_metadata")
        metabolite_de <- .pn_prepare_metabolite_de(
            differential[["metabolites"]], db_3, chem_props
        )
    }
    gene_lookup <- .pn_best_de_by_cluster(gene_de, "rampId")
    metabolite_lookup <- .pn_best_de_by_cluster(metabolite_de, "ramp_id")

    all_node_ids <- unique(unlist(lapply(prepared_edges, function(x) {
        c(x$src, x$dest)
    }), use.names = FALSE))
    annotated_metabolites <- stats::setNames(
        vector("list", length(pathway_names)), pathway_names
    )
    annotation_scores <- numeric()
    if ("metabolites" %in% analyteTypes) {
        annotated_metabolites <- .pn_annotated_metabolites_by_pathway(
            metabolite_annotations,
            pathway_names,
            pathway_rows,
            analytehaspathway,
            pathway
        )
        annotation_scores <- .pn_annotation_score_map(metabolite_annotations)
    }
    all_node_ids <- unique(c(
        all_node_ids,
        unlist(annotated_metabolites, use.names = FALSE)
    ))
    node_names <- .pn_name_map(source_df, all_node_ids, chem_props)
    clusters <- stringr::str_sort(
        unique(as.character(pathway_rows$Cluster_id)),
        numeric = TRUE
    )
    built <- .pn_build_network_collection(
        pathway_rows, pathway_names, prepared_edges, clusters,
        gene_lookup, metabolite_lookup, node_names, maxNodes,
        metaboliteDetection = metaboliteDetection,
        annotatedMetabolites = annotated_metabolites,
        annotationScores = annotation_scores,
        annotationScoreThreshold = annotationScoreThreshold,
        fixedCompoundIds = identity_input, geneIndex = gene_view$index
    )

    verbose_message("Extracting only plotted assay features ... ", verbose)
    gene_matrix <- if ("genes" %in% analyteTypes) {
        .pn_layer_data(object, stAssay, stAssayName)
    } else {
        NULL
    }
    metabolite_matrix <- if ("metabolites" %in% analyteTypes) {
        .pn_layer_data(object, smAssay, smAssayName)
    } else {
        NULL
    }
    spatial <- .pn_spatial_payload(
        object, ident, image, built$matrix_ids,
        gene_de, metabolite_de, gene_matrix, metabolite_matrix,
        maxSpatialPoints
    )

    if (is.null(colourPalette)) {
        colourPalette <- grDevices::colorRampPalette(
            c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921")
        )(100)
    }
    if (!is.character(colourPalette) || length(colourPalette) < 2L) {
        stop("colourPalette must contain at least two colours.")
    }

    measured_ids <- used_ids <- character()
    if ("genes" %in% analyteTypes) {
        measured_ids <- .pathway_gene_ids(
            rownames(gene_matrix),
            gene_view$pathway_index
        )
        used_ids <- .pathway_gene_ids(
            unique(gene_de$rampId),
            gene_view$pathway_index
        )
    }
    if ("metabolites" %in% analyteTypes) {
        measured_ids <- unique(c(
            measured_ids,
            unlist(annotated_metabolites, use.names = FALSE)
        ))
        used_ids <- unique(c(used_ids, metabolite_de$ramp_id))
    }
    pathway_coverage <- .pathway_coverage(gene_view$pathway_index, measured_ids, used_ids,
        types = ifelse(analyteTypes == "genes", "G", "C")
    )
    selected_ids <- pathway_rows$pathwayRampId
    if (is.null(selected_ids)) {
        selected_ids <- pathway_coverage$pathwayRampId[
            tolower(pathway_coverage$pathwayName) %in% tolower(pathway_rows$pathwayName)
        ]
    }
    pathway_coverage <- pathway_coverage[pathway_coverage$pathwayRampId %in% selected_ids, , drop = FALSE]
    payload <- list(
        title = paste("SpaMTP pathway networks -", ident),
        pathways = pathway_names,
        clusters = clusters,
        networks = built$networks,
        spatial = spatial,
        palette = colourPalette,
        fc_limit = .pn_fc_limit(differential),
        label_mode = labelMode,
        layout_mode = layoutMode,
        metadata = list(
            max_nodes = maxNodes,
            effect_label = if (any(vapply(differential, function(d) {
                all(c("mean_auc", "mean_cohen") %in% names(d))
            }, logical(1)))) {
                "Regional effect (workflow units)"
            } else {
                "Log2 fold change"
            },
            annotation = annotation_metadata,
            gene_mapping = attr(gene_de, "gene_mapping"),
            pathway_index = gene_view$pathway_index$provenance,
            pathway_coverage = pathway_coverage,
            annotation_score_threshold = annotationScoreThreshold,
            annotation_score_floor = annotationScoreFloor,
            annotation_score_ceiling = max(
                c(annotationScoreThreshold, annotation_scores),
                na.rm = TRUE
            ),
            metabolite_detection = metaboliteDetection,
            annotation_controls = "metabolites" %in% analyteTypes && !identity_input,
            generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
        )
    )
    html <- .pn_render_html(payload)

    file_stem <- gsub("[^A-Za-z0-9._-]+", "_", ident)
    return_name <- paste0(
        file_stem,
        "_",
        format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
        ".html"
    )
    full_path <- file.path(path, return_name)
    writeLines(html, full_path, useBytes = TRUE)
    verbose_message(paste0("Pathway network written to ", full_path), verbose)
    attr(full_path, "network_summary") <- .pn_network_summary(
        built$networks, clusters, pathway_names
    )
    attr(full_path, "pathway_coverage") <- pathway_coverage
    attr(full_path, "pathway_index") <- gene_view$pathway_index$provenance
    invisible(full_path)
}

# One row per non-empty region x pathway network drawn in the HTML view.
.pn_network_summary <- function(networks, clusters, pathways) {
    grid <- expand.grid(
        pathway = seq_along(pathways),
        region = seq_along(clusters)
    )
    rows <- Map(function(i, j) {
        network <- networks[[i]][[j]]
        data.frame(
            region = clusters[[i]], pathway = pathways[[j]],
            nodes = nrow(network$nodes), edges = nrow(network$links),
            measured_compounds = sum(network$nodes$annotated)
        )
    }, grid$region, grid$pathway)
    summary <- do.call(rbind, rows)
    summary <- summary[summary$nodes > 0L, , drop = FALSE]
    rownames(summary) <- NULL
    summary
}
