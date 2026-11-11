.spamtp_pathway_cache <- new.env(parent = emptyenv())

.pn_is_current_annotation_metadata <- function(metadata) {
    current_engines <- c(.spamtp_annotation_engine, "curated-fmp10-v2")
    is.list(metadata) && length(metadata$engine) == 1L &&
        metadata$engine %in% current_engines
}

.pn_assert_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) {
        stop(
            label,
            " is missing required column(s): ",
            paste(missing, collapse = ", ")
        )
    }
    invisible(x)
}

.pn_type_key <- function(x) {
    value <- tolower(trimws(as.character(x %||% "")))
    if (grepl("wiki", value, fixed = TRUE)) {
        return("wiki")
    }
    if (grepl("reactome", value, fixed = TRUE)) {
        return("reactome")
    }
    if (grepl("kegg", value, fixed = TRUE)) {
        return("kegg")
    }
    if (grepl("hmdb", value, fixed = TRUE)) {
        return("hmdb")
    }
    NA_character_
}

.pn_topology_signature <- function(topologySources) {
    paste(
        names(topologySources),
        lengths(topologySources),
        vapply(topologySources, function(x) {
            if (!length(x)) {
                return("")
            }
            paste(x[[1]]$id %||% "", x[[length(x)]]$id %||% "", sep = ":")
        }, character(1)),
        collapse = "|"
    )
}

.pn_topology_catalog <- function(topologySources, useCache = TRUE) {
    if (is.null(names(topologySources)) || any(!nzchar(names(topologySources)))) {
        stop("topologySources must be a named list.")
    }
    signature <- .pn_topology_signature(topologySources)
    if (isTRUE(useCache) &&
        identical(.spamtp_pathway_cache$signature, signature) &&
        is.data.frame(.spamtp_pathway_cache$catalog)) {
        return(.spamtp_pathway_cache$catalog)
    }

    parts <- Map(function(collection, source) {
        if (!length(collection)) {
            return(NULL)
        }
        titles <- names(collection)
        if (is.null(titles)) {
            titles <- vapply(
                collection,
                function(x) as.character(x$title %||% ""),
                character(1)
            )
        }
        ids <- vapply(
            collection,
            function(x) as.character(x$id %||% ""),
            character(1)
        )
        data.frame(
            source = source,
            index = seq_along(collection),
            name_key = tolower(trimws(titles)),
            id_key = tolower(trimws(ids)),
            stringsAsFactors = FALSE
        )
    }, topologySources, names(topologySources))
    catalog <- do.call(rbind, parts)
    rownames(catalog) <- NULL
    if (isTRUE(useCache)) {
        .spamtp_pathway_cache$signature <- signature
        .spamtp_pathway_cache$catalog <- catalog
    }
    catalog
}

.pn_select_pathway_rows <- function(
    regpathway, selectedPathways = NULL,
    topNPathways = 10L
) {
    .pn_assert_columns(
        regpathway,
        c("pathwayName", "sourceId", "type", "Cluster_id", "leadingEdge"),
        "regpathway"
    )
    if (!nrow(regpathway)) stop("regpathway contains no rows.")

    name_key <- tolower(trimws(as.character(regpathway$pathwayName)))
    id_key <- tolower(trimws(as.character(regpathway$sourceId)))
    if (is.null(selectedPathways)) {
        .pn_assert_columns(regpathway, "NES", "regpathway")
        score <- suppressWarnings(abs(as.numeric(regpathway$NES)))
        score[!is.finite(score)] <- 0
        totals <- rowsum(score, group = name_key, reorder = FALSE, na.rm = TRUE)
        ranked <- rownames(totals)[order(totals[, 1], decreasing = TRUE)]
        selected_key <- utils::head(ranked, max(1L, as.integer(topNPathways)))
        keep <- name_key %in% selected_key
    } else {
        selected_key <- tolower(trimws(as.character(selectedPathways)))
        selected_key <- selected_key[nzchar(selected_key)]
        if (!length(selected_key)) stop("selectedPathways contains no usable values.")
        keep <- name_key %in% selected_key | id_key %in% selected_key
    }
    result <- regpathway[keep, , drop = FALSE]
    if (!nrow(result)) {
        stop("None of selectedPathways matched pathwayName or sourceId in regpathway.")
    }
    result
}

.pn_resolve_topologies <- function(pathwayRows, catalog, topologySources) {
    unique_rows <- pathwayRows[!duplicated(tolower(pathwayRows$pathwayName)), , drop = FALSE]
    resolved <- vector("list", nrow(unique_rows))
    found <- logical(nrow(unique_rows))
    for (i in seq_len(nrow(unique_rows))) {
        source <- .pn_type_key(unique_rows$type[[i]])
        candidates <- catalog
        if (!is.na(source)) {
            candidates <- candidates[candidates$source == source, , drop = FALSE]
        }
        name_key <- tolower(trimws(as.character(unique_rows$pathwayName[[i]])))
        id_key <- tolower(trimws(as.character(unique_rows$sourceId[[i]])))
        hit <- which(candidates$name_key == name_key | candidates$id_key == id_key)
        if (length(hit)) {
            row <- candidates[hit[[1]], , drop = FALSE]
            resolved[[i]] <- topologySources[[row$source]][[row$index]]
            found[[i]] <- TRUE
        }
    }
    names(resolved) <- as.character(unique_rows$pathwayName)
    list(topologies = resolved, found = found, pathways = unique_rows)
}

.pn_prepare_edges <- function(topology, reactionTypes) {
    edge_sets <- c(
        "mixedEdges", "metabolEdges", "metabolPropEdges",
        "protEdges", "protPropEdges"
    )
    parts <- lapply(intersect(edge_sets, names(topology)), function(name) {
        edge <- topology[[name]]
        if (!is.data.frame(edge) || !nrow(edge) ||
            !all(c("src", "dest", "reaction_type") %in% names(edge))) {
            return(NULL)
        }
        edge <- edge[, intersect(c(
            "src", "dest", "directed", "reaction_type",
            "source_reaction_type"
        ), names(edge)), drop = FALSE]
        if (!"directed" %in% names(edge)) edge$directed <- NA_integer_
        if (!"source_reaction_type" %in% names(edge)) edge$source_reaction_type <- NA_character_
        edge
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) {
        return(.pn_empty_links())
    }

    edges <- unique(do.call(rbind, parts))
    edges$src <- as.character(edges$src)
    edges$dest <- as.character(edges$dest)
    keep <- !is.na(edges$src) & nzchar(edges$src) &
        !is.na(edges$dest) & nzchar(edges$dest)
    edges <- edges[keep, , drop = FALSE]
    if (!nrow(edges)) {
        return(.pn_empty_links())
    }

    reaction_index <- match(
        as.character(edges$reaction_type),
        as.character(reactionTypes$reaction_type)
    )
    lookup <- function(column, default) {
        value <- as.character(reactionTypes[[column]][reaction_index])
        value[is.na(value) | !nzchar(value)] <- default
        value
    }
    edges$reaction_name <- lookup("reaction_name", "Interaction")
    labelled <- !is.na(edges$source_reaction_type) & nzchar(edges$source_reaction_type)
    edges$reaction_name[labelled] <- as.character(edges$source_reaction_type[labelled])
    edges$linetype <- lookup("linetype", "solid")
    edges$arrowhead <- lookup("arrowhead", "arrow")
    edges$colour <- lookup("colour", "#64748b")
    # graphite direction code 2 is undirected; do not imply a causal arrow.
    edges$arrowhead[!is.na(edges$directed) & edges$directed == 2L] <- "none"

    degree <- table(c(edges$src, edges$dest))
    edges$weight <- as.integer(degree[edges$src]) + as.integer(degree[edges$dest])
    rownames(edges) <- NULL
    edges
}

.pn_empty_nodes <- function() {
    data.frame(
        id = character(), group = character(), expr = numeric(),
        name = character(), display = character(), detected = logical(),
        leading_edge = logical(), annotated = logical(),
        annotation_score = numeric(),
        stringsAsFactors = FALSE
    )
}

.pn_empty_links <- function() {
    data.frame(
        source = character(), target = character(), colour = character(),
        weight = integer(), style = character(), head = character(),
        reaction = character(), stringsAsFactors = FALSE
    )
}

.pn_prune_edges <- function(edges, detected, maxNodes) {
    if (!nrow(edges)) {
        return(edges)
    }
    nodes <- unique(c(edges$src, edges$dest))
    if (!is.finite(maxNodes)) {
        return(edges)
    }
    maxNodes <- as.integer(maxNodes)
    if (maxNodes <= 0L || length(nodes) <= maxNodes) {
        return(edges)
    }
    detected <- intersect(unique(as.character(detected)), nodes)
    incident <- edges$src %in% detected | edges$dest %in% detected
    neighbours <- unique(c(edges$src[incident], edges$dest[incident]))
    degree <- sort(table(c(edges$src, edges$dest)), decreasing = TRUE)
    cap <- max(maxNodes, length(detected))
    keep_nodes <- unique(c(detected, neighbours, names(degree)))
    keep_nodes <- utils::head(keep_nodes, cap)
    edges[edges$src %in% keep_nodes & edges$dest %in% keep_nodes, , drop = FALSE]
}

.pn_expand_values <- function(df, column, separator = ";") {
    values <- strsplit(as.character(df[[column]]), separator, fixed = TRUE)
    lengths <- lengths(values)
    keep <- lengths > 0L
    if (!any(keep)) {
        return(df[0, , drop = FALSE])
    }
    expanded <- df[rep(which(keep), lengths[keep]), , drop = FALSE]
    expanded[[column]] <- trimws(unlist(values[keep], use.names = FALSE))
    expanded[nzchar(expanded[[column]]) & !is.na(expanded[[column]]), , drop = FALSE]
}

.pn_normalise_source_id <- function(x) {
    x <- trimws(as.character(x))
    has_prefix <- grepl(":", x, fixed = TRUE)
    source <- ifelse(has_prefix, tolower(sub(":.*$", "", x)), "")
    value <- ifelse(has_prefix, sub("^[^:]+:", "", x), x)
    is_hmdb <- grepl("^HMDB", value, ignore.case = TRUE) | source == "hmdb"
    is_chebi <- grepl("^CHEBI", value, ignore.case = TRUE) | source == "chebi"
    is_lipid <- grepl("^(LM|LMPK)", value, ignore.case = TRUE) |
        source %in% c("lipidmaps", "lipidmap")
    value[is_chebi] <- sub("^CHEBI:", "", value[is_chebi], ignore.case = TRUE)
    result <- x
    result[is_hmdb] <- paste0("hmdb:", toupper(value[is_hmdb]))
    result[is_chebi] <- paste0("chebi:", value[is_chebi])
    result[is_lipid] <- paste0("lipidmaps:", value[is_lipid])
    result
}

.pn_normalise_de_list <- function(deList, analyteTypes) {
    if (!is.list(deList) || length(deList) != length(analyteTypes)) {
        stop("deList must contain exactly one data.frame for each analyte type.")
    }
    if (all(analyteTypes %in% names(deList))) {
        result <- deList[analyteTypes]
    } else {
        result <- stats::setNames(deList, analyteTypes)
    }
    for (type in analyteTypes) {
        de <- as.data.frame(result[[type]], stringsAsFactors = FALSE)
        if ("logFC" %in% names(de) && !"avg_log2FC" %in% names(de)) {
            names(de)[names(de) == "logFC"] <- "avg_log2FC"
        }
        if ("FDR" %in% names(de) && !"p_val_adj" %in% names(de)) {
            names(de)[names(de) == "FDR"] <- "p_val_adj"
        }
        if (!"p_val_adj" %in% names(de) && all(c("mean_auc", "mean_cohen") %in% names(de))) {
            # Native descriptive markers have effects, but no pixel P values.
            de$p_val_adj <- NA_real_
        }
        .pn_assert_columns(
            de, c("cluster", "gene", "avg_log2FC", "p_val_adj"),
            paste0("deList[['", type, "']]")
        )
        de$cluster <- as.character(de$cluster)
        de$gene <- as.character(de$gene)
        de$avg_log2FC <- suppressWarnings(as.numeric(de$avg_log2FC))
        de$p_val_adj <- suppressWarnings(as.numeric(de$p_val_adj))
        de$p_val_adj[!is.finite(de$p_val_adj)] <- Inf
        de$.cluster_key <- tolower(de$cluster)
        result[[type]] <- de
    }
    result
}

.pn_prepare_gene_de <- function(de, sourceTable, geneIndex = NULL) {
    if (!is.null(geneIndex)) {
        return(.gene_prepare_de(de, geneIndex, expandRecords = TRUE))
    }
    direct <- startsWith(toupper(de$gene), "RAMP_G_")
    direct_de <- de[direct, , drop = FALSE]
    direct_de$rampId <- toupper(direct_de$gene)

    symbol_de <- de[!direct, , drop = FALSE]
    mapped <- symbol_de[0, , drop = FALSE]
    mapped$rampId <- character()
    if (nrow(symbol_de)) {
        symbols <- unique(toupper(symbol_de$gene))
        gene_source <- sourceTable[
            startsWith(sourceTable$rampId, "RAMP_G_") &
                toupper(sourceTable$commonName) %in% symbols,
            c("rampId", "commonName"),
            drop = FALSE
        ]
        gene_source$gene_key <- toupper(gene_source$commonName)
        gene_source <- unique(gene_source[c("rampId", "gene_key")])
        symbol_de$gene_key <- toupper(symbol_de$gene)
        mapped <- merge(symbol_de, gene_source, by = "gene_key", sort = FALSE)
        mapped$gene_key <- NULL
    }
    result <- rbind(direct_de, mapped[, names(direct_de), drop = FALSE])
    if (!nrow(result)) {
        return(result)
    }
    result$commonName <- toupper(result$gene)
    result
}

.pn_prepare_metabolite_map <- function(db3, chemicalProperties = NULL) {
    db3 <- as.data.frame(db3, stringsAsFactors = FALSE)
    if (!"mz_name" %in% names(db3)) {
        .pn_assert_columns(db3, "observed_mz", "stored SpaMTP db_3 result")
        db3$mz_name <- paste0("mz-", db3$observed_mz)
    }
    adduct <- if ("Adduct" %in% names(db3)) as.character(db3$Adduct) else ""
    db3$Adduct <- adduct

    if ("ramp_id" %in% names(db3) &&
        any(grepl("^RAMP_C_", as.character(db3$ramp_id)))) {
        mapping <- db3
        mapping$ramp_id <- toupper(as.character(mapping$ramp_id))
        mapping$common_name <- if ("IsomerNames" %in% names(mapping)) {
            as.character(mapping$IsomerNames)
        } else if ("common_name" %in% names(mapping)) {
            as.character(mapping$common_name)
        } else {
            NA_character_
        }
    } else if ("Ramp_IDs" %in% names(db3) && any(nzchar(stats::na.omit(db3$Ramp_IDs)))) {
        mapping <- .pn_expand_values(db3, "Ramp_IDs")
        mapping$ramp_id <- mapping$Ramp_IDs
        mapping$common_name <- if ("IsomerNames" %in% names(mapping)) {
            as.character(mapping$IsomerNames)
        } else {
            NA_character_
        }
    } else {
        if (is.null(chemicalProperties)) {
            stop("Legacy metabolite annotations require the chemProps lookup table.")
        }
        id_column <- if ("Isomers_IDs" %in% names(db3)) "Isomers_IDs" else "Isomers"
        .pn_assert_columns(db3, id_column, "stored SpaMTP db_3 result")
        mapping <- .pn_expand_values(db3, id_column)
        mapping$chem_source_id <- .pn_normalise_source_id(mapping[[id_column]])
        chem_map <- unique(chemicalProperties[c("chem_source_id", "ramp_id")])
        mapping <- merge(mapping, chem_map, by = "chem_source_id", sort = FALSE)
    }
    mapping <- mapping[startsWith(mapping$ramp_id, "RAMP_C_"), , drop = FALSE]
    if (!nrow(mapping)) {
        return(data.frame(
            ramp_id = character(), mz_name = character(), Adduct = character(),
            common_name = character(), stringsAsFactors = FALSE
        ))
    }
    if (!"common_name" %in% names(mapping)) mapping$common_name <- NA_character_
    if (!is.null(chemicalProperties)) {
        name_map <- chemicalProperties[
            chemicalProperties$ramp_id %in% unique(mapping$ramp_id),
            c("ramp_id", "common_name"),
            drop = FALSE
        ]
        name_map <- name_map[!duplicated(name_map$ramp_id), , drop = FALSE]
        names(name_map)[names(name_map) == "common_name"] <- ".ramp_common_name"
        mapping <- merge(
            mapping,
            name_map,
            by = "ramp_id",
            all.x = TRUE,
            sort = FALSE
        )
        missing_name <- is.na(mapping$common_name) | !nzchar(mapping$common_name)
        mapping$common_name[missing_name] <- mapping$.ramp_common_name[missing_name]
        mapping$.ramp_common_name <- NULL
    }
    unique(mapping[c("ramp_id", "mz_name", "Adduct", "common_name")])
}

.pn_prepare_metabolite_de <- function(de, db3, chemicalProperties = NULL) {
    names(de)[names(de) == "gene"] <- "mz_name"
    mapping <- .pn_prepare_metabolite_map(db3, chemicalProperties)
    merge(mapping, de, by = "mz_name", sort = FALSE)
}

.pn_best_de_by_cluster <- function(de, idColumn) {
    if (is.null(de) || !nrow(de)) {
        return(list())
    }
    id <- as.character(de[[idColumn]])
    order_index <- order(de$.cluster_key, id, de$p_val_adj, na.last = TRUE)
    de <- de[order_index, , drop = FALSE]
    key <- paste(de$.cluster_key, de[[idColumn]], sep = "\r")
    de <- de[!duplicated(key), , drop = FALSE]
    split(de, de$.cluster_key, drop = TRUE)
}

.pn_leading_edge <- function(x) {
    value <- if (is.list(x)) unlist(x, use.names = FALSE) else as.character(x)
    value <- trimws(value)
    unique(value[!is.na(value) & nzchar(value)])
}

.pn_name_map <- function(sourceTable, nodeIds, chemicalProperties = NULL) {
    value <- sourceTable[
        sourceTable$rampId %in% nodeIds & !is.na(sourceTable$commonName) &
            nzchar(trimws(sourceTable$commonName)),
        c("rampId", "commonName"),
        drop = FALSE
    ]
    value <- value[!duplicated(value$rampId), , drop = FALSE]
    result <- stats::setNames(
        as.character(value$commonName), as.character(value$rampId)
    )
    if (!is.null(chemicalProperties)) {
        .pn_assert_columns(
            chemicalProperties, c("ramp_id", "common_name"), "chem_props"
        )
        missing_ids <- setdiff(nodeIds, names(result))
        chemical <- chemicalProperties[
            chemicalProperties$ramp_id %in% missing_ids &
                !is.na(chemicalProperties$common_name) &
                nzchar(trimws(chemicalProperties$common_name)),
            c("ramp_id", "common_name"),
            drop = FALSE
        ]
        chemical <- chemical[!duplicated(chemical$ramp_id), , drop = FALSE]
        result <- c(
            result,
            stats::setNames(
                as.character(chemical$common_name), as.character(chemical$ramp_id)
            )
        )
    }
    missing_ids <- setdiff(nodeIds, names(result))
    if (length(missing_ids) && "sourceId" %in% names(sourceTable)) {
        fallback <- sourceTable[
            sourceTable$rampId %in% missing_ids &
                !is.na(sourceTable$sourceId) & nzchar(trimws(sourceTable$sourceId)),
            intersect(c("rampId", "sourceId", "IDtype"), names(sourceTable)),
            drop = FALSE
        ]
        if (nrow(fallback)) {
            id_type <- if ("IDtype" %in% names(fallback)) {
                tolower(as.character(fallback$IDtype))
            } else {
                rep("", nrow(fallback))
            }
            preferred <- c(
                "gene_symbol", "entrez", "uniprot", "ensembl", "hmdb", "chebi",
                "kegg", "lipidmaps", "refmet"
            )
            priority <- match(id_type, preferred, nomatch = length(preferred) + 1L)
            fallback <- fallback[order(priority), , drop = FALSE]
            fallback <- fallback[!duplicated(fallback$rampId), , drop = FALSE]
            label <- as.character(fallback$sourceId)
            symbol <- tolower(as.character(fallback$IDtype %||% "")) == "gene_symbol"
            label[symbol] <- sub("^[^:]+:", "", label[symbol])
            result <- c(
                result,
                stats::setNames(label, as.character(fallback$rampId))
            )
        }
    }
    result
}

.pn_build_one_network <- function(
    edges, detected, clusterKey,
    geneLookup, metaboliteLookup,
    nameMap, maxNodes,
    leadingEdge = detected,
    annotated = character(),
    annotationScores = numeric()
) {
    detected <- unique(as.character(detected))
    detected <- detected[grepl("^RAMP_[CG]_", detected)]
    annotated <- unique(as.character(annotated))
    annotated <- annotated[grepl("^RAMP_C_", annotated)]
    retained <- unique(c(detected, annotated))
    edges <- .pn_prune_edges(edges, retained, maxNodes)
    if (!nrow(edges) && !length(retained)) {
        return(list(nodes = .pn_empty_nodes(), links = .pn_empty_links()))
    }
    # A current RaMP leading-edge entity can be absent from an older stored
    # topology. Keep it as an isolated detected node instead of silently
    # discarding the current annotation result.
    node_ids <- unique(c(edges$src, edges$dest, retained))
    group <- ifelse(startsWith(node_ids, "RAMP_G_"), "rna", "mets")
    expr <- rep(NA_real_, length(node_ids))
    node_name <- unname(nameMap[node_ids])
    node_name[is.na(node_name) | !nzchar(node_name)] <- node_ids[is.na(node_name) | !nzchar(node_name)]
    display <- node_name

    gene_rows <- geneLookup[[clusterKey]]
    gene_index <- which(group == "rna")
    if (length(gene_index) && !is.null(gene_rows) && nrow(gene_rows)) {
        match_index <- match(node_ids[gene_index], gene_rows$rampId)
        present <- !is.na(match_index)
        expr[gene_index[present]] <- gene_rows$avg_log2FC[match_index[present]]
        label <- gene_rows$gene[match_index[present]]
        if ("commonName" %in% names(gene_rows)) {
            symbols <- gene_rows$commonName[match_index[present]]
            use_symbol <- !is.na(symbols) & nzchar(symbols)
            label[use_symbol] <- symbols[use_symbol]
        }
        node_name[gene_index[present]] <- label
        display[gene_index[present]] <- label
    }

    met_rows <- metaboliteLookup[[clusterKey]]
    met_index <- which(group == "mets")
    if (length(met_index) && !is.null(met_rows) && nrow(met_rows)) {
        match_index <- match(node_ids[met_index], met_rows$ramp_id)
        present <- !is.na(match_index)
        target <- met_index[present]
        source <- match_index[present]
        expr[target] <- met_rows$avg_log2FC[source]
        valid_name <- !is.na(met_rows$common_name[source]) & nzchar(met_rows$common_name[source])
        node_name[target[valid_name]] <- met_rows$common_name[source[valid_name]]
        display[target] <- paste0(
            met_rows$mz_name[source], " [", met_rows$Adduct[source], "] - ",
            node_name[target]
        )
    }

    nodes <- data.frame(
        id = node_ids,
        group = group,
        expr = expr,
        name = node_name,
        display = display,
        detected = node_ids %in% detected,
        leading_edge = node_ids %in% leadingEdge,
        annotated = node_ids %in% annotated,
        annotation_score = unname(annotationScores[node_ids]),
        stringsAsFactors = FALSE
    )
    links <- if (nrow(edges)) {
        data.frame(
            source = edges$src,
            target = edges$dest,
            colour = edges$colour,
            weight = edges$weight,
            style = edges$linetype,
            head = edges$arrowhead,
            reaction = edges$reaction_name,
            stringsAsFactors = FALSE
        )
    } else {
        .pn_empty_links()
    }
    list(nodes = nodes, links = links)
}

.pn_annotation_score_map <- function(annotationMap) {
    if (is.null(annotationMap) || !nrow(annotationMap)) {
        return(numeric())
    }
    .pn_assert_columns(annotationMap, "ramp_id", "annotation_map")
    score <- if ("Score" %in% names(annotationMap)) {
        suppressWarnings(as.numeric(annotationMap$Score))
    } else {
        rep(NA_real_, nrow(annotationMap))
    }
    ids <- as.character(annotationMap$ramp_id)
    keep <- startsWith(ids, "RAMP_C_") & is.finite(score)
    if (!any(keep)) {
        return(numeric())
    }
    score <- tapply(score[keep], ids[keep], max, na.rm = TRUE)
    score[order(-score, names(score))]
}

.pn_annotated_metabolites_by_pathway <- function(
    annotationMap, pathwayNames, pathwayRows,
    analytePathway, pathwayTable
) {
    result <- stats::setNames(
        vector("list", length(pathwayNames)),
        pathwayNames
    )
    if (is.null(annotationMap) || !nrow(annotationMap)) {
        return(result)
    }
    .pn_assert_columns(annotationMap, "ramp_id", "annotation_map")
    .pn_assert_columns(
        analytePathway, c("rampId", "pathwayRampId"), "analytehaspathway"
    )
    .pn_assert_columns(
        pathwayTable, c("pathwayRampId", "pathwayName", "sourceId"), "pathway"
    )
    score <- if ("Score" %in% names(annotationMap)) {
        suppressWarnings(as.numeric(annotationMap$Score))
    } else {
        rep(0, nrow(annotationMap))
    }
    score[!is.finite(score)] <- -Inf
    annotationMap$.pathway_score <- score
    annotationMap <- annotationMap[
        order(-annotationMap$.pathway_score, annotationMap$ramp_id), ,
        drop = FALSE
    ]

    for (pathway_name in pathwayNames) {
        selected_rows <- pathwayRows[
            tolower(as.character(pathwayRows$pathwayName)) == tolower(pathway_name), ,
            drop = FALSE
        ]
        name_keys <- unique(tolower(c(
            pathway_name,
            as.character(selected_rows$pathwayName)
        )))
        source_keys <- unique(tolower(as.character(selected_rows$sourceId)))
        pathway_ids <- unique(pathwayTable$pathwayRampId[
            tolower(as.character(pathwayTable$pathwayName)) %in% name_keys |
                tolower(as.character(pathwayTable$sourceId)) %in% source_keys
        ])
        members <- unique(as.character(analytePathway$rampId[
            analytePathway$pathwayRampId %in% pathway_ids &
                startsWith(as.character(analytePathway$rampId), "RAMP_C_")
        ]))
        mapped <- annotationMap[
            annotationMap$ramp_id %in% members, ,
            drop = FALSE
        ]
        result[[pathway_name]] <- unique(as.character(mapped$ramp_id))
    }
    result
}

.pn_topology_leading_edge <- function(ids, edges, geneIndex = NULL) {
    genes <- which(startsWith(ids, "RAMP_G_"))
    if (is.null(geneIndex) || !length(genes)) return(ids)
    mapping <- .map_gene_index(ids[genes], geneIndex, verbose = FALSE)
    nodes <- unique(c(edges$src, edges$dest))
    resolved <- lapply(seq_along(genes), function(i) {
        equivalents <- strsplit(mapping$ramp_ids[i], ";", fixed = TRUE)[[1L]]
        present <- intersect(equivalents, nodes)
        if (length(present)) present else ids[genes[i]]
    })
    unique(c(ids[-genes], unlist(resolved, use.names = FALSE)))
}

.pn_build_network_collection <- function(
    pathwayRows, pathwayNames,
    preparedEdges, clusters,
    geneLookup, metaboliteLookup,
    nameMap, maxNodes,
    metaboliteDetection = "leading_edge",
    annotatedMetabolites = NULL,
    annotationScores = numeric(),
    annotationScoreThreshold = 0.05,
    fixedCompoundIds = FALSE,
    geneIndex = NULL
) {
    networks <- lapply(clusters, function(cluster) {
        cluster_key <- tolower(as.character(cluster))
        cluster_rows <- pathwayRows[
            tolower(as.character(pathwayRows$Cluster_id)) == cluster_key, ,
            drop = FALSE
        ]
        cluster_networks <- lapply(pathwayNames, function(pathwayName) {
            pathway_key <- tolower(pathwayName)
            is_present <- any(tolower(cluster_rows$pathwayName) == pathway_key)
            if (!is_present) {
                return(list(nodes = .pn_empty_nodes(), links = .pn_empty_links()))
            }
            pathway_cluster_rows <- cluster_rows[
                tolower(cluster_rows$pathwayName) == pathway_key, ,
                drop = FALSE
            ]
            leading_edge <- unique(unlist(
                lapply(pathway_cluster_rows$leadingEdge, .pn_leading_edge),
                use.names = FALSE
            ))
            leading_edge <- .pn_topology_leading_edge(
                leading_edge, preparedEdges[[pathwayName]], geneIndex
            )
            detected <- leading_edge
            if (identical(metaboliteDetection, "annotated")) {
                pathway_annotated <- annotatedMetabolites[[pathwayName]]
                pathway_scores <- unname(annotationScores[pathway_annotated])
                score_keep <- if (fixedCompoundIds) {
                    rep(TRUE, length(pathway_annotated))
                } else {
                    is.finite(pathway_scores) &
                        pathway_scores >= annotationScoreThreshold
                }
                detected <- unique(c(
                    detected,
                    pathway_annotated[score_keep]
                ))
            }
            .pn_build_one_network(
                preparedEdges[[pathwayName]], detected, cluster_key,
                geneLookup, metaboliteLookup, nameMap, maxNodes,
                leadingEdge = leading_edge,
                annotated = annotatedMetabolites[[pathwayName]],
                annotationScores = annotationScores
            )
        })
        names(cluster_networks) <- pathwayNames
        cluster_networks
    })
    names(networks) <- clusters
    matrix_ids <- unlist(lapply(networks, function(clusterNetworks) {
        lapply(clusterNetworks, function(network) {
            network$nodes$id[is.finite(network$nodes$expr) |
                network$nodes$leading_edge | network$nodes$annotated]
        })
    }), use.names = FALSE)
    list(networks = networks, matrix_ids = unique(as.character(matrix_ids)))
}

.pn_layer_data <- function(object, assay, layer) {
    .nativeExpression(object, assay, layer)
}

.pn_coordinate_columns <- function(coordinates) {
    lower <- tolower(names(coordinates))
    pick <- function(candidates) {
        index <- match(candidates, lower, nomatch = 0L)
        index <- index[index > 0L]
        if (!length(index)) {
            return(NA_integer_)
        }
        index[[1]]
    }
    x <- pick(c("x", "imagecol", "col", "pxl_col_in_fullres"))
    y <- pick(c("y", "imagerow", "row", "pxl_row_in_fullres"))
    if (is.na(x) || is.na(y)) {
        stop("Tissue coordinates must contain x/y or imagecol/imagerow columns.")
    }
    c(x = x, y = y)
}

.pn_matrix_expression <- function(matrix, featureNames, cellNames) {
    if (!length(featureNames)) {
        return(list())
    }
    feature_index <- match(tolower(featureNames), tolower(rownames(matrix)))
    cell_index <- if (!is.null(colnames(matrix)) && length(cellNames)) {
        match(cellNames, colnames(matrix))
    } else {
        seq_len(min(ncol(matrix), length(cellNames)))
    }
    present <- !is.na(cell_index)
    lapply(feature_index, function(feature) {
        if (is.na(feature)) {
            return(NULL)
        }
        values <- rep(NA_real_, length(cellNames))
        values[present] <- as.numeric(matrix[feature, cell_index[present], drop = TRUE])
        values
    })
}

.pn_spatial_payload <- function(
    object, ident, image, matrixIds,
    geneDe, metaboliteDe,
    geneMatrix = NULL, metaboliteMatrix = NULL,
    maxSpatialPoints = 50000L
) {
    object <- .nativeSingleSample(object)
    coordinates <- .nativeCoordinates(object)
    if (!is.null(image)) .nativeImage(object, image, unique(object$sample_id))
    columns <- .pn_coordinate_columns(coordinates)
    x <- suppressWarnings(as.numeric(coordinates[[columns[["x"]]]]))
    y <- suppressWarnings(as.numeric(coordinates[[columns[["y"]]]]))
    cell_column <- match("cell", tolower(names(coordinates)), nomatch = 0L)
    coordinate_cells <- if (cell_column > 0L) {
        as.character(coordinates[[cell_column]])
    } else {
        rownames(coordinates)
    }
    metadata <- .cellMetadata(object)
    assignment <- if (length(coordinate_cells) && all(coordinate_cells %in% rownames(metadata))) {
        as.character(metadata[coordinate_cells, ident])
    } else {
        as.character(metadata[[ident]])[seq_len(nrow(coordinates))]
    }
    keep <- is.finite(x) & is.finite(y) & !is.na(assignment)
    keep_index <- which(keep)
    if (is.finite(maxSpatialPoints)) maxSpatialPoints <- as.integer(maxSpatialPoints)
    if (is.finite(maxSpatialPoints) && maxSpatialPoints > 0L &&
        length(keep_index) > maxSpatialPoints) {
        keep_index <- keep_index[unique(round(seq(
            1, length(keep_index),
            length.out = maxSpatialPoints
        )))]
    }
    x <- x[keep_index]
    y <- y[keep_index]
    assignment <- assignment[keep_index]
    if (!length(coordinate_cells)) {
        available_cells <- colnames(geneMatrix %||% metaboliteMatrix)
        coordinate_cells <- available_cells[seq_len(nrow(coordinates))]
    }
    coordinate_cells <- coordinate_cells[keep_index]
    scale01 <- function(value) {
        range <- range(value, finite = TRUE)
        if (!is.finite(diff(range)) || diff(range) == 0) {
            return(rep(0.5, length(value)))
        }
        (value - range[[1]]) / diff(range)
    }
    span <- c(diff(range(x, finite = TRUE)), diff(range(y, finite = TRUE)))
    aspectRatio <- if (all(is.finite(span)) && all(span > 0)) {
        span[[1]] / span[[2]]
    } else {
        1
    }

    gene_ids <- matrixIds[startsWith(matrixIds, "RAMP_G_")]
    met_ids <- matrixIds[startsWith(matrixIds, "RAMP_C_")]
    gene_expression <- list()
    metabolite_expression <- list()
    if (length(gene_ids) && !is.null(geneMatrix) && !is.null(geneDe) && nrow(geneDe)) {
        best <- geneDe[order(geneDe$p_val_adj), , drop = FALSE]
        best <- best[!duplicated(best$rampId), , drop = FALSE]
        features <- best$gene[match(gene_ids, best$rampId)]
        values <- .pn_matrix_expression(geneMatrix, features, coordinate_cells)
        names(values) <- gene_ids
        gene_expression <- Filter(Negate(is.null), values)
    }
    if (length(met_ids) && !is.null(metaboliteMatrix) &&
        !is.null(metaboliteDe) && nrow(metaboliteDe)) {
        best <- metaboliteDe[order(metaboliteDe$p_val_adj), , drop = FALSE]
        best <- best[!duplicated(best$ramp_id), , drop = FALSE]
        features <- best$mz_name[match(met_ids, best$ramp_id)]
        values <- .pn_matrix_expression(
            metaboliteMatrix,
            features,
            coordinate_cells
        )
        names(values) <- met_ids
        metabolite_expression <- Filter(Negate(is.null), values)
    }
    list(
        coordinates = data.frame(x = scale01(x), y = scale01(y)),
        aspect_ratio = aspectRatio,
        clusters = assignment,
        expression = list(rna = gene_expression, mets = metabolite_expression),
        displayed_points = length(keep_index),
        total_points = sum(keep)
    )
}

.pn_fc_limit <- function(deList) {
    value <- unlist(
        lapply(deList, function(x) x$avg_log2FC),
        use.names = FALSE
    )
    value <- abs(value[is.finite(value)])
    if (!length(value)) {
        return(1)
    }
    max(1, as.numeric(stats::quantile(value, 0.98, names = FALSE, type = 7)))
}

.pn_template_path <- function() {
    installed <- system.file(
        "templates",
        "pathway_network.html",
        package = "SpaMTP"
    )
    if (nzchar(installed) && file.exists(installed)) {
        return(installed)
    }
    development <- file.path("inst", "templates", "pathway_network.html")
    if (file.exists(development)) {
        return(development)
    }
    stop("Cannot locate the pathway network HTML template.")
}

.pn_render_html <- function(payload, templatePath = .pn_template_path()) {
    template <- paste(
        readLines(templatePath, warn = FALSE, encoding = "UTF-8"),
        collapse = "\n"
    )
    library_path <- system.file("vendor", "d3", "d3.min.js", package = "SpaMTP")
    if (!nzchar(library_path)) library_path <- file.path("inst", "vendor", "d3", "d3.min.js")
    if (!file.exists(library_path)) stop("Cannot locate the bundled offline D3 library.")
    script <- paste(
        readLines(library_path, warn = FALSE, encoding = "UTF-8"),
        collapse = "\n"
    )
    parts <- strsplit(template, "__SPAMTP_D3_LIBRARY__", fixed = TRUE)[[1]]
    if (length(parts) != 2L) stop("Pathway network template is missing its D3 placeholder.")
    template <- paste0(parts[1], script, parts[2])
    payload$pathways <- I(as.character(payload$pathways))
    payload$clusters <- I(as.character(payload$clusters))
    payload$networks <- unname(lapply(payload$networks, function(cluster) {
        unname(cluster)
    }))
    payload$spatial$clusters <- I(as.character(payload$spatial$clusters))
    payload$spatial$expression <- lapply(payload$spatial$expression, function(group) {
        lapply(group, function(values) I(as.numeric(values)))
    })
    payload_json <- jsonlite::toJSON(
        payload,
        auto_unbox = TRUE, dataframe = "rows", null = "null",
        na = "null", digits = 6, POSIXt = "ISO8601"
    )
    payload_json <- gsub("</", "<\\/", payload_json, fixed = TRUE)
    placeholder <- "__SPAMTP_PAYLOAD__"
    if (!grepl(placeholder, template, fixed = TRUE)) {
        stop("Pathway network template is missing its payload placeholder.")
    }
    sub(placeholder, payload_json, template, fixed = TRUE)
}
