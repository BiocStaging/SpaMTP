#' Visualise Significant Pathways
#'
#' Displays the pathway analysis results form running the
#' 'fishersPathwayAnalysis()' function
#'
#' @param object A SpatialExperiment used for fishersPathwayAnalysis.
#' @param pathwayDf Dataframe containing the pathway enrichment results (output
#'   from SpaMTP::fishersPathwayAnalysis function).
#' @param assay Primary (`"main"`) or alternative experiment name. The default
#'   is `"main"` for the primary experiment.
#' @param assayName Expression assay name within the selected experiment
#'   (default = `"counts"`), read with `SummarizedExperiment::assay()`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param minN Integer value specifying the minimum number of analytes required
#'   to be present in a pathway (default = 3).
#' @param pValThreshold The p-val cutoff to keep the pathways generated from
#'   fisher exact test (default = "0.1").
#' @param method Character string defining the statistical method used to
#'   calculate hclust (default = "ward.D2").
#' @param verbose Boolean indicating whether to show informative messages. If
#'   FALSE these messages will be suppressed (default = TRUE).
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param ... The arguments pass to stats::hclust
#'
#' @return A combined gg, ggplot object with pathway and dendrogram
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' enrichment <- fishersPathwayAnalysis(list(metabolites = c("RAMP_C_1", "RAMP_C_2")),
#'     universe = list(metabolites = paste0("RAMP_C_", 1:3)),
#'     database = spamtp_example$database, minPathSize = 1,
#'     pathwayAllInfo = TRUE, verbose = FALSE
#' )
#' visualisePathways(spe, enrichment,
#'     assay = "main", minN = 1,
#'     pValThreshold = 1, database = spamtp_example$database, verbose = FALSE
#' )
visualisePathways <- function(
    object,
    pathwayDf,
    assay = "main",
    assayName = "counts",
    minN = 3,
    pValThreshold = 0.1,
    method = "ward.D2",
    verbose = TRUE,
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    ...
) {
    no_pathways_plot <- function(reason) {
        verbose_message(messageText = reason, verbose = verbose)
        ggplot2::ggplot() +
            ggplot2::annotate(
                "text",
                x = 0,
                y = 0,
                label = reason,
                size = 5
            ) +
            ggplot2::xlim(-1, 1) +
            ggplot2::ylim(-1, 1) +
            ggplot2::labs(title = "No pathways to visualise") +
            ggplot2::theme_void()
    }

    if (is.null(minN)) {
        stop("Incorrect minimum analyte number! `minN` must be set to a value > 1. Please adjust this value accordingly ...")
    }
    pathwayDf <- pathwayDf[which(pathwayDf$analytes_in_pathways >= minN), ]
    if (nrow(pathwayDf) == 0L) {
        return(no_pathways_plot(
            paste0(
                "No pathways contain at least ",
                minN,
                " detected analytes."
            )
        ))
    }

    pathwayDf$duplicate_pathways <- rep(NA_character_, nrow(pathwayDf))
    verbose_message(
        messageText = "Reducing synonymous pathways",
        verbose = verbose
    )
    index <- seq_len(nrow(pathwayDf))
    merged_pathways <- data.frame()
    pb <- .progressBar(length(index), verbose)
    while (length(index) != 0) {
        pattern <- stringr::str_extract(
            pathwayDf$pathway_name[index[1]],
            pattern = "[A-Z][A-Z]\\([a-z0-9]"
        )
        if (length(pattern) != 0 & !is.na(pattern)) {
            pattern <- sub("\\(", "\\\\(", pattern)
            name <- strsplit(pathwayDf$pathway_name[index[1]], pattern)[[1]][1]
            full_name <- paste0(name, pattern)
            frst_ind <- which(grepl(
                pattern = full_name,
                x = pathwayDf$pathway_name
            ))
            all_pathways <- pathwayDf[frst_ind, ]
            second_ind <- which(duplicated(all_pathways$p_val))
            index <- index[-which(index %in% frst_ind)]
            if (length(second_ind) > 0) {
                name <- stringr::str_trim(name, side = "right")

                unique_pathways <- all_pathways[-second_ind, ]
                # Collect the IDs of all pathways sharing each unique p_val
                unique_pathways$duplicate_pathways <- vapply(unique_pathways$p_val, function(pVal) {
                    duplicates <- which(all_pathways$p_val == pVal)
                    paste(all_pathways$pathway_id[duplicates], collapse = ", ")
                }, character(1), USE.NAMES = FALSE)

                unique_pathways$pathway_name <- name
                merged_pathways <- rbind(merged_pathways, unique_pathways)
            } else {
                merged_pathways <- rbind(merged_pathways, all_pathways)
            }
        } else {
            merged_pathways <- rbind(merged_pathways, pathwayDf[index[1], ])
            index <- index[-1]
        }
        .progressTick(pb, nrow(pathwayDf) - length(index))
    }
    .progressDone(pb)
    merged_pathways <- merged_pathways %>%
        dplyr::filter(.data$p_val <= pValThreshold) %>%
        dplyr::mutate(signif_at_005level = ifelse(
            .data$p_val <= 0.05,
            "Significant",
            "Non-significant"
        ))

    if (nrow(merged_pathways) == 0L) {
        return(no_pathways_plot(
            paste0(
                "No pathways passed the p-value threshold (",
                pValThreshold,
                ")."
            )
        ))
    }

    if (nrow(merged_pathways) > 1L) {
        database_resources <- .spamtp_db_bundle(
            c("analytehaspathway", "pathway"),
            database = database,
            version = databaseVersion,
            source = match.arg(databaseSource),
            localDir = databaseLocalDir
        )
        analytehaspathway <- database_resources$analytehaspathway
        pathway <- database_resources$pathway

        # Adding the dendrogram

        # data(pathway, package = "SpaMTP")
        # data(analytehaspathway, package = "SpaMTP")
        jaccard_matrix <- matrix(
            nrow = nrow(merged_pathways),
            ncol = nrow(merged_pathways),
            dimnames = list(merged_pathways$pathway_id, merged_pathways$pathway_id)
        )

        for (i in seq_len(max(0L, nrow(merged_pathways) - 1L))) {
            pathway_id_i <- merged_pathways$pathway_id[i]
            pathway_content_i <- unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_i)])])
            for (j in (i + 1):nrow(merged_pathways)) {
                pathway_id_j <- merged_pathways$pathway_id[j]
                pathway_content_j <- unique(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathway$pathwayRampId[which(pathway$sourceId == pathway_id_j)])])

                jc_simi <- if (length(union(pathway_content_i, pathway_content_j))) {
                    length(intersect(pathway_content_i, pathway_content_j)) / length(union(pathway_content_i, pathway_content_j))
                } else 0
                jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- jc_simi
            }
        }
        diag(jaccard_matrix) <- 1


        # Generate a dendrogram
        clustering <- hclust(as.dist(1 - jaccard_matrix), method = method, ...)
        merged_pathways <- merged_pathways[clustering$order, , drop = FALSE]
    }

    gg_bar1 <- with(
        merged_pathways,
        ggplot2::ggplot() +
            ggplot2::geom_bar(
                ggplot2::aes(
                    x = paste0(pathway_name, "(", pathway_id, ")"),
                    y = total_in_pathways,
                    color = "Total analytes in pathway"
                ),
                stat = "identity",
                fill = NA,
                position = "dodge"
            ) +
            ggplot2::geom_bar(
                ggplot2::aes(
                    x = paste0(pathway_name, "(", pathway_id, ")"),
                    y = analytes_in_pathways,
                    colour = "Analytes detected in pathway"
                ),
                stat = "identity",
                fill = "blue",
                position = "dodge"
            ) +
            ggplot2::geom_point(
                ggplot2::aes(
                    x = paste0(pathway_name, "(", pathway_id, ")"),
                    y = 1.05 * max(total_in_pathways),
                    size = p_val,
                    fill = signif_at_005level
                ),
                color = "black",
                shape = 21,
                position = ggplot2::position_nudge(y = 0.5)
            ) +
            ggplot2::scale_size_area(max_size = 5, name = "p value") +
            ggplot2::scale_fill_manual(
                values = c(
                    "Significant" = "green",
                    "Non-significant" = "red"
                ),
                name = ""
            ) +
            ggplot2::scale_colour_manual(
                values = c(
                    "Analytes detected in pathway" = "blue",
                    "Total analytes in pathway" = "grey"
                ),
                name = ""
            ) +
            ggplot2::theme_minimal() +
            ggplot2::coord_flip() +
            ggplot2::theme(
                legend.position = "left",
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()
            ) +
            ggplot2::labs(
                y = "Number of analytes in pathway",
                x = "Pathway names"
            )
    )

    gg_bar1 <- gg_bar1 + ggplot2::scale_x_discrete(limits = paste0(
        merged_pathways$pathway_name, "(", merged_pathways$pathway_id, ")"
    ))

    # Get the raster illustration for each pathway
    image_raster <- vector("list", nrow(merged_pathways))
    verbose_message(
        messageText = "Generating spatial pathway images ... ", verbose = verbose
    )
    pb <- .progressBar(nrow(merged_pathways), verbose)
    mass_matrix <- Matrix::t(.assayData(object, assay, assayName))
    coords <- .nativeCoordinates(object, requireNames = FALSE)
    for (z in seq_len(nrow(merged_pathways))) {
        mzs <- paste0("mz-", stringr::str_extract_all(
            merged_pathways$adduct_info[z], "\\d+\\.\\d+"
        )[[1]])
        mat_ind <- which(colnames(mass_matrix) %in% mzs)
        if (length(mat_ind)) {
            image_raster[z] <- list(.pathwayFeatureRaster(
                mass_matrix[, mat_ind, drop = FALSE], coords
            ))
        }
        .progressTick(pb, z)
    }
    .progressDone(pb)
    if (length(image_raster) > 0) {
        for (k in seq_along(image_raster)) {
            if (!is.null(image_raster[[k]])) {
                gg_bar1 <- gg_bar1 + ggplot2::annotation_custom(
                    grid::rasterGrob(
                        image_raster[[k]],
                        width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc")
                    ),
                    xmin = k - 0.5,
                    xmax = k + 0.5,
                    ymin = -3.5,
                    ymax = -0.3
                )
            }
        }
    }
    gg_bar1 <- gg_bar1 + ggplot2::ylim(
        if (any(lengths(image_raster))) -3.5 else 0,
        max(merged_pathways$total_in_pathways) + 5
    )
    # gg_bar1

    if (nrow(merged_pathways) == 1L) {
        verbose_message(
            messageText = "\nOnly one pathway passed filtering; returning the pathway plot without a dendrogram.",
            verbose = verbose
        )
        return(gg_bar1)
    }

    hc <- as.dendrogram(clustering)
    segment_hc <- with(
        ggdendro::segment(ggdendro::dendro_data(hc)),
        data.frame(
            x = y,
            y = x,
            xend = yend,
            yend = xend
        )
    )
    pos_table <- with(
        ggdendro::dendro_data(hc)$labels,
        data.frame(
            y_center = x,
            gene = as.character(label),
            height = 1
        )
    )
    axis_limits <- with(
        pos_table,
        c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
    ) + 0.1 * c(-1, 1)

    plt_dendr <- ggplot2::ggplot(segment_hc) +
        ggplot2::geom_segment(ggplot2::aes(
            x = .data$x,
            y = .data$y,
            xend = .data$xend,
            yend = .data$yend
        )) +
        ggplot2::scale_x_continuous(
            expand = c(0, 0.5),
            limits = c(0, max(segment_hc$x, segment_hc$xend))
        ) +
        ggplot2::scale_y_continuous(
            expand = c(0, 0),
            breaks = pos_table$y_center,
            labels = pos_table$gene,
            limits = axis_limits,
            position = "right"
        ) +
        ggplot2::labs(
            x = "Clustering height",
            y = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

    # # Combine the dendrogram and bar plot
    # combined_plot <- grid.arrange(gg_bar1,ggplotify::as.ggplot(dendro), widths = c(5, 1), nrow =1)

    combined_plot <- cowplot::plot_grid(gg_bar1,
        plt_dendr,
        align = "h",
        rel_widths = c(6, 1)
    )
    verbose_message(messageText = "\nDone", verbose = verbose)
    # Show the plot
    return(combined_plot)
}


#' Plot significantly enriched pathways per region
#'
#' Visualisation of Set Enrichment Analysis Results from
#' `SpaMTP::findRegionalPathways()`.
#'
#' @param regpathway A dataframe generated by the
#'   `SpaMTP::findRegionalPathways()` function, containing identified regional
#'   pathways.
#' @param identColumn A character string specifying the column name of the
#'   `regpathway` dataframe containing the idents or clusters to compare
#'   (default = "Cluster_id").
#' @param selectedPathways A character vector specifying the names or IDs of
#'   pathways to be included in the analysis (e.g., `c("Amino acid metabolism",
#'   "WP1902", "Aspartate and asparagine metabolism")`). This argument is not
#'   case-sensitive (default = NULL).
#' @param sigCutoff Within-family BH-adjusted P-value cutoff for the point
#'   outline. NULL uses 0.05. Display selection does not recompute adjustment.
#' @param numDisplay An integer specifying the number of pathways to display in
#'   the plot. If set to null will plot the smaller of either 10 pathways or the
#'   number of unique pathways provided in 'regpathway' (default = NULL).
#' @param textSize Base text size in points. NULL uses 11 points. Pathway
#'   labels and legend entries use a slightly smaller size.
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for pathway lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param labelWidth Maximum target line width for pathway names, in
#'   characters.
#'   Defaults to 32. NULL keeps names on one line. Stable pathway identifiers
#'   are displayed on a separate line in either case.
#'
#' @return A `ggplot` object representing the set enrichment analysis results.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' db <- spamtp_example$database
#' regional <- data.frame(
#'     pathwayRampId = c("P1", "P2", "P3"),
#'     pathwayName = db$pathway$pathwayName, Cluster_id = "left",
#'     NES = c(1.5, -1.2, 0.6), pval = c(0.01, 0.05, 0.3),
#'     padj = c(0.03, 0.075, 0.3), size = c(4, 4, 3)
#' )
#' plotRegionalPathways(regional, database = db)
plotRegionalPathways <- function(
    regpathway,
    identColumn = "Cluster_id",
    selectedPathways = NULL,
    sigCutoff = NULL,
    numDisplay = NULL,
    textSize = NULL,
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    labelWidth = 32
) {
    if (!is.null(labelWidth) && (length(labelWidth) != 1L ||
        !is.numeric(labelWidth) || !is.finite(labelWidth) || labelWidth < 1)) {
        stop("labelWidth must be NULL or a positive number.", call. = FALSE)
    }

    resources <- .spamtp_db_bundle(c("analytehaspathway", "pathway"),
        database = database, version = databaseVersion, source = match.arg(databaseSource),
        localDir = databaseLocalDir
    )
    required <- c(
        "pathwayRampId",
        "pathwayName",
        identColumn,
        "NES",
        "size",
        "padj"
    )
    if (!all(required %in% names(regpathway))) {
        stop(
            "Regional pathway plots require pathwayRampId, pathwayName, NES, size, padj and the region column.",
            call. = FALSE
        )
    }
    if (!nrow(regpathway)) stop("No pathways to plot.", call. = FALSE)
    regpathway <- as.data.frame(regpathway)
    if (!is.null(selectedPathways)) {
        keep <- regpathway$pathwayRampId %in% selectedPathways |
            regpathway$pathwayName %in% selectedPathways
        if ("sourceId" %in% names(regpathway)) {
            keep <- keep | regpathway$sourceId %in% selectedPathways
        }
        regpathway <- regpathway[keep, , drop = FALSE]
        if (!nrow(regpathway)) stop("No selected pathways were found.", call. = FALSE)
    }
    importance <- tapply(
        abs(regpathway$NES),
        regpathway$pathwayRampId,
        function(z) max(z, na.rm = TRUE)
    )
    ids <- head(
        names(sort(importance, decreasing = TRUE)),
        numDisplay %||% 10L
    )
    if (!length(ids)) stop(
        "No selected pathways have estimable enrichment scores.",
        call. = FALSE
    )
    membership <- split(
        resources$analytehaspathway$rampId,
        resources$analytehaspathway$pathwayRampId
    )
    sets <- lapply(ids, function(id) unique(membership[[id]]))
    distance <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    if (length(ids) > 1L) {
        for (i in seq_len(length(ids) - 1L)) {
            for (j in (i + 1L):length(ids)) {
                total <- length(union(sets[[i]], sets[[j]]))
                distance[i, j] <- distance[j, i] <- if (total) 1 - length(intersect(sets[[i]], sets[[j]])) / total else 1
            }
        }
    }
    ordering <- if (length(ids) > 1L) ids[stats::hclust(stats::as.dist(distance), method = "average")$order] else ids
    d <- regpathway[regpathway$pathwayRampId %in% ids, , drop = FALSE]
    display_names <- if (is.null(labelWidth)) {
        d$pathwayName
    } else {
        stringr::str_wrap(d$pathwayName, width = labelWidth)
    }
    labels <- stats::setNames(
        paste0(display_names, "\n", d$pathwayRampId),
        d$pathwayRampId
    )
    d$.pathway <- factor(d$pathwayRampId, levels = rev(ordering))
    d$.region <- as.character(d[[identColumn]])
    d$.FDR <- factor(ifelse(
        !is.na(d$padj) & d$padj <= (sigCutoff %||% .05),
        "Within-family BH FDR",
        "Above cutoff / unavailable"
    ))
    limit <- max(abs(d$NES), na.rm = TRUE)
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$.region, .data$.pathway)) +
        ggplot2::geom_point(
            ggplot2::aes(
                fill = .data$NES,
                size = .data$size,
                colour = .data$.FDR
            ),
            shape = 21, stroke = 0.8, show.legend = TRUE
        ) +
        ggplot2::scale_fill_gradient2(
            low = "#3b5b92",
            mid = "#f7f7f3",
            high = "#b73e45",
            limits = c(-limit, limit)
        ) +
        ggplot2::scale_colour_manual(
            values = c(
                "Within-family BH FDR" = "#132f40",
                "Above cutoff / unavailable" = "#d0d7dc"
            ),
            limits = c("Within-family BH FDR", "Above cutoff / unavailable"),
            labels = c(
                paste0("BH-adjusted P <= ", sigCutoff %||% .05),
                "Above cutoff /\nunavailable"
            )
        ) +
        ggplot2::scale_size_area(max_size = 5.5) +
        ggplot2::scale_y_discrete(labels = labels[rev(ordering)]) +
        ggplot2::labs(
            x = identColumn, y = NULL, fill = "NES", size = "Tested members", colour = "Outline",
            subtitle = "Pathways ordered by membership overlap"
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_colourbar(order = 1),
            size = ggplot2::guide_legend(order = 2),
            colour = ggplot2::guide_legend(
                order = 3,
                override.aes = list(size = 3.5, fill = "white")
            )
        ) +
        ggplot2::theme_minimal(base_size = textSize %||% 11) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
            axis.text.y = ggplot2::element_text(
                size = ggplot2::rel(.85),
                lineheight = 1.05
            ),
            legend.position = "right", legend.box = "vertical",
            legend.text = ggplot2::element_text(size = ggplot2::rel(.85)),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(12, 12, 12, 12)
        )
    attr(p, "pathway_distance") <- distance
    attr(p, "selected_pathways") <- ids
    attr(p, "significance_column") <- "padj"
    p
}

# Use up to three principal components; one feature produces a grayscale image.
# Empty grid locations are transparent, so irregular samples retain their shape.
.pathwayFeatureRaster <- function(values, coords) {
    values <- as.matrix(values)
    if (!ncol(values) || !nrow(values)) return(NULL)
    if (!.allFiniteAssay(values)) stop("Pathway intensities must be finite.", call. = FALSE)
    rank <- min(3L, ncol(values), max(1L, nrow(values) - 1L))
    components <- if (nrow(values) > 1L) {
        stats::prcomp(values, rank. = rank)$x[, seq_len(rank), drop = FALSE]
    } else matrix(0, nrow = 1L, ncol = rank)
    scaled <- apply(components, 2, function(x) {
        span <- diff(range(x))
        if (span == 0) rep(0.5, length(x)) else (x - min(x)) / span
    })
    scaled <- matrix(scaled, nrow = nrow(values), ncol = rank)
    channels <- scaled[, rep(seq_len(rank), length.out = 3L), drop = FALSE]
    x <- match(coords$x, sort(unique(coords$x)))
    y <- match(coords$y, sort(unique(coords$y), decreasing = TRUE))
    pixels <- matrix("transparent", nrow = max(y), ncol = max(x))
    pixels[cbind(y, x)] <- grDevices::rgb(channels[, 1], channels[, 2], channels[, 3])
    as.raster(pixels)
}
