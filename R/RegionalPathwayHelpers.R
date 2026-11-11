# Development assistance: see inst/COPYRIGHTS.

.standardizeRegionalDE <- function(data, index) {
    columns <- colnames(data)
    valid <- any(c("avg_log2FC", "logFC") %in% columns) &&
        any(c("p_val_adj", "FDR") %in% columns) &&
        all(c("cluster", "gene") %in% columns)
    if (!valid) {
        stop(
            "DE dataframe [", index, "] requires 'cluster', 'gene', ",
            "'avg_log2FC' or 'logFC', and 'p_val_adj' or 'FDR'.",
            call. = FALSE
        )
    }
    columns[columns == "logFC"] <- "avg_log2FC"
    columns[columns == "FDR"] <- "p_val_adj"
    colnames(data) <- columns
    data
}

.regionalRegulation <- function(ranks, ids, prefix) {
    selected <- names(ranks) %in% ids & grepl(prefix, names(ranks))
    paste0(ifelse(ranks[which(selected)] >= 0, "\u2191", "\u2193"),
        collapse = ";"
    )
}

.regionalLeadingEdge <- function(
    ids, analyteTypes, metaboliteDe, geneSource, ranks, geneResults
) {
    ids <- unique(unlist(ids))
    metabolites <- genes <- NULL
    if ("metabolites" %in% analyteTypes) {
        metabolites <- metaboliteDe[
            which(metaboliteDe$ramp_id %in% ids),
        ] %>%
            dplyr::mutate(adduct_info = paste0(
                .data$observed_mz, "[", .data$Adduct, "]"
            )) %>%
            dplyr::filter(!duplicated(.data$adduct_info))
    }
    if ("genes" %in% analyteTypes) {
        genes <- geneSource[which(
            geneSource$rampId %in% ids & grepl("RAMP_G", geneSource$rampId)
        ), ]
    }
    data.frame(
        adduct_info = paste0(metabolites$adduct_info, collapse = ";"),
        leadingEdge_metabolites = paste0(
            sub(";.*", "", metabolites$IsomerNames),
            collapse = ";"
        ),
        leadingEdge_metabolites_id = paste0(
            metabolites$chem_source_id,
            collapse = ";"
        ),
        leadingEdge_genes = paste0(genes$commonName, collapse = ";"),
        met_regulation = if ("metabolites" %in% analyteTypes) {
            .regionalRegulation(ranks, ids, "RAMP_C")
        } else "",
        rna_regulation = if (geneResults) {
            .regionalRegulation(ranks, ids, "RAMP_G")
        } else ""
    )
}
