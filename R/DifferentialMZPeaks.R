#### SpaMTP Differential Metabolite Analysis Functions ########################################################################################################################################################################################


#' Pools Bioconductor experiment into random pools for pseudo-bulking.
#'
#' Runs pooling of a SpaMTP dataset to generate pseudo-replicates for each
#' unique identity provided.
#' This function is used by `findAllDEMs()`.
#' These random pools are technical partitions, not independent biological
#' replicates; resulting p-values must not be interpreted as biological
#' replication. For population inference aggregate independent samples instead.
#'
#' @param dataFilt A Bioconductor experiment containing count values for
#'   pooling.
#' @param idents A character string defining the idents column to pool the data
#'   against.
#' @param n An integer defining the amount of pseudo-replicates to generate for
#'   each sample (default = 3).
#' @param assay Primary (`"main"`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment. Values
#'   must be finite and non-negative for technical pooling.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param seed Numeric value used to set the seed for reproducible randomisation
#'   (default = 1234).
#'
#' @returns A SingleCellExperiment object which contains pooled
#'   (n)-pseudo-replicate counts data based on the Bioconductor experiment input
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' pooled <- runPooling(spe,
#'     idents = "region", n = 3, assay = "main",
#'     assayName = "counts", verbose = FALSE
#' )
#' dim(pooled)
runPooling <- function(
    dataFilt,
    idents,
    n,
    assay,
    assayName,
    seed = 1234,
    verbose = TRUE
) {
    .validateFeatureCount(n, "n")
    cell_metadata <- .cellMetadata(dataFilt)
    if (length(idents) != 1L || !idents %in% colnames(cell_metadata) ||
        anyNA(cell_metadata[[idents]])) {
        stop(
            "idents must name a colData column without missing values.",
            call. = FALSE
        )
    }
    samples <- unique(cell_metadata[[idents]])
    verbose_message(paste("Partitioning each group into", n, "technical pools."),
        verbose = verbose
    )
    for (i in seq_along(samples)) {
        wo <- which(cell_metadata[[idents]] == samples[i])
        if (length(wo) < n) {
            stop("Every group must have at least n pixels.", call. = FALSE)
        }
        pooled_ids <- withr::with_seed(
            seed + i,
            sample(rep(seq_len(n), length.out = length(wo)))
        )
        cell_metadata[wo, "orig.ident2"] <- paste0("group", i, "_pool", pooled_ids)
    }
    expression <- .assayData(dataFilt, assay, assayName)
    if (!nrow(expression) || !ncol(expression) || !.allFiniteAssay(expression) ||
        any(expression < 0)) {
        stop(
            "Pooling needs non-empty, finite, non-negative intensities.",
            call. = FALSE
        )
    }
    # edgeR materializes its input. Bound each call to about one million values
    # so sparse MSI input is not densified in its entirety.
    blockSize <- max(1L, floor(1e6 / ncol(expression)))
    blocks <- split(
        seq_len(nrow(expression)),
        ceiling(seq_len(nrow(expression)) / blockSize)
    )
    pooledCounts <- do.call(rbind, lapply(blocks, function(rows) {
        edgeR::sumTechReps(
            expression[rows, , drop = FALSE],
            ID = cell_metadata$orig.ident2
        )
    }))
    pools <- colnames(pooledCounts)
    pooledMetadata <- cell_metadata[match(pools, cell_metadata$orig.ident2), , drop = FALSE]
    for (column in colnames(cell_metadata)) {
        for (index in seq_along(pools)) {
            values <- cell_metadata[[column]][cell_metadata$orig.ident2 == pools[index]]
            if (length(unique(values)) != 1L) pooledMetadata[[column]][index] <- NA
        }
    }
    rownames(pooledMetadata) <- pools
    pooledMetadata$ncells <- as.integer(table(factor(cell_metadata$orig.ident2, levels = pools)))
    summed <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = pooledCounts),
        rowData = S4Vectors::DataFrame(.featureMetadata(dataFilt, assay)),
        colData = S4Vectors::DataFrame(pooledMetadata)
    )
    S4Vectors::metadata(summed)$spamtp_pooling <- list(
        group = idents, pools_per_group = n, seed = seed,
        replication = "technical partitions, not biological replicates"
    )
    return(summed)
}


#' Run limma contrasts on log-transformed technical pools
#'
#' Worker function for calculating differentially abundant metabolites per
#' pooling group.
#' Used by the compatibility `method = "technical"` branch of [findAllDEMs()].
#' Fits limma to `log2(mean intensity per pixel + 1)` and uses `limma::treat()`;
#' an edgeR DGEList stores the values and results, but no edgeR count-model
#' test or voom transformation is performed. Technical pools are not independent
#' biological replicates.
#'
#' @param pooledData A SingleCellExperiment object which contains the pooled
#'   pseudo-replicate data.
#' @param data The original SummarizedExperiment, including SingleCellExperiment
#'   or SpatialExperiment, supplying pixel groups and feature annotations.
#' @param ident Column in `colData(data)` and `colData(pooledData)` defining
#'   groups.
#' @param outputDir Existing directory for output CSV tables, or NULL to omit
#'   file output.
#' @param runName A character string defining the title of this DE analysis
#'   (will be used when saving DEMs to .csv file).
#' @param foldChangeThreshold Fold-change ratio, a finite number at least 1,
#'   passed as `log2(foldChangeThreshold)` to `limma::treat()`. For example,
#'   1.2 tests against a 20 percent change and 2 tests against a twofold change.
#' @param annotationColumn Column in the selected experiment's `rowData()`
#'   containing feature annotations, or NULL to omit annotations.
#' @param assay Primary (`"main"`) or alternative experiment name in `data`.
#'   Pooled expression is always read from `counts(pooledData)`.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param returnIndividual Boolean value defining whether to return a list of
#'   individual edgeR objects for each designated ident. If FALSE, one merged
#'   edgeR object will be returned (default = FALSE).
#'
#' @returns A modified edgeR object which contains the relative pseudo-bulking
#'   analysis outputs, including a DEMs data.frame with a list of differential
#'   expressed m/z metabolites
#' @export
#'
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' # Random pools are technical subdivisions of pixels.
#' pooled <- runPooling(spe,
#'     idents = "region", n = 3, assay = "main",
#'     assayName = "counts", verbose = FALSE
#' )
#' result <- runDE(pooled, spe,
#'     ident = "region", outputDir = NULL,
#'     runName = "example", foldChangeThreshold = 1.2,
#'     annotationColumn = "all_IsomerNames", assay = "main", verbose = FALSE
#' )
#' head(result$DEMs)
runDE <- function(
    pooledData,
    data,
    ident,
    outputDir,
    runName,
    foldChangeThreshold,
    annotationColumn,
    assay,
    returnIndividual = FALSE,
    verbose = TRUE
) {
    .validateFoldChangeThreshold(foldChangeThreshold)
    cellMetadata <- .cellMetadata(data)
    verbose_message(
        messageText = paste(
            "Running limma DE Analysis for ",
            runName,
            " -> with samples [",
            paste(unique(unlist(cellMetadata[[ident]])), collapse = ", "),
            "]"
        ),
        verbose = verbose
    )

    annotation_result <- list()

    for (condition in unique(cellMetadata[[ident]])) {
        # Create groups
        groups <- SingleCellExperiment::colData(pooledData)[[ident]]
        groups <- ifelse(groups == condition, "Comp_A", "Comp_B")

        # Extract continuous expression data (e.g., intensity matrix)
        expression_data <- SingleCellExperiment::counts(pooledData) # Or the assay holding your continuous data

        y <- edgeR::DGEList(
            SingleCellExperiment::counts(pooledData),
            samples = SingleCellExperiment::colData(pooledData)$orig.ident2,
            group = groups
        )

        y$samples$condition <- groups
        y$samples$ident <- as.character(SummarizedExperiment::colData(pooledData)[[ident]])


        # Compare average intensity per pixel, not pool totals. Region sizes
        # otherwise appear as fold changes even for identical distributions.
        pool_sizes <- SummarizedExperiment::colData(pooledData)$ncells
        if (length(pool_sizes) != ncol(expression_data) ||
            anyNA(pool_sizes) || any(!is.finite(pool_sizes) | pool_sizes <= 0)) {
            stop("Pooled data require a positive ncells value for every pool.",
                call. = FALSE)
        }
        expression_data <- log1p(.divideColumns(expression_data, pool_sizes)) / log(2)

        # Create design matrix
        design <- model.matrix(~groups)
        design[, 2] <- 1 - design[, 2] # To match your original contrast logic

        # Fit linear model
        fit <- limma::lmFit(expression_data, design)

        # Empirical Bayes moderation
        fit <- limma::eBayes(fit, robust = TRUE)

        # Use treat for log fold change threshold testing if desired
        res <- limma::treat(fit, lfc = log2(foldChangeThreshold), robust = TRUE)

        decisions <- limma::decideTests(res)
        all_decisions <- stats::setNames(
            as.integer(decisions[, ncol(decisions)]), rownames(decisions)
        )

        res_table <- limma::topTreat(
            res,
            coef = ncol(fit$design),
            n = nrow(expression_data)
        )

        res_table$regulate <- dplyr::recode(
            as.character(all_decisions[rownames(res_table)]),
            "0" = "Normal",
            "1" = "Up",
            "-1" = "Down"
        )

        # Order by p-value or FDR
        de_group_limma <- res_table[order(res_table$adj.P.Val), ]

        # Rename adj.P.Val to FDR
        colnames(de_group_limma) <- ifelse(
            colnames(de_group_limma) == "adj.P.Val",
            "FDR",
            colnames(de_group_limma)
        )

        # Add gene/metabolite names
        de_group_limma$gene <- rownames(de_group_limma)

        # Add annotations if requested
        if (!is.null(annotationColumn)) {
            annotation.data <- .featureMetadata(data, assay)
            if (!(annotationColumn %in% colnames(annotation.data))) {
                stop("The annotation column does not exist in the assay feature metadata.")
            } else {
                annotation.data_subset <- annotation.data[rownames(de_group_limma), ]
                de_group_limma$annotations <- annotation.data_subset[[annotationColumn]]
            }
        }

        # Write CSV output if directory specified
        if (!is.null(outputDir)) {
            utils::write.csv(
                de_group_limma,
                file.path(outputDir, paste0(condition, "_", runName, ".csv"))
            )
        }

        # Store results
        y$expression <- expression_data
        y$DEMs <- de_group_limma
        annotation_result[[condition]] <- y

        verbose_message(
            messageText = paste("Analysis complete for condition:", condition),
            verbose = verbose
        )
    }


    if (returnIndividual) {
        annotation_result <- lapply(names(annotation_result), function(x) {
            annotation_result[[x]]$DEMs$cluster <- x
            annotation_result[[x]]
        })

        return(annotation_result)
    } else {
        edger <- edgeR::DGEList(
            counts = annotation_result[[1]]$counts,
            samples = annotation_result[[1]]$samples
        )
        edger$samples$group <- edger$samples$ident
        edger$samples$condition <- NULL

        dems <- lapply(names(annotation_result), function(x) {
            annotation_result[[x]]$DEMs$cluster <- x
            rownames(annotation_result[[x]]$DEMs) <- NULL
            annotation_result[[x]]$DEMs
        })

        combined_dems <- do.call(rbind, dems)
        rownames(combined_dems) <- seq_along(combined_dems$cluster)

        edger$expression <- annotation_result[[1]]$expression
        edger$DEMs <- combined_dems
        return(edger)
    }
}


#' Finds differentially expressed m/z values/metabolites between all comparison
#' groups.
#'
#' @param data A SingleCellExperiment, including SpatialExperiment, for
#'   `markers` or `replicate` mode. Technical pooling also accepts a plain
#'   SummarizedExperiment.
#' @param ident Column name in `SummarizedExperiment::colData(data)` defining
#'   the groups to compare. Technical mode requires this column to be a factor.
#' @param n An integer that defines the number of pseudo-replicates (pools) per
#'   sample (default = 3).
#' @param foldChangeThreshold Technical-mode fold-change ratio, a finite number
#'   at least 1, passed as `log2(foldChangeThreshold)` to `limma::treat()`.
#'   The default 1.2 tests against a 20 percent change; 2 tests a twofold change.
#'   Not used to filter marker or biological-replicate results.
#' @param outputDir A character string defining the directory path for all
#'   output files to be stored. This path must a new directory. Else, set to
#'   NULL as default.
#' @param runName A character string defining the title of this DE analysis
#'   that will be used when saving DEMs to .csv file (default = 'findAllDEMs').
#' @param annotationColumn Optional column in the selected experiment's
#' `rowData()`, such as `"all_IsomerNames"` from `annotateSM()` (default NULL).
#' @param assay Primary (`"main"`) or alternative experiment name. The default
#'   is `"main"` for the primary experiment.
#' @param assayName Expression assay name within the selected experiment
#'   (default = `"counts"`). Choose the value scale for the analysis: technical
#'   pooling mode requires finite non-negative unlogged values; marker and
#'   replicate modes use the supplied values directly, commonly `"logcounts"`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param returnIndividual In technical mode, return one result per group
#'   instead of a merged result. Unsupported in marker/replicate modes.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param seed Numeric value used to set the seed for reproducible randomisation
#'   (default = 1234).
#'
#' @param method `markers` uses scran effect-size summaries without pixel P
#' values; `replicate` fits configured contrasts on biological-replicate means.
#'   The default is `"markers"`. The opt-in `technical` mode fits limma to
#'   log2(1 + mean intensity per pixel) in random pools. These pools are
#'   technical subdivisions and cannot establish biological replication.
#' @param replicate Complete colData field identifying biological replicates.
#' @param contrasts Named list of numerator, denominator and paired settings,
#'   required for replicate inference. At least three replicates per arm or
#'   three complete pairs are required.
#' @param minRegionSize Minimum labelled observations in a marker group.
#' @param spatialBlocks Number of spatial blocks per sample for a
#'   leave-one-block
#'   sensitivity analysis of marker effects, or zero to disable. Blocks use
#' coordinates only. They are not biological replicates or confidence intervals.
#' @returns A result list, not a modified experiment. Technical mode returns
#'   pooled values in `counts` and differential results in `DEMs` (or a list
#'   per group when `returnIndividual = TRUE`). Marker mode returns `DEMs`,
#'   `scores`, `expression`, group eligibility and provenance. Replicate mode
#'   returns `group_means` and named `tests` with status and, when eligible,
#'   limma tables and the biological-replicate design. Marker mode may return
#'   `status = "skipped"` with a reason when too few groups are eligible.
#' @export
#'
#' @examples
#' ## Paired treatment/control design with four simulated patients
#' design <- expand.grid(
#'     pixel = seq_len(3), condition = c("control", "treated"),
#'     patient_id = paste0("patient", seq_len(4))
#' )
#' expression <- withr::with_seed(1234, matrix(
#'     stats::rnorm(6 * nrow(design), mean = 5, sd = 0.3),
#'     nrow = 6,
#'     dimnames = list(
#'         paste0("feature", seq_len(6)),
#'         paste0("observation", seq_len(nrow(design)))
#'     )
#' ))
#' treated <- design$condition == "treated"
#' expression[1:2, treated] <- expression[1:2, treated] + 1.2
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(logcounts = expression),
#'     colData = S4Vectors::DataFrame(design)
#' )
#' fit <- findAllDEMs(sce,
#'     ident = "condition", assay = "main",
#'     assayName = "logcounts", method = "replicate", replicate = "patient_id",
#'     contrasts = list(treated_vs_control = list(
#'         numerator = "treated", denominator = "control", paired = TRUE
#'     )),
#'     verbose = FALSE
#' )
#' head(fit$tests$treated_vs_control$table)
findAllDEMs <- function(
    data, ident, n = 3, foldChangeThreshold = 1.2, outputDir = NULL, runName = "findAllDEMs", annotationColumn = NULL, assay = "main", assayName = "counts", returnIndividual = FALSE, verbose = TRUE, seed = 1234,
    method = c("markers", "replicate", "technical"), replicate = NULL,
    contrasts = list(), minRegionSize = 3, spatialBlocks = 0
) {
    .validateFoldChangeThreshold(foldChangeThreshold)
    method <- match.arg(method)
    if (method != "technical") {
        if (!is.null(outputDir) || returnIndividual) {
            stop(
                "For markers/replicate modes, export the returned tables explicitly.",
                call. = FALSE
            )
        }
        .requireExperiment(data, "SingleCellExperiment")
        if (length(ident) != 1L || !ident %in% colnames(SummarizedExperiment::colData(data))) {
            stop("ident must name a colData column.", call. = FALSE)
        }
        if (method == "replicate") {
            return(.replicateRegionContrasts(
                .assayData(data, assay, assayName), as.data.frame(SummarizedExperiment::colData(data)),
                ident, replicate, contrasts
            ))
        }
        return(.nativeRegionMarkers(
            data, ident, assay, assayName, minRegionSize,
            spatialBlocks, seed, annotationColumn
        ))
    }
    warning(
        "Random pixel pools are technical subdivisions, not biological replicates. Use method='markers' for region characterization or method='replicate' for biological inference.",
        call. = FALSE
    )

    if (!(is.null(outputDir))) {
        if (dir.exists(outputDir)) {
            warning("Please supply a directory path that doesn't already exist")
            stop("dir.exists(outputDir) = TRUE")
        } else {
            dir.create(outputDir)
        }
    }

    if (!is.factor(.cellMetadata(data)[[ident]])) {
        stop("ident provided is not a factor! Please convert the ident column using `factor()` ...")
    }

    # Step 1: Run Pooling to split each unique ident into 'n' number of pseudo-replicate pools
    pooled_data <- runPooling(
        data,
        ident,
        n = n,
        assay = assay,
        assayName = assayName,
        verbose = verbose,
        seed = seed
    )

    # Step 2: Run EdgeR to calculate differentially expressed m/z peaks
    DEM_results <- runDE(
        pooled_data,
        data,
        ident = ident,
        outputDir = outputDir,
        runName = runName,
        foldChangeThreshold = foldChangeThreshold,
        annotationColumn = annotationColumn,
        assay = assay,
        verbose = verbose,
        returnIndividual = returnIndividual
    )

    # Returns an EDGEr object which contains the pseudo-bulk counts in $counts and DEMs in $DEMs
    return(DEM_results)
}

.validateFoldChangeThreshold <- function(value) {
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 1) {
        stop("foldChangeThreshold must be one finite fold-change ratio >= 1.",
            call. = FALSE)
    }
    invisible(value)
}


#' Heatmap of Differentially Expressed Metabolites
#'
#' Generates a heatmap of DEMs generated from edgeR analysis run using
#' `findAllDEMs()`.
#' This function uses `pheatmap` to plot data.
#'
#' @param edgerOutput A list containing outputs from edgeR analysis (from
#'   findAllDEMs()). This includes pseudo-bulked counts and DEMs.
#' @param n A numeric integer that defines the number of UP and DOWN regulated
#'   peaks to plot (default = 5).
#' @param onlyPos Boolean indicating if only positive markers should be
#'   returned (default = FALSE).
#' @param fdrThreshold Numeric value that defines the FDR threshold to use for
#'   defining most significant results (default = 0.05).
#' @param logFCThreshold Numeric value that defines the logFC threshold to use
#'   for filtering significant results (default = 0.5).
#' @param orderBy Ranking column: `FDR`, `logFC`, `mean_auc` or `mean_cohen`.
#'   Descriptive native marker results default to `mean_auc`; other inputs use
#'   `FDR`. FDR ranking is unavailable for descriptive marker results.
#' @param scale A character string indicating if the values should be centered
#'   and scaled in either the row direction or the column direction, or none.
#'   Corresponding values are "row", "column" and "none"
#' @param color A vector of colors used in heatmap (default =
#'   grDevices::colorRampPalette(c("navy", "white", "red"))(50)).
#' @param clusterCols Boolean value determining if columns should be clustered
#'   or hclust object (default = FALSE).
#' @param clusterRows Boolean value determining if rows should be clustered or
#'   hclust object (default = TRUE).
#' @param fontsizeRow A numeric value defining the fontsize of rownames
#'   (default = 15).
#' @param fontsizeCol A numeric value defining the fontsize of colnames
#'   (default = 15).
#' @param cutreeCols A numeric value defining the number of clusters the
#'   columns are divided into, based on the hierarchical clustering(using
#'   cutree), if cols are not clustered, the argument is ignored (default = 9).
#' @param silent Boolean value indicating if the plot should not be draw
#'   (default = TRUE).
#' @param plotAnnotationsColumn Character string indicating the column name
#'   that contains the metabolite annotations to plot. Annotations = TRUE must
#'   be used in findAllDEMs() for edgeR output to include annotations. If
#'   plotAnnotationsColumn = NULL, m/z vaues will be plotted (default = NULL).
#' @param saveToPath Character string defining the full filepath and name of
#'   the plot to be saved as.
#' @param plotSaveWidth Integer value representing the width of the saved pdf
#'   plot (default = 20).
#' @param plotSaveHeight Integer value representing the height of the saved
#'   pdf plot (default = 20).
#' @param nlabelsToShow Numeric value defining the number of annotations to
#'   show per m/z (default = NULL).
#' @param annotationColors List for specifying annotation_row and
#'   annotation_col track colors manually. Check pheatmap R-Package
#'   documentation for details. If set to 'NA', default coloring will be used
#'   (default = NA).
#'
#' @returns A heatmap plot of significantly differentially expressed metabolites
#'   defined in the edgeR ouput object.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' markers <- findAllDEMs(spe,
#'     ident = "region", method = "markers",
#'     minRegionSize = 3, spatialBlocks = 0, verbose = FALSE
#' )
#' if (requireNamespace("pheatmap", quietly = TRUE)) {
#'     heatmap <- demsHeatmap(markers,
#'         n = 2, orderBy = "mean_auc",
#'         logFCThreshold = 0, silent = TRUE
#'     )
#'     heatmap
#' }
#'
#' # demsHeatmap(DEMs)
demsHeatmap <- function(
    edgerOutput,
    n = 5,
    onlyPos = FALSE,
    fdrThreshold = 0.05,
    logFCThreshold = 0.5,
    orderBy = "FDR",
    scale = "row",
    color = grDevices::colorRampPalette(c("navy", "white", "red"))(50),
    clusterCols = FALSE,
    clusterRows = TRUE,
    fontsizeRow = 15,
    fontsizeCol = 15,
    cutreeCols = 9,
    silent = TRUE,
    plotAnnotationsColumn = NULL,
    saveToPath = NULL,
    plotSaveWidth = 20,
    plotSaveHeight = 20,
    nlabelsToShow = NULL,
    annotationColors = NULL
) {
    .validateFeatureCount(n, "n")
    degs <- edgerOutput$DEMs
    markers <- inherits(edgerOutput, "spamtp_markers")
    if (markers && missing(orderBy)) orderBy <- "mean_auc"
    if (!orderBy %in% c("FDR", "logFC", "mean_auc", "mean_cohen")) {
        stop("Unknown marker ranking column.", call. = FALSE)
    }
    if (orderBy == "FDR" && !"FDR" %in% names(degs)) {
        stop(
            "These descriptive markers have no FDR. Rank by mean_auc or mean_cohen.",
            call. = FALSE
        )
    }
    if (!markers) degs <- degs[!is.na(degs$FDR) & degs$FDR < fdrThreshold, , drop = FALSE]
    selected <- lapply(unique(degs$cluster), function(g) {
        d <- degs[degs$cluster == g, , drop = FALSE]
        positive <- d[d$logFC > logFCThreshold & !is.na(d$logFC), , drop = FALSE]
        negative <- d[d$logFC < -logFCThreshold & !is.na(d$logFC), , drop = FALSE]
        pos <- head(
            positive[order(
                if (orderBy == "FDR") positive$FDR else -positive[[orderBy]],
                positive$gene,
                na.last = TRUE
            ), , drop = FALSE],
            n
        )
        neg <- if (onlyPos) NULL else head(
            negative[order(negative[[orderBy]], negative$gene, na.last = TRUE), , drop = FALSE],
            n
        )
        rbind(pos, neg)
    })
    df <- do.call(rbind, selected)
    if (!nrow(df)) stop(
        "No features meet the heatmap selection criteria.",
        call. = FALSE
    )
    # Native marker / replicate results supply exactly the matrix used by the
    # analysis. Do not re-normalize already transformed values through edgeR CPM.
    expression <- edgerOutput$expression
    if (is.null(expression)) expression <- edgeR::cpm(edgerOutput, log = TRUE)
    genes <- unique(df$gene)
    mtx <- as.matrix(expression[genes, , drop = FALSE])
    if (scale == "row") {
        variable <- apply(mtx, 1, stats::sd) > 0
        mtx <- mtx[variable, , drop = FALSE]
        genes <- genes[variable]
    }
    if (!nrow(mtx)) stop("No variable heatmap rows remain.", call. = FALSE)
    col_annot <- data.frame(
        region = edgerOutput$samples$ident,
        row.names = colnames(mtx)
    )
    labels <- genes
    if (!is.null(plotAnnotationsColumn)) {
        if (!plotAnnotationsColumn %in% names(df)) stop("Unknown annotation column.", call. = FALSE)
        annotation <- as.character(df[[plotAnnotationsColumn]][match(genes, df$gene)])
        if (!is.null(nlabelsToShow)) annotation <- labels_to_show(annotation, n = nlabelsToShow)
        good <- !is.na(annotation) & nzchar(annotation)
        labels[good] <- paste0(annotation[good], " [", genes[good], "]")
    }
    colors <- if (is.null(annotationColors)) {
        NA
    } else if (is.list(annotationColors) &&
        "region" %in% names(annotationColors)) {
        annotationColors
    } else {
        list(region = unlist(annotationColors))
    }
    .requireSuggested("pheatmap")
    p <- pheatmap::pheatmap(mtx,
        scale = scale, color = color,
        cluster_cols = if (ncol(mtx) < 2L) FALSE else clusterCols,
        annotation_col = col_annot, cluster_rows = if (nrow(mtx) < 2L) FALSE else clusterRows,
        labels_row = labels, fontsize_row = fontsizeRow, fontsize_col = fontsizeCol,
        cutree_cols = if (identical(clusterCols, FALSE)) NA else min(cutreeCols, ncol(mtx)),
        silent = silent, annotation_colors = colors
    )
    p$selected_features <- genes
    p$expression <- mtx
    p$selection <- df
    if (!is.null(saveToPath)) savePheatmapAsPDF(
        p,
        saveToPath,
        plotSaveWidth,
        plotSaveHeight
    )
    p
}


#' Saves a demsHeatmap as a PDF
#'
#' @param pheatmap A pheatmap plot object that is being saved.
#' @param filename Character string defining the full filepath and name of the
#'   plot to be saved as.
#' @param width Integer value representing the width of the saved pdf plot
#'   (default = 20).
#' @param height Integer value representing the height of the saved pdf plot
#'   (default = 20).
#'
#' @return The generated PDF path, invisibly.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("pheatmap", quietly = TRUE)) {
#'     heatmap <- pheatmap::pheatmap(matrix(1:12, 3), silent = TRUE)
#'     output <- tempfile(fileext = ".pdf")
#'     savePheatmapAsPDF(heatmap, output, width = 5, height = 4)
#'     file.info(output)$size
#'     unlink(output)
#' }
savePheatmapAsPDF <- function(pheatmap, filename, width = 20, height = 20) {
    output_file <- if (grepl("[.]pdf$", filename, ignore.case = TRUE)) filename else paste0(filename, ".pdf")
    grDevices::pdf(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(pheatmap$gtable)
    invisible(output_file)
}


########################################################################################################################################################################################################################
