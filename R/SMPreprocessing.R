#' Normalize spatial metabolomics intensity data
#'
#' SingleCellExperiment methods store the result in the standard `normcounts` or
#' `logcounts` assays. Cardinal input delegates TIC normalization to Cardinal
#' and remains file-backed. Convert Seurat inputs explicitly before analysis.
#' Cardinal normalization follows Cardinal's deferred-processing model;
#' queued steps are executed by Cardinal::process(), binSpaMTP(), or conversion
#' of an aligned experiment to SpatialExperiment.
#'
#' @param data A SingleCellExperiment (including SpatialExperiment), Cardinal
#'   MSImagingArrays/MSImagingExperiment.
#' @param normalisation.type `"TIC"` (default) or `"RC"` divides each pixel's
#'   intensities by its total and multiplies by `scale.factor`, storing
#'   `normcounts`. `"LogNormalize"` additionally applies `log2(1 + x)`
#'   and stores `logcounts`. Cardinal input supports only `"TIC"`.
#' @param scale.factor Positive target total before log transformation. When
#'   NULL, native TIC uses the number of features and RC/LogNormalize use
#'   10000. Set to 1e6 explicitly for counts/intensities per million. Zero-total
#'   pixels remain zero. Cardinal requires NULL and uses its own TIC scaling.
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param libSizes Optional pixel totals used instead of the column sums of
#'   `assayName`: one non-negative number per pixel, or the name of a numeric
#'   `colData()` column. Use this when the assay holds a feature subset but the
#'   full-spectrum total ion current (TIC) was recorded, so that subsetting does
#'   not change the normalization denominator.
#'
#' @details Sparse input stays sparse while at most half of the normalized
#'   values are nonzero; otherwise the result is stored as an ordinary
#'   `matrix`, which needs less memory than a dense-valued `dgCMatrix`.
#'   `DelayedArray` input (for example an `HDF5Array`) returns a delayed
#'   result, so file-backed intensities are not loaded into memory.
#'
#' @return An object of the same container family with normalized values.
#' @export
#'
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(c(1, 0, 3, 2, 2, 0),
#'         nrow = 2,
#'         dimnames = list(c("mz-100", "mz-200"), paste0("p", 1:3))
#'     )),
#'     spatialCoords = cbind(x = 1:3, y = 1)
#' )
#' spe <- normalizeSMData(spe, "LogNormalize", verbose = FALSE)
#' SummarizedExperiment::assay(spe, "logcounts")
#'
#' ## Normalize a feature subset by recorded full-spectrum totals
#' spe$tic <- c(10, 20, 5)
#' spe <- normalizeSMData(spe, "RC", libSizes = "tic", verbose = FALSE)
#' SummarizedExperiment::assay(spe, "normcounts")
methods::setGeneric(
    "normalizeSMData",
    function(
        data,
                     normalisation.type = "TIC",
                     scale.factor = NULL,
                     assay = "main",
                     assayName = "counts",
                     verbose = TRUE,
                     libSizes = NULL
    ) {
        methods::standardGeneric("normalizeSMData")
    }
)

.normalizationType <- function(normalisation.type) {
    match.arg(normalisation.type, c("TIC", "LogNormalize", "RC"))
}

.normalizationTotals <- function(data, counts, libSizes) {
    if (is.null(libSizes)) {
        return(as.numeric(.columnTotals(counts)))
    }
    if (is.character(libSizes)) {
        if (length(libSizes) != 1L || !libSizes %in% names(SummarizedExperiment::colData(data))) {
            stop("libSizes must name one colData() column.", call. = FALSE)
        }
        libSizes <- SummarizedExperiment::colData(data)[[libSizes]]
    }
    if (!is.numeric(libSizes) || length(libSizes) != ncol(counts) ||
        any(!is.finite(libSizes)) || any(libSizes < 0)) {
        stop(
            "libSizes must hold one finite, non-negative total per pixel.",
            call. = FALSE
        )
    }
    as.numeric(libSizes)
}

#' @rdname normalizeSMData
#' @export
methods::setMethod(
    "normalizeSMData",
    "SingleCellExperiment",
    function(
        data,
                     normalisation.type = "TIC",
                     scale.factor = NULL,
                     assay = "main",
                     assayName = "counts",
                     verbose = TRUE,
                     libSizes = NULL
    ) {
        normalisation.type <- .normalizationType(normalisation.type)
        counts <- .assayData(data, assay = assay, layer = assayName)
        if (!nrow(counts) || !ncol(counts) || anyNA(counts) ||
            any(is.infinite(counts)) || any(counts < 0)) {
            stop("Normalization requires non-empty, finite, non-negative intensities.",
                call. = FALSE
            )
        }
        librarySizes <- .normalizationTotals(data, counts, libSizes)
        if (is.null(scale.factor)) {
            scale.factor <- if (identical(normalisation.type, "TIC")) {
                nrow(counts)
            } else {
                10000
            }
        }
        if (!is.numeric(scale.factor) || length(scale.factor) != 1L ||
            !is.finite(scale.factor) || scale.factor <= 0) {
            stop(
                "scale.factor must be one positive finite number.",
                call. = FALSE
            )
        }
        experiment <- .experimentForAssay(data, assay)
        sizeFactors <- librarySizes / scale.factor
        sizeFactors[!is.finite(sizeFactors) | sizeFactors <= 0] <- 1
        transform <- if (identical(normalisation.type, "LogNormalize")) {
            "log"
        } else {
            "none"
        }
        normalized <- .divideColumns(counts, sizeFactors)
        if (identical(transform, "log")) {
            normalized <- log1p(normalized) / log(2)
        }
        normalized <- .compactAssay(normalized)
        outputName <- if (identical(transform, "log")) "logcounts" else "normcounts"
        SummarizedExperiment::assay(
            experiment,
            outputName,
            withDimnames = FALSE
        ) <- normalized
        if (methods::is(experiment, "SingleCellExperiment")) {
            SingleCellExperiment::sizeFactors(experiment) <- sizeFactors
        }
        S4Vectors::metadata(experiment)$spamtp_normalization <- list(
            method = normalisation.type,
            log_base = if (identical(transform, "log")) 2 else NA_real_,
            source_assay = assayName,
            output_assay = outputName,
            scale_factor = scale.factor,
            library_sizes = if (is.null(libSizes)) "column sums" else if (is.character(libSizes)) libSizes else "supplied"
        )
        verbose_message(
            paste0(
                "Stored ",
                normalisation.type,
                " values in assay `",
                outputName,
                "`."
            ),
            verbose = verbose
        )
        .replaceExperiment(data, experiment, assay)
    }
)

#' @rdname normalizeSMData
#' @export
methods::setMethod(
    "normalizeSMData",
    "MSImagingExperiment",
    function(
        data,
                     normalisation.type = "TIC",
                     scale.factor = NULL,
                     assay = "main",
                     assayName = "counts",
                     verbose = TRUE,
                     libSizes = NULL
    ) {
        .normalizeCardinal(
            data,
            normalisation.type,
            scale.factor,
            verbose,
            libSizes
        )
    }
)

.normalizeCardinal <- function(
    data, normalisation.type, scale.factor, verbose,
    libSizes = NULL
) {
    normalisation.type <- .normalizationType(normalisation.type)
    if (!identical(normalisation.type, "TIC") || !is.null(scale.factor) ||
        !is.null(libSizes)) {
        stop(
            "Cardinal input supports TIC with Cardinal's default scaling only. ",
            "For RC, LogNormalize, a custom scale.factor or libSizes, bin the Cardinal ",
            "object and convert it to SpatialExperiment first.",
            call. = FALSE
        )
    }
    Cardinal::normalize(data, method = "tic", verbose = verbose)
}

#' @rdname normalizeSMData
#' @export
methods::setMethod(
    "normalizeSMData", "MSImagingArrays",
    function(
        data, normalisation.type = "TIC", scale.factor = NULL,
                     assay = "main", assayName = "counts", verbose = TRUE, libSizes = NULL
    ) {
        .normalizeCardinal(
            data,
            normalisation.type,
            scale.factor,
            verbose,
            libSizes
        )
    }
)

#' @rdname normalizeSMData
#' @export
methods::setMethod(
    "normalizeSMData",
    "ANY",
    function(
        data,
                     normalisation.type = "TIC",
                     scale.factor = NULL,
                     assay = "main",
                     assayName = "counts",
                     verbose = TRUE,
                     libSizes = NULL
    ) {
        .requireExperiment(data, "SingleCellExperiment")
    }
)

#' Performs TMM normalization between categories based on a specified ident
#'
#' This function is mainly used for normalising a merged Bioconductor experiment
#' containing multiple samples.
#'
#' @param combined.obj Bioconductor experiment that contains groups being
#'   normalized.
#' @param ident Character string defining the column name or ident group to
#'   normalize between.
#' @param refIdent Character string specifying one class/group type to use as a
#'   reference for TMM normalisation.
#' @param normalisation.type Character string defining the normalization method
#'   to run. Options are either c("CPM", "TIC", "LogNormalize") which represent
#'   counts per million (CPM), Total Ion Current (TIC) normalization or Log
#'   Normalization, respectively (default = "CPM").
#' @param CPM.scale.factor Numeric value that sets the scale factor for
#'   pixel/spot level normalization. Following normalization the total intensity
#'   value across each pixel will equal this value (default = 1e6).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#'
#' @return Bioconductor experiment with count values normalised and corrected
#'   for between categories
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' normalized <- tmmNormalize(spe, ident = "region", refIdent = "left")
#' SummarizedExperiment::assayNames(normalized)
tmmNormalize <- function(
    combined.obj,
    ident,
    refIdent,
    normalisation.type = "CPM",
    CPM.scale.factor = 1e6,
    assay = "main",
    assayName = "counts") {
    .requireExperiment(combined.obj, "SummarizedExperiment")
    normalisation.type <- match.arg(
        normalisation.type,
        c("CPM", "TIC", "LogNormalize")
    )
    if (length(CPM.scale.factor) != 1L || !is.finite(CPM.scale.factor) ||
        CPM.scale.factor <= 0) {
        stop(
            "CPM.scale.factor must be one positive finite number.",
            call. = FALSE
        )
    }
    metadata <- .cellMetadata(combined.obj)
    if (!ident %in% colnames(metadata)) {
        stop(
            "Grouping column `",
            ident,
            "` was not found in colData.",
            call. = FALSE
        )
    }
    groups <- factor(metadata[[ident]])
    if (anyNA(groups)) {
        stop("The grouping column must not contain missing values.",
            call. = FALSE
        )
    }
    if (nlevels(groups) <= 1L) {
        stop("TMM normalization requires at least two groups.", call. = FALSE)
    }
    if (!refIdent %in% levels(groups)) {
        stop("`refIdent` is not present in the grouping column.", call. = FALSE)
    }
    counts <- .assayData(combined.obj, assay, assayName)
    if (!nrow(counts) || anyNA(counts) ||
        any(is.infinite(counts)) || any(counts < 0)) {
        stop("TMM requires finite, non-negative intensities.", call. = FALSE)
    }
    groupIndices <- split(seq_len(ncol(counts)), groups)
    pseudoBulk <- vapply(
        groupIndices,
        function(index) Matrix::rowSums(counts[, index, drop = FALSE]),
        numeric(nrow(counts))
    )
    rownames(pseudoBulk) <- rownames(counts)
    if (any(colSums(pseudoBulk) <= 0)) {
        stop(
            "Every TMM group must have positive total intensity.",
            call. = FALSE
        )
    }
    factors <- edgeR::normLibSizes(
        pseudoBulk,
        method = "TMM",
        refColumn = match(refIdent, colnames(pseudoBulk))
    )
    librarySizes <- Matrix::colSums(counts)
    target <- if (identical(normalisation.type, "TIC")) nrow(counts) else CPM.scale.factor
    multipliers <- target / (librarySizes * factors[as.character(groups)])
    multipliers[!is.finite(multipliers)] <- 0
    normalized <- counts %*% Matrix::Diagonal(x = multipliers)
    outputName <- if (identical(normalisation.type, "LogNormalize")) "logcounts" else "normcounts"
    if (identical(normalisation.type, "LogNormalize")) {
        normalized <- log1p(normalized) / log(2)
    }
    dimnames(normalized) <- dimnames(counts)
    experiment <- .experimentForAssay(combined.obj, assay)
    SummarizedExperiment::assay(experiment, outputName) <- normalized
    if (methods::is(experiment, "SingleCellExperiment")) {
        SingleCellExperiment::sizeFactors(experiment) <- ifelse(
            multipliers > 0, 1 / multipliers, 1
        )
    }
    S4Vectors::metadata(experiment)$tmm <- list(
        log_base = if (identical(normalisation.type, "LogNormalize")) 2 else NA_real_,
        group = ident,
        reference = refIdent,
        factors = factors,
        target = target
    )
    return(.replaceExperiment(combined.obj, experiment, assay))
}


########### Pre-processing plots #############


#' Helper function for QC plots by generating intensity count data
#'
#' @param data A SingleCellExperiment, including SpatialExperiment.
#' @param group.by Name of the `colData()` column to group by (default = NULL).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param bottom.cutoff Numeric value defining the percent of data to exclude
#'   for the lower end of the distribution. A bottom.cutoff = 0.05 will remove
#'   the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for
#'   the upper end of the distribution. A top.cutoff = 0.05 will remove the top
#'   5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis
#'   values (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default =
#'   FALSE).

#'
#' @return A data.frame containing the relative transformed and sum counts
#'   required for various QC plots
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' statPlot(spe, group.by = "region")
statPlot <- function(
    data,
    group.by = NULL,
    assay = "main",
    assayName = "counts",
    bottom.cutoff = NULL,
    top.cutoff = NULL,
    log.data = FALSE,
    verbose = FALSE
) {
    assayMatrix <- .assayData(data, assay = assay, layer = assayName)
    metadata <- .cellMetadata(data)
    groups <- if (is.null(group.by)) {
        factor(rep("data", ncol(assayMatrix)))
    } else {
        if (!group.by %in% colnames(metadata)) {
            stop(
                "Grouping column `",
                group.by,
                "` was not found.",
                call. = FALSE
            )
        }
        factor(metadata[[group.by]])
    }
    data_list <- split(seq_len(ncol(assayMatrix)), groups)
    df <- data.frame(mz = rownames(assayMatrix))
    rownames(df) <- df$mz

    for (dataset in names(data_list)) {
        indices <- data_list[[dataset]]
        df[[dataset]] <- Matrix::rowSums(
            assayMatrix[, indices, drop = FALSE]
        )
    }

    df2 <- tidyr::pivot_longer(
        df,
        cols = names(data_list),
        names_to = "var",
        values_to = "x"
    )

    df2 <- df2 %>% dplyr::arrange(.data$x)

    if (!(is.null(bottom.cutoff))) {
        df2 <- df2 %>%
            dplyr::group_by(.data$var) %>%
            dplyr::mutate(bottom_cutoff = stats::quantile(.data$x, bottom.cutoff))

        df2 <- df2 %>%
            dplyr::group_by(.data$var) %>%
            dplyr::filter(.data$x >= .data$bottom_cutoff)

        verbose_message(
            message_text = paste0(
                "Removing bottom ",
                bottom.cutoff * 100,
                "% of datapoints"
            ),
            verbose = verbose
        )
    }

    if (!(is.null(top.cutoff))) {
        df2 <- df2 %>%
            dplyr::group_by(.data$var) %>%
            dplyr::mutate(top_cutoff = quantile(.data$x, 1 - top.cutoff))

        df2 <- df2 %>%
            dplyr::group_by(.data$var) %>%
            dplyr::filter(.data$x <= .data$top_cutoff)

        verbose_message(
            message_text = paste0("Removing top ", top.cutoff * 100, "% of datapoints"),
            verbose = verbose
        )
    }

    if (log.data) {
        df2$x <- log10(df2$x + 1)
    }

    return(df2)
}


#' Generates a ridge plot of spatial metabolic intensity data
#'
#' @param data Bioconductor experiment containing the metabolomic intensity
#'   data.
#' @param group.by Name of the `colData()` column to group by (default = NULL).
#' @param mzs Vector of characters defining which features (m/z's) to label on
#'   the plot. If `NULL` no features will be labeled (default = NULL).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param title Character string of the plot title (default = "RidgePlot").
#' @param x.lab Character string of the x-axis label (default = "intensity").
#' @param y.lab Character string of the y-axis label (default = "var").
#' @param bottom.cutoff Numeric value defining the percent of data to exclude
#'   for the lower end of the distribution. A bottom.cutoff = 0.05 will remove
#'   the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for
#'   the upper end of the distribution. A top.cutoff = 0.05 will remove the top
#'   5% of data point (default = NULL).
#' @param bins number of bins to group
#' @param log.data Boolean value indicating whether to log transform the y-axis
#'   values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This
#'   vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default =
#'   FALSE).
#'
#' @return A `ggplot` ridge-plot object.
#'
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' mzRidgePlot(spe, group.by = "region", mzs = c("ion1", "ion3"))
mzRidgePlot <- function(
    data,
    group.by = NULL,
    mzs = NULL,
    assay = "main",
    assayName = "counts",
    title = "RidgePlot",
    x.lab = "intensity",
    y.lab = "var",
    bottom.cutoff = NULL,
    top.cutoff = NULL,
    bins = 1000,
    log.data = FALSE,
    cols = NULL,
    verbose = FALSE
) {
    data <- statPlot(
        data = data,
        group.by = group.by,
        assay = assay,
        assayName = assayName,
        bottom.cutoff = bottom.cutoff,
        top.cutoff = top.cutoff,
        log.data = log.data,
        verbose = verbose
    )

    ridge_plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(y = .data$var, x = .data$x, fill = .data$var)
    ) +
        ggridges::geom_density_ridges(
            alpha = 0.6,
            stat = "binline",
            bins = bins
        ) +
        ggridges::theme_ridges() +
        ggplot2::labs(title = title, x = x.lab, y = y.lab)

    if (!(is.null(cols))) {
        ridge_plot <- ridge_plot + ggplot2::scale_fill_manual(values = cols)
    }

    if (!(is.null(mzs))) {
        # gets feature data
        label_data <- data[data$mz %in% mzs, ]
        label_data$Features <- factor(label_data$mz, levels = mzs)

        # Plot feature points for both features
        ridge_plot <- ridge_plot +
            ggplot2::geom_point(
                data = label_data,
                ggplot2::aes(
                    x = .data$x,
                    y = .data$var,
                    color = .data$Features
                ),
                size = 3,
                shape = 19
            ) +
            ggplot2::guides(
                fill = ggplot2::guide_legend(override.aes = list(shape = NA)),
                color = ggplot2::guide_legend(override.aes = list(size = 3))
            )
    }

    return(ridge_plot)
}


#' Generates a violin plot of spatial metabolic intensity data
#'
#' @param data Bioconductor experiment containing the metabolomic intensity
#'   data.
#' @param group.by Name of the `colData()` column to group by (default = NULL).
#' @param mzs Vector of characters defining which features (m/z's) to label on
#'   the plot. If `NULL` no features will be labeled (default = NULL).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param title Character string of the plot title (default = "Expression").
#' @param x.lab Character string of the x-axis label (default = "var").
#' @param y.lab Character string of the y-axis label (default = "intensity").
#' @param show.points Boolean value describing whether to show each individual
#'   data point (default = TRUE).
#' @param bottom.cutoff Numeric value defining the percent of data to exclude
#'   for the lower end of the distribution. A bottom.cutoff = 0.05 will remove
#'   the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for
#'   the upper end of the distribution. A top.cutoff = 0.05 will remove the top
#'   5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis
#'   values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This
#'   vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default =
#'   FALSE).
#'
#' @return A `ggplot` violin-plot object.
#'
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' if (requireNamespace("scater", quietly = TRUE)) {
#'     mzViolinPlot(spe, group.by = "region", mzs = c("ion1", "ion3"))
#' }
mzViolinPlot <- function(
    data,
    group.by = NULL,
    mzs = NULL,
    assay = "main",
    assayName = "counts",
    title = "Expression",
    x.lab = "var",
    y.lab = "intensity",
    show.points = TRUE,
    bottom.cutoff = NULL,
    top.cutoff = NULL,
    log.data = FALSE,
    cols = NULL,
    verbose = FALSE
) {
    if (!is.null(mzs)) {
        .requireExperiment(data, "SingleCellExperiment")
        experiment <- .experimentForAssay(data, assay)
        if (!methods::is(experiment, "SingleCellExperiment")) {
            experiment <- methods::as(experiment, "SingleCellExperiment")
        }
        SummarizedExperiment::colData(experiment) <- SummarizedExperiment::colData(data)
        features <- vapply(
            mzs,
            function(feature) {
                if (feature %in% rownames(experiment)) feature else findNearestMZ(data, feature, assay)
            },
            character(1)
        )
        if (!assayName %in% SummarizedExperiment::assayNames(experiment)) {
            stop("Assay `", assayName, "` was not found.", call. = FALSE)
        }
        .requireSuggested("scater")
        pointFunction <- if (isTRUE(show.points)) ggplot2::geom_point else NULL
        plot <- scater::plotExpression(
            experiment,
            features = features,
            x = group.by,
            exprs_values = assayName,
            point_fun = pointFunction
        ) +
            ggplot2::labs(title = title, x = x.lab, y = y.lab)
        if (!is.null(cols)) {
            plot <- plot + ggplot2::scale_fill_manual(values = cols)
        }
        return(plot)
    }

    data <- statPlot(
        data = data,
        group.by = group.by,
        assay = assay,
        assayName = assayName,
        bottom.cutoff = bottom.cutoff,
        top.cutoff = top.cutoff,
        log.data = log.data,
        verbose = verbose
    )

    violin_plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(x = .data$var, y = .data$x, fill = .data$var)
    ) +
        ggplot2::geom_violin() +
        ggplot2::theme_classic() +
        ggplot2::labs(title = title, x = x.lab, y = y.lab)

    if (show.points) {
        violin_plot <- violin_plot + ggplot2::geom_point(
            position = ggplot2::position_jitter(width = 0.2),
            alpha = 0.5
        )
    }

    if (!(is.null(cols))) {
        violin_plot <- violin_plot + ggplot2::scale_fill_manual(values = cols)
    }

    if (!(is.null(mzs))) {
        # gets feature data
        label_data <- data[data$mz %in% mzs, ]
        label_data$Features <- factor(label_data$mz, levels = mzs)

        # Plot feature points for both features
        violin_plot <- violin_plot +
            ggplot2::geom_point(
                data = label_data,
                ggplot2::aes(
                    x = .data$var,
                    y = .data$x,
                    color = .data$Features
                ),
                size = 3,
                shape = 19
            ) +
            ggplot2::guides(
                fill = ggplot2::guide_legend(override.aes = list(shape = NA)),
                color = ggplot2::guide_legend(override.aes = list(size = 3))
            )
    }


    return(violin_plot)
}


#' Generates a Boxplot of spatial metabolic intensity data
#'
#' @param data Bioconductor experiment containing the metabolomic intensity
#'   data.
#' @param group.by Name of the `colData()` column to group by (default = NULL).
#' @param mzs Vector of characters defining which features (m/z's) to label on
#'   the plot. If `NULL` no features will be labeled (default = NULL).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param title Character string of the plot title (default = "BoxPlot").
#' @param x.lab Character string of the x-axis label (default = "var").
#' @param y.lab Character string of the y-axis label (default = "intensity").
#' @param show.points Boolean value describing whether to show each individual
#'   data point (default = TRUE).
#' @param bottom.cutoff Numeric value defining the percent of data to exclude
#'   for the lower end of the distribution. A bottom.cutoff = 0.05 will remove
#'   the bottom 5% of data point (default = NULL).
#' @param top.cutoff Numeric value defining the percent of data to exclude for
#'   the upper end of the distribution. A top.cutoff = 0.05 will remove the top
#'   5% of data point (default = NULL).
#' @param log.data Boolean value indicating whether to log transform the y-axis
#'   values (default = FALSE).
#' @param cols Vector of strings defining the colours to use for plotting. This
#'   vector should match the length of unique groups (default = NULL).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default =
#'   FALSE).
#'
#' @return A `ggplot` box-plot object.
#'
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' mzBoxPlot(spe, group.by = "region", mzs = c("ion1", "ion3"))
mzBoxPlot <- function(
    data,
    group.by = NULL,
    mzs = NULL,
    assay = "main",
    assayName = "counts",
    title = "BoxPlot",
    x.lab = "var",
    y.lab = "intensity",
    show.points = TRUE,
    bottom.cutoff = NULL,
    top.cutoff = NULL,
    log.data = FALSE,
    cols = NULL,
    verbose = FALSE
) {
    data <- statPlot(
        data = data,
        group.by = group.by,
        assay = assay,
        assayName = assayName,
        bottom.cutoff = bottom.cutoff,
        top.cutoff = top.cutoff,
        log.data = log.data,
        verbose = verbose
    )

    box_plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(x = .data$var, y = .data$x, fill = .data$var)
    ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_classic() +
        ggplot2::labs(title = title, x = x.lab, y = y.lab)

    if (show.points) {
        box_plot <- box_plot + ggplot2::geom_point(
            position = ggplot2::position_jitter(width = 0.2),
            alpha = 0.5
        )
    }

    if (!(is.null(cols))) {
        box_plot <- box_plot + ggplot2::scale_fill_manual(values = cols)
    }

    if (!(is.null(mzs))) {
        # gets feature data
        label_data <- data[data$mz %in% mzs, ]
        label_data$Features <- factor(label_data$mz, levels = mzs)

        # Plot feature points for both features
        box_plot <- box_plot +
            ggplot2::geom_point(
                data = label_data,
                ggplot2::aes(
                    x = .data$var,
                    y = .data$x,
                    color = .data$Features
                ),
                size = 3,
                shape = 19
            ) +
            ggplot2::guides(
                fill = ggplot2::guide_legend(override.aes = list(shape = NA)),
                color = ggplot2::guide_legend(override.aes = list(size = 3))
            )
    }

    return(box_plot)
}
