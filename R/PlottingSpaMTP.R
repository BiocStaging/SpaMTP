# MSI plotting and feature queries for native Bioconductor containers.

#' Find the measured features nearest to target m/z values
#'
#' @param data A SummarizedExperiment, including SingleCellExperiment or
#'   SpatialExperiment, with feature m/z values in `rowData()` or feature names.
#' @param mz Numeric vector of target m/z values.
#' @param assay Primary (`"main"` or NULL) or alternative experiment name.
#'   The mass lookup uses the selected experiment's `rowData()`, not an
#'   expression layer; see [experimentAccess].
#'
#' @returns A character vector with the feature identifier of the nearest
#'   measured m/z for each target, not a numeric mass. Ties select the first
#'   feature in the stored order.
#' @export
#'
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(1, 3, 2,
#'         dimnames = list(c("mz-100.01", "mz-200.02", "mz-300.03"), c("p1", "p2"))
#'     )),
#'     rowData = S4Vectors::DataFrame(mz = c(100.01, 200.02, 300.03)),
#'     spatialCoords = cbind(x = 1:2, y = 1)
#' )
#' findNearestMZ(spe, c(200, 299.9))
findNearestMZ <- function(data, mz, assay = NULL) {
    featureMetadata <- .featureMetadata(data, assay)
    featureNames <- rownames(featureMetadata)
    if (is.null(featureNames)) {
        featureNames <- rownames(.assayData(data, assay, "counts"))
    }
    numbers <- .massValues(featureMetadata, featureNames)
    if (!is.numeric(mz) || !length(mz) || any(!is.finite(mz))) {
        stop("`mz` must contain finite numeric values.", call. = FALSE)
    }
    featureNames[vapply(
        mz,
        function(target) which.min(abs(numbers - target)),
        integer(1)
    )]
}


#' Sum selected metabolite features per pixel
#'
#' The sum is stored in colData without modifying the expression assays.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param mzs Character feature IDs to sum. Duplicates count once.
#' @param assay Primary or alternative experiment.
#' @param assayName Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param binName Output colData column name. Existing columns are protected.
#' @return The updated SpatialExperiment.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' binned <- binMetabolites(spe, c("ion1", "ion2"), binName = "combined")
#' head(binned$combined)
binMetabolites <- function(
    data, mzs, assay = "main", assayName = "counts",
    binName = "Binned_Metabolites"
) {
    data <- .nativeSpatialObject(data)
    values <- .nativeExpression(data, assay, assayName)
    if (!length(mzs) || anyNA(mzs) || any(!mzs %in% rownames(values))) {
        stop("Supply existing feature IDs.", call. = FALSE)
    }
    if (length(binName) != 1L || is.na(binName) || !nzchar(binName) ||
        binName %in% colnames(SummarizedExperiment::colData(data))) {
        stop("Choose a new colData column name.", call. = FALSE)
    }
    SummarizedExperiment::colData(data)[[binName]] <-
        as.numeric(Matrix::colSums(values[unique(mzs), , drop = FALSE]))
    data
}

#' Convert point symbols to squares without changing coordinates or scales
#'
#' For combined feature panels, request combine=FALSE from the original plot
#' function and pass the resulting list. This function changes point symbols,
#' not tissue geometry or the physical pixel footprint.
#' @param plot A ggplot with point layers, or a list of such plots.
#' @return The ggplot or list, retaining images, coordinates, scales and labels.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' pixelPlot(imageMZPlot(spe, "ion1"))
pixelPlot <- function(plot) {
    if (is.list(plot) && !inherits(plot, "ggplot")) {
        return(lapply(plot, pixelPlot))
    }
    if (!inherits(plot, "ggplot")) stop("Supply a ggplot or list of ggplots.", call. = FALSE)
    points <- which(vapply(
        plot$layers, function(layer) inherits(layer$geom, "GeomPoint"),
        logical(1)
    ))
    if (!length(points)) {
        stop(
            "No point layers; use combine=FALSE for multi-panel feature plots.",
            call. = FALSE
        )
    }
    for (i in points) {
        layer <- plot$layers[[i]]
        filled <- !is.null(layer$mapping$fill %||% plot$mapping$fill)
        parameters <- utils::modifyList(
            layer$aes_params,
            list(shape = if (filled) 22 else 15)
        )
        plot$layers[[i]] <- ggplot2::ggproto(NULL, layer, aes_params = parameters)
    }
    plot
}


#' Plot MSI features from a Bioconductor spatial container
#'
#' These functions read the primary experiment or a named altExp without
#' modifying the object. Numeric m/z queries select the nearest measured mass;
#' character feature IDs select exact rows. A mass tolerance sums all peaks
#' within that distance of the selected measured mass. Overlapping annotation
#' windows count each peak once. Multiple samples are plotted separately.
#'
#' @param object A SpatialExperiment or aligned Cardinal MSImagingExperiment.
#' @param mzs Numeric m/z values or character feature IDs.
#' @param plusminus Optional non-negative absolute m/z tolerance.
#' @param fov Optional sample IDs from colData(object)$sample_id.
#' @param cols Gradient colours.
#' @param size Point size.
#' @param minCutoff,maxCutoff One cutoff or one per query; numeric or
#'   quantile strings such as "q10". NA uses the observed range.
#' @param splitBy Optional colData column for separate panels within samples.
#' @param alpha Point opacity, or a two-value range mapped to expression.
#' @param borderColor,borderSize Point border colour and width.
#' @param darkBackground Use a black plot background.
#' @param crop Crop to the selected pixels.
#' @param cells Optional pixel IDs to display.
#' @param scale Colour limits: "feature" across samples, "all" across all
#'   queries, or "none" separately for each panel.
#' @param axes Show axes.
#' @param combine Combine panels; FALSE returns a list of ggplots.
#' @param coordFixed Use a fixed coordinate aspect ratio.
#' @param assay Primary experiment ("main") or altExp name.
#' @param assayName Expression assay name, e.g. "counts" or "logcounts".
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param plotPixel Draw square rather than circular point symbols.
#' @param ... Former Seurat segmentation, molecule and blending arguments
#'   are unsupported and cause an error; convert and plot these separately.
#' @return A ggplot, or a list of panels when combine is FALSE.
#' @export
#' @examples
#' x <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = rbind(a = 1:4, b = 4:1)),
#'     rowData = S4Vectors::DataFrame(mz = c(100, 200)),
#'     colData = S4Vectors::DataFrame(row.names = paste0("p", 1:4)),
#'     spatialCoords = cbind(x = 1:4, y = c(0, 1, 0, 1))
#' )
#' imageMZPlot(x, mzs = 100)
imageMZPlot <- function(
    object, mzs, plusminus = NULL, fov = NULL,
    cols = c("lightgrey", "firebrick1"), size = 0.5,
    minCutoff = NA, maxCutoff = NA, splitBy = NULL, alpha = 1,
    borderColor = "white", borderSize = 0, darkBackground = TRUE,
    crop = FALSE, cells = NULL, scale = c("feature", "all", "none"),
    axes = FALSE, combine = TRUE, coordFixed = TRUE,
    assay = "main", assayName = "counts", plotPixel = FALSE, ...
) {
    if (length(list(...))) stop(
        "Seurat-specific plotting arguments are unsupported.",
        call. = FALSE
    )
    object <- .nativeSpatialObject(object)
    values <- .mzPlotValues(object, mzs, assay, assayName, plusminus)
    .plotNativeExpression(object, values,
        sampleIds = fov, colors = cols,
        pointSize = size, minCutoff = minCutoff, maxCutoff = maxCutoff,
        splitBy = splitBy, alpha = alpha, borderColor = borderColor,
        stroke = borderSize %||% 0, dark = darkBackground, crop = crop,
        cells = cells, scale = match.arg(scale), axes = axes, combine = combine,
        coordFixed = coordFixed, pixel = plotPixel
    )
}

#' Plot annotated metabolites without an optical image
#'
#' All matching peaks are summed per metabolite. A peak shared by overlapping
#' mass-tolerance windows is counted once. No temporary metadata is written.
#' @inheritParams imageMZPlot
#' @param metabolites Metabolite names to search in feature annotations.
#' @param columnName Annotation column in rowData.
#' @param plotExact Match complete semicolon-delimited names
#'   (case-insensitive).
#' @return A ggplot or list of panels.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' imageMZAnnotationPlot(spe, "Glucose")
imageMZAnnotationPlot <- function(
    object, metabolites, plusminus = NULL, fov = NULL,
    cols = c("lightgrey", "firebrick1"), size = 0.5,
    minCutoff = NA, maxCutoff = NA, splitBy = NULL, alpha = 1,
    borderColor = "white", borderSize = 0, darkBackground = TRUE,
    crop = FALSE, cells = NULL, scale = c("feature", "all", "none"),
    axes = FALSE, combine = TRUE, coordFixed = TRUE,
    assay = "main", assayName = "counts", columnName = "all_IsomerNames",
    plotExact = TRUE, plotPixel = FALSE, ...
) {
    if (length(list(...))) stop(
        "Seurat-specific plotting arguments are unsupported.",
        call. = FALSE
    )
    object <- .nativeSpatialObject(object)
    values <- .mzPlotValues(
        object, NULL, assay, assayName, plusminus,
        metabolites, columnName, plotExact
    )
    .plotNativeExpression(object, values,
        sampleIds = fov, colors = cols,
        pointSize = size, minCutoff = minCutoff, maxCutoff = maxCutoff,
        splitBy = splitBy, alpha = alpha, borderColor = borderColor,
        stroke = borderSize %||% 0, dark = darkBackground, crop = crop,
        cells = cells, scale = match.arg(scale), axes = axes, combine = combine,
        coordFixed = coordFixed, pixel = plotPixel
    )
}

#' Plot MSI features with optional optical images
#'
#' Images come from imgData and use its scaleFactor to place the raster in
#' spatialCoords units. Each sample/image is a separate panel. With images=NULL
#' all samples are plotted without an image.
#' @inheritParams imageMZPlot
#' @param images Optional image IDs in imgData.
#' @param keepScale Colour scaling: "feature", "all" or "none"; NULL means
#'   "none".
#' @param ncol Number of columns in combined plots.
#' @param ptSizeFactor Point size.
#' @param imageAlpha Optical image opacity.
#' @param stroke Point border width.
#' @param interactive Return Plotly widgets rather than static plots.
#' @param information Deprecated Seurat hover data; non-NULL is unsupported.
#' @return A ggplot or Plotly widget; a list when combine=FALSE.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' spatialMZPlot(spe, c("ion1", "ion3"))
spatialMZPlot <- function(
    object, mzs, plusminus = NULL, images = NULL, crop = TRUE,
    assay = "main", assayName = "counts", keepScale = "feature",
    minCutoff = NA, maxCutoff = NA, ncol = NULL, combine = TRUE,
    ptSizeFactor = 1.6, alpha = c(1, 1), imageAlpha = 1,
    stroke = 0.25, interactive = FALSE, information = NULL) {
    if (!is.null(information)) {
        stop(
            "information is a retired Seurat hover-data argument.",
            call. = FALSE
        )
    }
    object <- .nativeSpatialObject(object)
    values <- .mzPlotValues(object, mzs, assay, assayName, plusminus)
    .plotNativeExpression(object, values,
        images = images, crop = crop,
        scale = keepScale %||% "none", minCutoff = minCutoff,
        maxCutoff = maxCutoff, ncol = ncol, combine = combine,
        pointSize = ptSizeFactor, alpha = alpha, imageAlpha = imageAlpha,
        stroke = stroke, interactive = interactive
    )
}

#' Plot annotated metabolites with optional optical images
#'
#' Matching peaks are summed using the same rules as imageMZAnnotationPlot.
#' @inheritParams spatialMZPlot
#' @inheritParams imageMZAnnotationPlot
#' @return A ggplot or Plotly widget; a list when combine=FALSE.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' spatialMZAnnotationPlot(spe, "Glucose")
spatialMZAnnotationPlot <- function(
    object, metabolites, plusminus = NULL, images = NULL, crop = TRUE,
    assay = "main", assayName = "counts", keepScale = "feature",
    minCutoff = NA, maxCutoff = NA, ncol = NULL, combine = TRUE,
    ptSizeFactor = 1.6, alpha = c(1, 1), imageAlpha = 1,
    stroke = 0.25, interactive = FALSE, information = NULL,
    columnName = "all_IsomerNames", plotExact = TRUE) {
    if (!is.null(information)) {
        stop(
            "information is a retired Seurat hover-data argument.",
            call. = FALSE
        )
    }
    object <- .nativeSpatialObject(object)
    values <- .mzPlotValues(
        object, NULL, assay, assayName, plusminus,
        metabolites, columnName, plotExact
    )
    .plotNativeExpression(object, values,
        images = images, crop = crop,
        scale = keepScale %||% "none", minCutoff = minCutoff,
        maxCutoff = maxCutoff, ncol = ncol, combine = combine,
        pointSize = ptSizeFactor, alpha = alpha, imageAlpha = imageAlpha,
        stroke = stroke, interactive = interactive
    )
}


#' Plot mean mass spectra from a native spatial experiment
#'
#' Feature masses come from rowData, so arbitrary feature IDs are supported.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param groupBy Optional colData column defining overlaid groups.
#' @param splitBy Optional colData column defining separate panels; mutually
#'   exclusive with group.by.
#' @param cols Optional group colours.
#' @param assay Primary or alternative experiment.
#' @param assayName Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param labelAnnotations Display annotations instead of m/z labels.
#' @param annotationColumn Annotation column in rowData.
#' @param mzLabels Numeric m/z values to label using nearest measured peaks.
#' @param metaboliteLabels Metabolite names to label; mutually exclusive with
#'   mz.labels.
#' @param xlab,ylab Axis labels.
#' @param massRange Optional two-value displayed mass range.
#' @param yLim Optional y-axis limits.
#' @param labelCex Text size.
#' @param labelAdj,labelOffset Vertical and horizontal text justification.
#' @param labelCol Text colour.
#' @param nlabelsToShow Maximum annotations displayed per feature.
#' @return A ggplot.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' massIntensityPlot(spe, groupBy = "region")
massIntensityPlot <- function(
    data, groupBy = NULL, splitBy = NULL, cols = NULL, assay = "main",
    assayName = "counts", labelAnnotations = FALSE, annotationColumn = "all_IsomerNames",
    mzLabels = NULL, metaboliteLabels = NULL, xlab = "m/z", ylab = "intensity",
    massRange = NULL, yLim = NULL, labelCex = 5, labelAdj = -1,
    labelOffset = 0, labelCol = "#eb4034", nlabelsToShow = NULL
) {
    data <- .nativeSpatialObject(data)
    expression <- .nativeExpression(data, assay, assayName)
    metadata <- .featureMetadata(data, assay)
    masses <- .massValues(metadata, rownames(expression))
    if (!is.null(groupBy) && !is.null(splitBy)) {
        stop("Choose groupBy or splitBy, not both.", call. = FALSE)
    }
    if (!is.null(mzLabels) && !is.null(metaboliteLabels)) {
        stop("Choose mzLabels or metaboliteLabels, not both.", call. = FALSE)
    }
    grouping <- groupBy %||% splitBy
    if (is.null(grouping)) {
        groups <- rep("Sample mean", ncol(data))
    } else {
        if (length(grouping) != 1L || !grouping %in% colnames(SummarizedExperiment::colData(data))) {
            stop("Grouping column was not found in colData.", call. = FALSE)
        }
        groups <- as.character(SummarizedExperiment::colData(data)[[grouping]])
        if (anyNA(groups)) stop("Grouping contains missing values.", call. = FALSE)
    }
    selected <- integer()
    if (!is.null(mzLabels)) {
        if (!is.numeric(mzLabels) || any(!is.finite(mzLabels))) {
            stop("mzLabels must be finite numeric masses.", call. = FALSE)
        }
        selected <- vapply(
            mzLabels,
            function(mass) which.min(abs(masses - mass)),
            integer(1)
        )
    }
    if (!is.null(metaboliteLabels)) {
        selected <- unique(unlist(lapply(metaboliteLabels, function(metabolite) {
            rows <- searchAnnotations(data, metabolite, assay,
                searchExact = TRUE,
                columnName = annotationColumn
            )
            if (!nrow(rows)) stop("No annotation for: ", metabolite, call. = FALSE)
            match(rows$mz_names, rownames(expression))
        })))
    }
    labels <- rep(NA_character_, nrow(expression))
    if (length(selected)) {
        if (labelAnnotations) {
            if (is.null(annotationColumn) || !annotationColumn %in% names(metadata)) {
                stop(
                    "Annotation column was not found in rowData.",
                    call. = FALSE
                )
            }
            labels[selected] <- as.character(metadata[[annotationColumn]][selected])
            if (!is.null(nlabelsToShow)) {
                labels[selected] <- labels_to_show(labels[selected], n = nlabelsToShow)
            }
        } else {
            labels[selected] <- paste0("m/z ", format(masses[selected], trim = TRUE))
        }
    }
    frames <- lapply(unique(groups), function(group) {
        data.frame(
            mass = masses, intensity = as.numeric(Matrix::rowMeans(
                expression[, groups == group, drop = FALSE]
            )),
            variable = group, annotation = labels
        )
    })
    frame <- do.call(rbind, frames)
    if (any(!is.finite(frame$intensity))) stop("Expression must be finite.", call. = FALSE)
    if (!is.null(massRange)) {
        if (length(massRange) != 2L || any(!is.finite(massRange)) || diff(massRange) < 0) {
            stop(
                "massRange must contain increasing finite limits.",
                call. = FALSE
            )
        }
        frame <- frame[frame$mass >= massRange[1L] & frame$mass <= massRange[2L], ]
    }
    annotations <- frame[!is.na(frame$annotation), , drop = FALSE]
    if (is.null(splitBy)) {
        annotations <- annotations[order(annotations$mass, -annotations$intensity), ]
        annotations <- annotations[!duplicated(annotations$mass), ]
    }
    plot <- ggplot2::ggplot(
        frame,
        ggplot2::aes(
            x = .data$mass,
            y = .data$intensity,
            colour = .data$variable
        )
    ) +
        ggplot2::geom_line() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = xlab, y = ylab, colour = grouping %||% "Group")
    if (nrow(annotations)) {
        plot <- plot + ggplot2::geom_text(
            data = annotations, ggplot2::aes(label = .data$annotation), size = labelCex,
            vjust = labelAdj, hjust = labelOffset, colour = labelCol
        )
    }
    if (!is.null(splitBy)) plot <- plot + ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y")
    if (!is.null(cols)) plot <- plot + ggplot2::scale_colour_manual(values = cols)
    if (!is.null(yLim)) plot <- plot + ggplot2::coord_cartesian(ylim = yLim)
    plot
}


#' Compare two spatial features in 3D
#'
#' Expression comes from the primary or alternative experiment; numeric
#' colData columns can also be plotted. Coordinates and images use
#' SpatialExperiment accessors. Only one sample is displayed at a time.
#' @param data A SpatialExperiment or aligned Cardinal experiment.
#' @param features One or two feature IDs or numeric colData column names.
#'   A single feature is repeated for comparison across two experiments.
#' @param assays One or two primary/altExp names.
#' @param assayNames One or two expression assay names.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param betweenLayerHeight Distance between layers.
#' @param names Optional layer labels.
#' @param size Marker size.
#' @param colPalette One or two Plotly colour palettes.
#' @param xAxisLabel,yAxisLabel,zAxisLabel Axis labels.
#' @param showXTicks,showYTicks,showZTicks Show axis tick labels.
#' @param showImage Optional image ID from imgData.
#' @param plotHeight,plotWidth Widget dimensions.
#' @param downscaleImage Optional positive integer raster sampling stride.
#' @param sampleId Sample ID; required when data contains multiple samples.
#' @return A Plotly widget.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' plot3DFeature(spe,
#'     features = c("GENE1", "ion1"),
#'     assays = c("transcriptome", "main")
#' )
plot3DFeature <- function(
    data, features, assays = c("transcriptome", "main"), assayNames = "counts",
    betweenLayerHeight = 100, names = NULL, size = 3, colPalette = "Reds",
    xAxisLabel = "x", yAxisLabel = "y", zAxisLabel = "z",
    showXTicks = FALSE, showYTicks = FALSE, showZTicks = FALSE,
    showImage = NULL, plotHeight = 800, plotWidth = 1500,
    downscaleImage = NULL, sampleId = NULL
) {
    data <- .nativeSingleSample(data, sampleId)
    arguments <- list(features, assays, assayNames, colPalette)
    if (any(!lengths(arguments) %in% c(1L, 2L))) {
        stop(
            "Supply one or two features, assays, assay names and palettes.",
            call. = FALSE
        )
    }
    features <- rep(features, length.out = 2L)
    assays <- rep(assays, length.out = 2L)
    assayNames <- rep(assayNames, length.out = 2L)
    palettes <- rep(as.list(colPalette), length.out = 2L)
    if (!is.null(names) && !length(names) %in% c(1L, 2L)) {
        stop("Supply one or two layer labels.", call. = FALSE)
    }
    labels <- if (is.null(names)) {
        paste(features, assays, sep = " | ")
    } else {
        rep(names, length.out = 2L)
    }
    coordinates <- .nativeCoordinates(data)
    metadata <- .cellMetadata(data)
    plot <- plotly::plot_ly(height = plotHeight, width = plotWidth)
    for (i in seq_len(2L)) {
        if (features[i] %in% colnames(metadata)) {
            values <- metadata[[features[i]]]
        } else {
            expression <- .nativeExpression(data, assays[i], assayNames[i])
            if (!features[i] %in% rownames(expression)) stop("Unknown feature: ", features[i], call. = FALSE)
            values <- as.numeric(expression[features[i], ])
        }
        if (!is.numeric(values) || any(!is.finite(values))) {
            stop("Features must contain finite numeric values.", call. = FALSE)
        }
        plot <- plotly::add_trace(plot,
            x = coordinates$x, y = coordinates$y,
            z = rep((2L - i) * betweenLayerHeight, nrow(coordinates)),
            text = coordinates$cell, type = "scatter3d", mode = "markers",
            name = labels[i], marker = list(
                color = values, size = size,
                colorscale = palettes[[i]], showscale = TRUE,
                colorbar = list(x = (i - 1L) * 1.1, title = list(text = labels[i]))
            )
        )
    }
    if (!is.null(showImage)) {
        image <- .nativeImage(data, showImage, unique(data$sample_id))
        stride <- downscaleImage %||% 1L
        if (length(stride) != 1L || !is.finite(stride) || stride < 1 || stride != floor(stride)) {
            stop("downscaleImage must be a positive integer.", call. = FALSE)
        }
        rows <- seq.int(1L, nrow(image$raster), by = stride)
        columns <- seq.int(1L, ncol(image$raster), by = stride)
        # Match the row-major colour vector returned by the raster subset.
        grid <- expand.grid(column = columns, row = rows)
        colours <- as.vector(image$raster[rows, columns, drop = FALSE])
        plot <- plotly::add_trace(plot,
            x = (grid$column - 0.5) / image$scaleFactor,
            y = (grid$row - 0.5) / image$scaleFactor,
            z = rep(-betweenLayerHeight, nrow(grid)),
            type = "scatter3d", mode = "markers", name = showImage,
            marker = list(color = colours, size = size, showscale = FALSE)
        )
    }
    plotly::layout(plot, scene = list(
        aspectmode = "data",
        xaxis = list(title = xAxisLabel, showticklabels = showXTicks),
        yaxis = list(
            title = yAxisLabel, showticklabels = showYTicks,
            autorange = "reversed"
        ),
        zaxis = list(title = zAxisLabel, showticklabels = showZTicks)
    ))
}


#' Interactively inspect mass windows in a spatial experiment
#'
#' The app sums every measured feature within the selected mass window without
#' changing the input assays or annotations. The mass-axis controls use rowData
#' masses, never mean intensities. No browser is launched by this function.
#' @param obj A SpatialExperiment or aligned Cardinal experiment.
#' @param assay Primary or alternative experiment.
#' @param assayName Expression assay.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param image Optional image ID from imgData.
#' @param sampleId Optional sample ID. NULL displays separate sample panels.
#' @return A Shiny application object; use shiny::runApp to launch.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' app <- interactiveSpatialPlot(spe)
#' if (interactive()) shiny::runApp(app)
interactiveSpatialPlot <- function(
    obj, assay = "main", assayName = "counts",
    image = NULL, sampleId = NULL
) {
    obj <- if (is.null(sampleId)) {
        .nativeSpatialObject(obj)
    } else {
        .nativeSingleSample(obj, sampleId)
    }
    expression <- .nativeExpression(obj, assay, assayName)
    masses <- .massValues(.featureMetadata(obj, assay), rownames(expression))
    if (!length(masses)) stop("No features to plot.", call. = FALSE)
    means <- as.numeric(Matrix::rowMeans(expression))
    if (any(!is.finite(means))) stop("Expression must be finite.", call. = FALSE)
    centre <- stats::median(masses)
    width <- max(diff(range(masses)) / 100, 0.001)
    ui <- shiny::fluidPage(
        shiny::titlePanel("Interactive Spatial Intensity Plot"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::numericInput("curve_center", "m/z value",
                    value = centre,
                    min = min(masses), max = max(masses)
                ),
                shiny::radioButtons("bin_mode", "Window units",
                    choices = c("Absolute (m/z)" = "mz", "Relative (ppm)" = "ppm"),
                    selected = "mz"
                ),
                shiny::numericInput(
                    "curve_width",
                    "Half-width",
                    value = width,
                    min = 0
                ),
                shiny::numericInput(
                    "spot_size",
                    "Point size",
                    value = 1.6,
                    min = 0
                ),
                shiny::numericInput(
                    "x_min",
                    "Minimum m/z",
                    value = min(masses)
                ),
                shiny::numericInput("x_max", "Maximum m/z", value = max(masses))
            ),
            shiny::mainPanel(
                plotly::plotlyOutput("plot"),
                shiny::verbatimTextOutput("selected_indices"),
                shiny::plotOutput("spatial_plot")
            )
        )
    )
    server <- function(input, output, session) {
        window <- shiny::reactive({
            shiny::req(
                input$curve_center,
                !is.null(input$curve_width),
                input$bin_mode
            )
            .mzWindow(input$curve_center, input$curve_width, input$bin_mode)
        })
        selectedRows <- shiny::reactive({
            bounds <- window()
            which(masses >= bounds[1L] & masses <= bounds[2L])
        })
        selectedValues <- shiny::reactive({
            as.numeric(Matrix::colSums(expression[selectedRows(), , drop = FALSE]))
        })
        output$plot <- plotly::renderPlotly({
            bounds <- window()
            plot <- plotly::plot_ly(
                x = masses,
                y = means,
                type = "bar",
                name = "Mean intensity"
            )
            plotly::layout(plot,
                xaxis = list(title = "m/z", range = c(input$x_min, input$x_max)),
                yaxis = list(title = "Mean intensity"),
                shapes = list(list(
                    type = "rect", x0 = bounds[1L], x1 = bounds[2L],
                    y0 = 0, y1 = 1, yref = "paper", fillcolor = "lightblue",
                    opacity = 0.3, line = list(width = 0)
                ))
            )
        })
        output$selected_indices <- shiny::renderText({
            if (!length(selectedRows())) {
                return("No features in this mass window.")
            }
            paste(
                "Selected features:",
                paste(rownames(expression)[selectedRows()], collapse = ", ")
            )
        })
        output$spatial_plot <- shiny::renderPlot({
            shiny::validate(shiny::need(
                length(selectedRows()) > 0,
                "No selected features."
            ))
            values <- matrix(selectedValues(),
                nrow = 1L,
                dimnames = list(
                    paste0("m/z ", input$curve_center, " window"),
                    colnames(obj)
                )
            )
            .plotNativeExpression(
                obj,
                values,
                images = image,
                pointSize = input$spot_size
            )
        })
    }
    shiny::shinyApp(ui, server)
}

.mzWindow <- function(centre, width, units) {
    if (length(centre) != 1L || !is.finite(centre) || length(width) != 1L ||
        !is.finite(width) || width < 0 || !units %in% c("mz", "ppm")) {
        stop(
            "Supply a finite m/z centre and non-negative window width.",
            call. = FALSE
        )
    }
    if (units == "ppm") {
        if (centre <= 0) stop("ppm windows require a positive mass.", call. = FALSE)
        width <- width * centre * 1e-6
    }
    c(centre - width, centre + width)
}
