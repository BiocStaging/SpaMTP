#' Convert an object to a SpatialExperiment
#'
#' `asSpatialExperiment()` is the Bioconductor-native conversion entry point
#' used by SpaMTP. Continuous or unaligned raw spectra should remain in a
#' Cardinal `MSImagingArrays` or `MSImagingExperiment`; conversion is intended
#' for data with a
#' common feature-by-pixel matrix, typically after peak alignment or binning.
#'
#' @param x An object to convert.
#' @param coordinates Column names in `colData(x)` containing x/y coordinates
#'   when converting a `SummarizedExperiment`.
#' @param ... Arguments passed to the class-specific method.
#'
#' @return A `SpatialExperiment`.
#' @importClassesFrom Cardinal MSImagingArrays MSImagingExperiment
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' msi <- spatialExperimentToCardinal(spe)
#' asSpatialExperiment(msi)
methods::setGeneric(
    "asSpatialExperiment",
    function(x, ...) methods::standardGeneric("asSpatialExperiment")
)

#' Convert a SpatialExperiment to a Cardinal experiment
#'
#' @param x An object to convert.
#' @param ... Arguments passed to the class-specific method.
#'
#' @return A Cardinal `MSImagingExperiment`, or unchanged `MSImagingArrays`
#'   input.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' asCardinal(spe)
methods::setGeneric(
    "asCardinal",
    function(x, ...) methods::standardGeneric("asCardinal")
)

.requireOptionalPackage <- function(package, purpose) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(
            "Package `", package, "` is required to ", purpose,
            ". Install it with BiocManager::install(\"", package, "\").",
            call. = FALSE
        )
    }
}

.massValues <- function(rowData, featureNames = NULL) {
    rowData <- as.data.frame(rowData)
    for (column in c("mz", "raw_mz")) {
        if (column %in% colnames(rowData)) {
            values <- suppressWarnings(as.numeric(rowData[[column]]))
            if (length(values) && all(is.finite(values))) {
                return(values)
            }
        }
    }
    values <- suppressWarnings(as.numeric(sub("^mz-", "", featureNames)))
    if (!length(values) || any(!is.finite(values))) {
        stop(
            "Feature m/z values must be stored in `rowData(x)$mz`, ",
            "`rowData(x)$raw_mz`, or feature names beginning with `mz-`.",
            call. = FALSE
        )
    }
    values
}

.pixelNamesFromCardinal <- function(x, coordinates) {
    identifiers <- paste(coordinates[, "x"], coordinates[, "y"], sep = "_")
    runs <- as.character(Cardinal::run(x))
    if (length(unique(runs)) > 1L) {
        identifiers <- paste(identifiers, runs, sep = "-")
    }
    make.unique(identifiers)
}

#' Convert an aligned Cardinal experiment to SpatialExperiment
#'
#' The intensity matrix is stored as an assay, Cardinal `featureData()` becomes
#' `rowData()`, `pixelData()` becomes `colData()`, and `coord()` becomes
#' `spatialCoords()`. This conversion materialises a matter-backed intensity
#' matrix; raw continuous data should be binned with Cardinal before conversion.
#'
#' @param x A Cardinal `MSImagingExperiment` with a shared feature axis.
#' @param assayName Name of the assay in the returned object.
#' @param materialize Whether a non-matrix Cardinal assay may be materialised.
#'
#' @return A `SpatialExperiment`.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' msi <- spatialExperimentToCardinal(spe)
#' converted <- cardinalToSpatialExperiment(msi)
#' dim(converted)
cardinalToSpatialExperiment <- function(
    x,
    assayName = "counts",
    materialize = TRUE
) {
    if (methods::is(x, "MSImagingArrays")) {
        stop("Unaligned MSImagingArrays must be binned first with binSpaMTP() or Cardinal::bin().",
            call. = FALSE
        )
    }
    if (!methods::is(x, "MSImagingExperiment")) {
        stop("`x` must be a Cardinal MSImagingExperiment.", call. = FALSE)
    }
    x <- Cardinal::process(x, verbose = FALSE)
    spectra <- Cardinal::spectra(x)
    if (!is.matrix(spectra) && !inherits(spectra, "Matrix")) {
        if (!isTRUE(materialize)) {
            stop(
                "The Cardinal intensity data are file-backed. Bin the data in Cardinal ",
                "or set `materialize = TRUE` explicitly before conversion.",
                call. = FALSE
            )
        }
        spectra <- as.matrix(spectra)
    }

    featureData <- as.data.frame(Cardinal::featureData(x))
    pixelData <- as.data.frame(Cardinal::pixelData(x))
    coordinates <- as.matrix(Cardinal::coord(x))
    if (!all(c("x", "y") %in% colnames(coordinates))) {
        stop(
            "Cardinal coordinates must contain `x` and `y` columns.",
            call. = FALSE
        )
    }
    coordinates <- coordinates[, c("x", "y"), drop = FALSE]
    pixelNames <- .pixelNamesFromCardinal(x, coordinates)
    mz <- .massValues(featureData, rownames(spectra))
    featureNames <- make.unique(paste0("mz-", format(mz, digits = 15, trim = TRUE)))

    rownames(spectra) <- featureNames
    colnames(spectra) <- pixelNames
    rownames(featureData) <- featureNames
    rownames(pixelData) <- pixelNames
    rownames(coordinates) <- pixelNames
    if (!"mz" %in% colnames(featureData)) {
        featureData$mz <- mz
    }
    featureData$raw_mz <- mz
    featureData$mz_names <- featureNames

    runs <- as.character(Cardinal::run(x))
    if (!length(runs)) {
        runs <- rep("sample01", ncol(spectra))
    }
    object <- SpatialExperiment::SpatialExperiment(
        assays = stats::setNames(list(spectra), assayName),
        rowData = S4Vectors::DataFrame(featureData),
        colData = S4Vectors::DataFrame(pixelData),
        spatialCoords = coordinates,
        sample_id = runs,
        metadata = S4Vectors::metadata(x)
    )
    S4Vectors::metadata(object)$spamtp_source <- list(
        class = "MSImagingExperiment",
        converted = Sys.time()
    )
    object
}

#' @rdname asSpatialExperiment
#' @export
methods::setMethod(
    "asSpatialExperiment",
    "MSImagingExperiment",
    function(x, ...) cardinalToSpatialExperiment(x, ...)
)

#' @rdname asSpatialExperiment
#' @export
methods::setMethod(
    "asSpatialExperiment", "MSImagingArrays",
    function(x, ...) cardinalToSpatialExperiment(x, ...)
)

#' @rdname asSpatialExperiment
#' @export
methods::setMethod(
    "asSpatialExperiment",
    "SpatialExperiment",
    function(x, ...) {
        .checkDots(..., fn = "asSpatialExperiment", formalNames = "x")
        x
    }
)

#' @rdname asSpatialExperiment
#' @export
methods::setMethod(
    "asSpatialExperiment",
    "SummarizedExperiment",
    function(x, coordinates = c("x", "y"), ...) {
        .checkDots(
            ...,
            fn = "asSpatialExperiment",
            formalNames = c("x", "coordinates")
        )
        columnData <- as.data.frame(SummarizedExperiment::colData(x))
        if (length(coordinates) != 2L || anyNA(coordinates) || anyDuplicated(coordinates)) {
            stop(
                "coordinates must name two distinct x/y columns in colData.",
                call. = FALSE
            )
        }
        if (!all(coordinates %in% colnames(columnData))) {
            stop(
                "`colData(x)` must contain coordinate columns: ",
                paste(coordinates, collapse = ", "), ".",
                call. = FALSE
            )
        }
        samples <- if ("sample_id" %in% colnames(columnData)) {
            as.character(columnData$sample_id)
        } else {
            rep("sample01", ncol(x))
        }
        spatialCoordinates <- as.matrix(columnData[, coordinates, drop = FALSE])
        if (!is.numeric(spatialCoordinates) || any(!is.finite(spatialCoordinates))) {
            stop(
                "Coordinate columns must contain finite numeric values.",
                call. = FALSE
            )
        }
        colnames(spatialCoordinates) <- c("x", "y")
        object <- SpatialExperiment::SpatialExperiment(
            assays = SummarizedExperiment::assays(x),
            rowData = SummarizedExperiment::rowData(x),
            colData = SummarizedExperiment::colData(x),
            spatialCoords = spatialCoordinates,
            sample_id = samples,
            metadata = S4Vectors::metadata(x)
        )
        if (methods::is(x, "SingleCellExperiment")) {
            SingleCellExperiment::reducedDims(object) <- SingleCellExperiment::reducedDims(x)
            SingleCellExperiment::altExps(object) <- SingleCellExperiment::altExps(x)
            SingleCellExperiment::mainExpName(object) <- SingleCellExperiment::mainExpName(x)
            SingleCellExperiment::colPairs(object) <- SingleCellExperiment::colPairs(x)
            SingleCellExperiment::rowPairs(object) <- SingleCellExperiment::rowPairs(x)
        }
        object
    }
)

#' Convert a SpatialExperiment to Cardinal
#'
#' @param x A `SpatialExperiment` containing aligned MSI features.
#' @param assayName Name of the assay to convert.
#'
#' @return A Cardinal `MSImagingExperiment`.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' msi <- spatialExperimentToCardinal(spe)
#' Cardinal::mz(msi)
spatialExperimentToCardinal <- function(x, assayName = "counts") {
    if (!methods::is(x, "SpatialExperiment")) {
        stop("`x` must be a SpatialExperiment.", call. = FALSE)
    }
    assayNames <- SummarizedExperiment::assayNames(x)
    if (!assayName %in% assayNames) {
        stop(
            "Assay `", assayName, "` was not found. Available assays: ",
            paste(assayNames, collapse = ", "), ".",
            call. = FALSE
        )
    }
    featureData <- as.data.frame(SummarizedExperiment::rowData(x))
    pixelData <- as.data.frame(SummarizedExperiment::colData(x))
    coordinates <- as.data.frame(SpatialExperiment::spatialCoords(x))
    if (!all(c("x", "y") %in% colnames(coordinates))) {
        stop(
            "`spatialCoords(x)` must contain `x` and `y` columns.",
            call. = FALSE
        )
    }
    mz <- .massValues(featureData, rownames(x))
    featureExtras <- featureData[setdiff(colnames(featureData), c("mz", "raw_mz"))]
    featureFrame <- do.call(
        Cardinal::MassDataFrame,
        c(list(mz = mz), as.list(featureExtras))
    )

    runs <- if ("sample_id" %in% colnames(pixelData)) {
        factor(pixelData$sample_id)
    } else if ("run" %in% colnames(pixelData)) {
        pixelData$run
    } else {
        factor(rep("run0", ncol(x)))
    }
    pixelExtras <- pixelData[
        setdiff(colnames(pixelData), c("x", "y", "run", "sample_id"))
    ]
    pixelFrame <- do.call(
        Cardinal::PositionDataFrame,
        c(
            list(coord = coordinates[, c("x", "y"), drop = FALSE], run = runs),
            as.list(pixelExtras)
        )
    )
    intensity <- as.matrix(SummarizedExperiment::assay(x, assayName))
    dimnames(intensity) <- NULL
    Cardinal::MSImagingExperiment(
        spectraData = intensity,
        featureData = featureFrame,
        pixelData = pixelFrame,
        metadata = S4Vectors::metadata(x)
    )
}

#' @rdname asCardinal
#' @export
methods::setMethod(
    "asCardinal",
    "SpatialExperiment",
    function(x, ...) spatialExperimentToCardinal(x, ...)
)

#' @rdname asCardinal
#' @export
methods::setMethod(
    "asCardinal",
    "MSImagingExperiment",
    function(x, ...) {
        .checkDots(..., fn = "asCardinal", formalNames = "x")
        x
    }
)

#' @rdname asCardinal
#' @export
methods::setMethod(
    "asCardinal", "MSImagingArrays",
    function(x, ...) {
        .checkDots(..., fn = "asCardinal", formalNames = "x")
        x
    }
)

methods::setAs(
    "MSImagingExperiment",
    "SpatialExperiment",
    function(from) cardinalToSpatialExperiment(from)
)

methods::setAs(
    "SpatialExperiment",
    "MSImagingExperiment",
    function(from) spatialExperimentToCardinal(from)
)

#' Store a paired transcriptome as an alternative experiment
#'
#' The transcriptome must contain the same spatial pixels as the metabolomics
#' object. It is stored with `SingleCellExperiment::altExp()` rather than as a
#' second Seurat assay.
#'
#' @param x A `SpatialExperiment` containing the primary MSI experiment.
#' @param transcriptome A matrix or `SummarizedExperiment` with matching
#'   columns.
#' @param name Name used in `altExpNames(x)`.
#' @param assayName Assay name used when `transcriptome` is a matrix.
#' @param ... Reserved for class-specific methods.
#'
#' @return The updated `SpatialExperiment`.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' rna <- matrix(1:48, 2, dimnames = list(c("geneA", "geneB"), colnames(spe)))
#' paired <- addTranscriptome(spe, rna, name = "extraRNA")
#' SingleCellExperiment::altExpNames(paired)
methods::setGeneric(
    "addTranscriptome",
    function(x, transcriptome, name = "transcriptome", assayName = "counts", ...) {
        methods::standardGeneric("addTranscriptome")
    }
)

#' @rdname addTranscriptome
#' @export
methods::setMethod(
    "addTranscriptome",
    "SpatialExperiment",
    function(x, transcriptome, name = "transcriptome", assayName = "counts", ...) {
        .checkDots(
            ...,
            fn = "addTranscriptome",
            formalNames = names(formals(addTranscriptome))
        )
        if (is.matrix(transcriptome) || inherits(transcriptome, "Matrix")) {
            transcriptome <- SingleCellExperiment::SingleCellExperiment(
                assays = stats::setNames(list(transcriptome), assayName)
            )
        }
        if (!methods::is(transcriptome, "SummarizedExperiment")) {
            stop(
                "`transcriptome` must be a matrix or SummarizedExperiment.",
                call. = FALSE
            )
        }
        if (is.null(colnames(x)) || is.null(colnames(transcriptome)) ||
            anyNA(colnames(x)) || anyNA(colnames(transcriptome)) ||
            anyDuplicated(colnames(x)) || anyDuplicated(colnames(transcriptome))) {
            stop("Both experiments must have unique, non-missing column names for pixel matching.",
                call. = FALSE
            )
        }
        order <- match(colnames(x), colnames(transcriptome))
        if (anyNA(order) || length(order) != ncol(transcriptome)) {
            stop(
                "The transcriptome must contain exactly the same paired pixels as `x`.",
                call. = FALSE
            )
        }
        transcriptome <- transcriptome[, order, drop = FALSE]
        SingleCellExperiment::altExp(x, name, withColData = FALSE) <- transcriptome
        x
    }
)

#' Add an image using SpatialExperiment infrastructure
#'
#' This is a small SpaMTP generic around `SpatialExperiment::addImg()`. Images
#' are stored in `imgData()` and may remain path-backed by setting
#' `load = FALSE`.
#'
#' @param x A `SpatialExperiment`.
#' @param imageSource Local PNG/JPEG path or supported URL.
#' @param scaleFactor Image scale factor.
#' @param sampleId Sample identifier represented by the image.
#' @param imageId Image identifier.
#' @param load Whether to load the raster immediately.
#' @param ... Reserved for class-specific methods.
#'
#' @return The updated object.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' image_file <- tempfile(fileext = ".png")
#' png::writePNG(matrix(seq(0, 1, length.out = 24), 4, 6), image_file)
#' spe <- addSpatialImage(spe, image_file, imageId = "optical")
#' SpatialExperiment::imgData(spe)
#' unlink(image_file)
methods::setGeneric(
    "addSpatialImage",
    function(
        x,
                     imageSource,
                     scaleFactor = 1,
                     sampleId = unique(x$sample_id)[[1L]],
                     imageId = tools::file_path_sans_ext(basename(imageSource)),
                     load = FALSE,
                     ...
    ) {
        methods::standardGeneric("addSpatialImage")
    }
)

#' @rdname addSpatialImage
#' @export
methods::setMethod(
    "addSpatialImage",
    "SpatialExperiment",
    function(
        x,
                     imageSource,
                     scaleFactor = 1,
                     sampleId = unique(x$sample_id)[[1L]],
                     imageId = tools::file_path_sans_ext(basename(imageSource)),
                     load = FALSE,
                     ...
    ) {
        .checkDots(
            ...,
            fn = "addSpatialImage",
            formalNames = names(formals(addSpatialImage))
        )
        SpatialExperiment::addImg(
            x,
            imageSource = imageSource,
            scaleFactor = scaleFactor,
            sample_id = sampleId,
            image_id = imageId,
            load = load
        )
    }
)

#' Plot spatial features from a Bioconductor container
#'
#' This S4 plotting entry point reads feature values from an assay and pixel
#' positions from `SpatialExperiment::spatialCoords()`. Features are looked up
#' in the primary experiment first and then in the paired alternative
#' experiments (`SingleCellExperiment::altExps()`) that contain `assayName`, so
#' mass features and paired transcripts can be shown side by side. Each panel
#' has its own colour scale; the y axis is reversed to follow image
#' orientation.
#'
#' @param x A spatial experiment.
#' @param features Feature names to plot.
#' @param assayName Assay containing the feature values.
#' @param labels Optional panel titles corresponding to `features`.
#' @param pointSize Point size.
#' @param alpha Point opacity.
#' @param minCutoff,maxCutoff Numeric cutoffs, or strings such as `"q05"` and
#'   `"q95"`, recycled over features.
#' @param combine Whether to combine multiple panels.
#' @param ncol Number of columns in a combined plot.
#' @param legendTitle Title of the colour bar (default `assayName`).
#' @param ... Reserved for class-specific methods.
#'
#' @return A `ggplot`, a combined plot, or a list of plots.
#' @export
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = rbind(`mz-100` = 1:6, `mz-200` = 6:1)),
#'     spatialCoords = cbind(x = rep(1:3, 2), y = rep(1:2, each = 3))
#' )
#' plotSpatialFeature(spe, c("mz-100", "mz-200"), labels = c("m/z 100", "m/z 200"))
methods::setGeneric(
    "plotSpatialFeature",
    function(
        x,
                     features,
                     assayName = "counts",
                     labels = NULL,
                     pointSize = 1.6,
                     alpha = 1,
                     minCutoff = NA,
                     maxCutoff = NA,
                     combine = TRUE,
                     ncol = NULL,
                     legendTitle = assayName,
                     ...
    ) {
        methods::standardGeneric("plotSpatialFeature")
    }
)

.spatialCutoff <- function(values, cutoff, default) {
    if (!length(cutoff) || is.na(cutoff[[1L]])) {
        return(default)
    }
    cutoff <- cutoff[[1L]]
    if (is.character(cutoff) && grepl("^q[0-9]+([.][0-9]+)?$", cutoff)) {
        probability <- as.numeric(sub("^q", "", cutoff)) / 100
        return(as.numeric(stats::quantile(values, probability, na.rm = TRUE)))
    }
    cutoff <- suppressWarnings(as.numeric(cutoff))
    if (!is.finite(cutoff)) default else cutoff
}

# Shared look of spatial maps: axes kept, grid removed, legend below the map.
.spatialPlotTheme <- function() {
    ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            legend.position = "bottom",
            plot.title = ggplot2::element_text(size = 11, face = "bold")
        )
}

.plotSpatialValues <- function(
    x,
    values,
    labels = names(values),
    pointSize = 1.6,
    alpha = 1,
    minCutoff = NA,
    maxCutoff = NA,
    combine = TRUE,
    ncol = NULL,
    legendTitle = NULL
) {
    coordinates <- .nativeCoordinates(x, requireNames = FALSE)
    if (!is.list(values)) {
        values <- list(values)
    }
    if (length(labels) != length(values)) {
        stop("`labels` must have one value per plotted feature.", call. = FALSE)
    }
    minCutoff <- rep(minCutoff, length.out = length(values))
    maxCutoff <- rep(maxCutoff, length.out = length(values))
    plots <- Map(
        function(value, label, lower, upper) {
            if (length(value) != nrow(coordinates)) {
                stop(
                    "Every feature must have one value per spatial pixel.",
                    call. = FALSE
                )
            }
            lower <- .spatialCutoff(value, lower, min(value, na.rm = TRUE))
            upper <- .spatialCutoff(value, upper, max(value, na.rm = TRUE))
            value <- pmax(lower, pmin(upper, as.numeric(value)))
            frame <- data.frame(coordinates, intensity = value)
            plot <- ggplot2::ggplot(frame, ggplot2::aes(
                x = .data$x, y = .data$y, colour = .data$intensity
            )) +
                ggplot2::geom_point(size = pointSize, alpha = alpha) +
                ggplot2::coord_fixed() +
                ggplot2::scale_y_reverse() +
                ggplot2::scale_colour_viridis_c() +
                ggplot2::guides(colour = ggplot2::guide_colourbar(
                    title.position = "top", barwidth = grid::unit(3.5, "cm"),
                    barheight = grid::unit(0.3, "cm")
                )) +
                ggplot2::labs(title = label, colour = legendTitle %||% label) +
                .spatialPlotTheme()
            if (length(unique(frame$sample_id)) > 1L) {
                plot <- plot + ggplot2::facet_wrap(~sample_id)
            }
            plot
        },
        values,
        labels,
        minCutoff,
        maxCutoff
    )
    if (!isTRUE(combine)) {
        return(unname(plots))
    }
    if (length(plots) == 1L) {
        return(plots[[1L]])
    }
    cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv")
}

# Values for each feature from the primary experiment or, failing that, the
# one alternative experiment that holds both the feature and the assay.
.spatialFeatureValues <- function(x, features, assayName) {
    experiments <- list(main = x)
    if (inherits(x, "SingleCellExperiment")) {
        experiments <- c(experiments, as.list(SingleCellExperiment::altExps(x)))
    }
    experiments <- Filter(
        function(e) assayName %in% SummarizedExperiment::assayNames(e),
        experiments
    )
    if (!length(experiments)) {
        stop("Assay `", assayName, "` was not found.", call. = FALSE)
    }
    values <- lapply(features, function(feature) {
        found <- names(experiments)[vapply(
            experiments,
            function(e) feature %in% rownames(e), logical(1)
        )]
        if (!length(found)) {
            stop("Features were not found: ", feature, ".", call. = FALSE)
        }
        if (!"main" %in% found && length(found) > 1L) {
            stop("Feature `", feature, "` occurs in several alternative experiments: ",
                paste(found, collapse = ", "), ".",
                call. = FALSE
            )
        }
        source <- if ("main" %in% found) "main" else found
        as.numeric(SummarizedExperiment::assay(
            experiments[[source]],
            assayName
        )[feature, ])
    })
    names(values) <- features
    values
}

#' @rdname plotSpatialFeature
#' @export
methods::setMethod(
    "plotSpatialFeature",
    "SpatialExperiment",
    function(
        x,
                     features,
                     assayName = "counts",
                     labels = NULL,
                     pointSize = 1.6,
                     alpha = 1,
                     minCutoff = NA,
                     maxCutoff = NA,
                     combine = TRUE,
                     ncol = NULL,
                     legendTitle = assayName,
                     ...
    ) {
        .checkDots(
            ...,
            fn = "plotSpatialFeature",
            formalNames = names(formals(plotSpatialFeature))
        )
        values <- .spatialFeatureValues(x, features, assayName)
        .plotSpatialValues(
            x,
            values = values,
            labels = labels %||% features,
            pointSize = pointSize,
            alpha = alpha,
            minCutoff = minCutoff,
            maxCutoff = maxCutoff,
            combine = combine,
            ncol = ncol,
            legendTitle = legendTitle
        )
    }
)

#' @rdname plotSpatialFeature
#' @export
methods::setMethod(
    "plotSpatialFeature",
    "MSImagingExperiment",
    function(
        x,
        features,
        assayName = "counts",
        labels = NULL,
        pointSize = 1.6,
        alpha = 1,
        minCutoff = NA,
        maxCutoff = NA,
        combine = TRUE,
        ncol = NULL,
        legendTitle = assayName,
        ...
    ) {
        plotSpatialFeature(
            asSpatialExperiment(x, assayName = assayName),
            features = features,
            assayName = assayName,
            labels = labels,
            pointSize = pointSize,
            alpha = alpha,
            minCutoff = minCutoff,
            maxCutoff = maxCutoff,
            combine = combine,
            ncol = ncol,
            legendTitle = legendTitle,
            ...
        )
    }
)
