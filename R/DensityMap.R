# Development assistance: see inst/COPYRIGHTS.

#' Export an interactive spatial detection-density viewer
#'
#' Select a mass in the spectrum or feature selector to inspect the spatial
#' density of its nonzero pixels. The spatial kernel counts detections without
#' intensity weights. The spectrum overlay is an intensity-weighted mass
#' density. This exploratory display is not a spatial significance test.
#' Plotly is supplied by the R package and d3 by SpaMTP. The HTML and its
#' accompanying `density_map_files` directory work offline; keep them together
#' when copying the viewer.
#' @param object A SpatialExperiment or aligned Cardinal experiment.
#' @param assay Primary or alternative experiment.
#' @param assayName Non-negative intensity assay. See [experimentAccess] for
#'   the experiment/matrix distinction.
#' @param folder Output directory. Defaults to the session's temporary
#'   directory; supply a persistent directory to retain the viewer.
#' @param sampleId Sample ID; required for multiple samples.
#' @param ... Reserved; unused arguments cause an error.
#' @return Invisibly, the path to mzs_density_map.html.
#' @export
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = rbind(ion1 = 1:6, ion2 = c(0, 2, 0, 3, 0, 1))),
#'     rowData = S4Vectors::DataFrame(mz = c(100, 200)),
#'     colData = S4Vectors::DataFrame(row.names = paste0("p", 1:6)),
#'     spatialCoords = cbind(x = rep(1:3, 2), y = rep(1:2, each = 3))
#' )
#' folder <- tempfile("density-")
#' viewer <- densityMap(spe, folder = folder)
#' basename(viewer)
#' unlink(folder, recursive = TRUE)
densityMap <- function(
    object, assay = "main", assayName = "counts",
    folder = tempdir(), sampleId = NULL, ...
) {
    if (length(list(...))) stop("Unused densityMap arguments.", call. = FALSE)
    object <- .nativeSingleSample(object, sampleId)
    payload <- .densityMapData(object, assay, assayName)
    spectrum <- .densityMapSpectrum(payload)
    script <- readLines(system.file("htmlwidgets", "density-map.js",
        package = "SpaMTP"
    ), warn = FALSE)
    json <- jsonlite::toJSON(payload,
        digits = NA, na = "null",
        dataframe = "rows", auto_unbox = TRUE
    )
    json <- gsub("<", "\\u003c", json, fixed = TRUE)
    widget <- htmlwidgets::onRender(spectrum, paste(script, collapse = "\n"),
        data = json
    )
    widget$dependencies <- c(widget$dependencies, list(htmltools::htmlDependency(
        name = "spamtp-d3", version = "7.9.0",
        src = system.file("vendor", "d3", package = "SpaMTP"),
        script = "d3.min.js"
    )))
    if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE)) {
        stop("Cannot create output directory.", call. = FALSE)
    }
    output <- file.path(folder, "mzs_density_map.html")
    htmlwidgets::saveWidget(widget, output,
        selfcontained = FALSE,
        libdir = "density_map_files", title = "Spatial detection density"
    )
    invisible(output)
}

.densityMapData <- function(object, assay, assayName) {
    expression <- .nativeExpression(object, assay, assayName)
    if (!nrow(expression)) stop("No features to export.", call. = FALSE)
    annotated <- .featureMetadata(object, assay)
    masses <- .massValues(annotated, rownames(expression))
    coordinates <- .nativeCoordinates(object)[, c("x", "y"), drop = FALSE]
    locations <- lapply(seq_len(nrow(expression)), function(i) {
        values <- as.numeric(expression[i, ])
        if (!.allFiniteAssay(values) || any(values < 0)) {
            stop(
                "Density export requires finite non-negative intensities.",
                call. = FALSE
            )
        }
        unname(as.matrix(coordinates[values > 0, , drop = FALSE]))
    })
    list(
        features = data.frame(
            mz = masses, feature = rownames(expression),
            mean = as.numeric(Matrix::rowMeans(expression)),
            annotation = annotated$all_IsomerNames %||%
                annotated$annotation %||% rep("", nrow(expression)),
            identifiers = annotated$all_Isomers_IDs %||% rep("", nrow(expression))
        ),
        locations = locations, xlim = range(coordinates$x), ylim = range(coordinates$y)
    )
}

.densityMapSpectrum <- function(payload) {
    features <- payload$features
    spectrum <- plotly::plot_ly(features,
        x = ~mz, y = ~mean, type = "bar",
        customdata = seq_len(nrow(features)) - 1L, name = "Mean intensity",
        height = 350
    )
    if (sum(features$mean) > 0 && nrow(features) > 1L) {
        kde <- stats::density(features$mz,
            weights = features$mean / sum(features$mean),
            bw = 0.05
        )
        kde$y <- kde$y * max(features$mean) / max(kde$y)
        spectrum <- plotly::add_lines(spectrum,
            x = kde$x, y = kde$y,
            name = "Weighted mass density", inherit = FALSE, hoverinfo = "skip"
        )
    }
    plotly::layout(spectrum,
        xaxis = list(title = "m/z"),
        yaxis = list(title = "Mean intensity"), showlegend = FALSE
    )
}
